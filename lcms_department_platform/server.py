#!/usr/bin/env python3
"""Department LC-MS upload, queue, and report portal.

This is a lightweight single-machine MVP built with the Python standard
library. It accepts RAW/mzML uploads, processes one Peak-first comparison task
at a time, and serves the existing Peak-first report UI from each task output.
"""

from __future__ import annotations

import argparse
import base64
import calendar
import gzip
import html
import io
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import threading
import tempfile
import time
import traceback
import uuid
import zipfile
from dataclasses import dataclass, field
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, quote, unquote, urlparse


WORKSPACE = Path(__file__).resolve().parents[1]
LCMS_FEATURE_DIR = WORKSPACE / "lcms_feature_mvp"
MAX_UPLOAD_BYTES = 8 * 1024 * 1024 * 1024
MAX_FASTA_BYTES = 2 * 1024 * 1024
MAX_STRUCTURE_BYTES = 40 * 1024 * 1024
SPECTRUM_SUFFIXES = {".raw", ".mzml"}
FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".faa", ".txt"}
STRUCTURE_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif"}
DEFAULT_TOP_N_TIC_PEAKS = 80
DEFAULT_TOP_N_MZ = 40
DEFAULT_TOP_N_CHANGED_MZ = 15
DEFAULT_MAX_SPECTRUM_POINTS_PER_SCAN = 500
DEFAULT_MAX_PEAKS_PER_SCAN = 500
if str(LCMS_FEATURE_DIR) not in sys.path:
    sys.path.insert(0, str(LCMS_FEATURE_DIR))

from serve_peak_first_compare import (  # noqa: E402
    VENDOR_DIR,
    pdb_structure_text,
    read_artifact,
    read_bootstrap,
    read_features,
    read_msms_identifications,
    replace_features,
    xic_payload,
)
from core.lcms_msms import read_fasta  # noqa: E402
from run_peak_first_compare import PEAK_FIRST_TEMPLATE  # noqa: E402
try:  # Support both ``python server.py`` and package-level test imports.
    from .mcp_server import MCPApplication  # type: ignore[import-not-found]  # noqa: E402
except ImportError:  # pragma: no cover - exercised by direct script execution
    from mcp_server import MCPApplication  # noqa: E402


def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


class TaskCanceled(Exception):
    """Internal control-flow exception for a user-requested cancellation."""


def safe_id(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value.strip())
    return cleaned.strip("_") or f"task_{int(time.time())}"


def safe_filename(value: str) -> str:
    name = Path(value.replace("\\", "/")).name
    cleaned = "".join(ch if ch.isalnum() or ch in "._-()[] " else "_" for ch in name).strip()
    return cleaned or f"upload_{int(time.time())}.raw"


def sample_name_value(value: str, fallback: str) -> str:
    """Normalize a user-facing sample label while keeping it readable."""
    cleaned = " ".join(str(value or "").split())
    return cleaned[:160] or fallback


def unique_path(directory: Path, filename: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / safe_filename(filename)
    if not path.exists():
        return path
    stem, suffix = path.stem, path.suffix
    index = 1
    while True:
        candidate = directory / f"{stem}_{index}{suffix}"
        if not candidate.exists():
            return candidate
        index += 1


def validate_fasta_file(path: Path) -> dict[str, int]:
    """Validate a user-supplied FASTA reference file."""
    size = path.stat().st_size
    if size <= 0:
        raise ValueError("FASTA file is empty")
    if size > MAX_FASTA_BYTES:
        raise ValueError("FASTA file exceeds the 2 MB limit")
    chains = read_fasta(path)
    if not chains:
        raise ValueError("FASTA contains no usable protein sequence")
    return {str(name): len(sequence) for name, sequence in chains.items()}


def validate_structure_file(path: Path) -> None:
    size = path.stat().st_size
    if size <= 0:
        raise ValueError("structure file is empty")
    if size > MAX_STRUCTURE_BYTES:
        raise ValueError("structure file exceeds the 40 MB limit")
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        raise ValueError("structure file is empty")
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        looks_like_structure = "data_" in text[:1000] or "_atom_site." in text
    else:
        looks_like_structure = any(token in text[:20000] for token in ("ATOM", "HETATM", "HEADER", "MODEL"))
    if not looks_like_structure:
        raise ValueError("file does not look like a PDB or mmCIF structure")


def structure_payload_from_file(path: Path, label: str, source_meta: dict[str, object] | None = None) -> dict[str, object]:
    suffix = path.suffix.lower()
    return {
        "available": True,
        "format": "cif" if suffix in {".cif", ".mmcif"} else "pdb",
        "text": path.read_text(encoding="utf-8", errors="replace"),
        "label": label,
        "source_url": "",
        "source_meta": source_meta or {"kind": "uploaded"},
    }


def read_json_file(path: Path, default: object) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else default
    except (OSError, json.JSONDecodeError):
        return default


def upsert_sqlite_artifact(db_path: Path, key: str, payload: object) -> None:
    with sqlite3.connect(db_path) as connection:
        connection.execute(
            """
            INSERT INTO peak_first_artifacts (artifact_key, payload_json, updated_at)
            VALUES (?, ?, CURRENT_TIMESTAMP)
            ON CONFLICT(artifact_key) DO UPDATE SET
                payload_json = excluded.payload_json,
                updated_at = excluded.updated_at
            """,
            (key, json.dumps(payload, ensure_ascii=False, separators=(",", ":"))),
        )
        connection.commit()


def update_bootstrap_metadata(db_path: Path, metadata: dict[str, object]) -> None:
    """Merge task-level execution metadata into the report bootstrap artifact."""
    bootstrap = read_artifact(db_path, "bootstrap")
    if not isinstance(bootstrap, dict):
        raise ValueError("bootstrap artifact is not a JSON object")
    merged = dict(bootstrap.get("analysis_metadata") or {})
    merged.update(metadata)
    bootstrap["analysis_metadata"] = merged
    upsert_sqlite_artifact(db_path, "bootstrap", bootstrap)


def parse_content_disposition(value: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for item in value.split(";"):
        item = item.strip()
        if "=" not in item:
            continue
        key, raw = item.split("=", 1)
        fields[key.strip().lower()] = raw.strip().strip('"')
    return fields


def parse_multipart(content_type: str, body: bytes) -> tuple[dict[str, str], list[tuple[str, bytes]]]:
    marker = "boundary="
    if marker not in content_type:
        raise ValueError("request must be multipart/form-data")
    boundary = content_type.split(marker, 1)[1].split(";", 1)[0].strip().strip('"')
    if not boundary:
        raise ValueError("multipart boundary is missing")
    delimiter = b"--" + boundary.encode("utf-8")
    fields: dict[str, str] = {}
    files: list[tuple[str, bytes]] = []
    for part in body.split(delimiter):
        part = part.strip(b"\r\n")
        if not part or part == b"--":
            continue
        if part.endswith(b"--"):
            part = part[:-2].rstrip(b"\r\n")
        if b"\r\n\r\n" not in part:
            continue
        raw_headers, payload = part.split(b"\r\n\r\n", 1)
        headers: dict[str, str] = {}
        for line in raw_headers.decode("utf-8", errors="replace").split("\r\n"):
            if ":" in line:
                key, value = line.split(":", 1)
                headers[key.strip().lower()] = value.strip()
        disposition = parse_content_disposition(headers.get("content-disposition", ""))
        name = disposition.get("name", "")
        filename = disposition.get("filename", "")
        payload = payload.rstrip(b"\r\n")
        if filename:
            files.append((safe_filename(filename), payload))
        elif name:
            fields[name] = payload.decode("utf-8", errors="replace").strip()
    return fields, files


class MultipartReader:
    """Small pushback reader used to stream multipart file bodies to disk."""

    def __init__(self, stream: Any):
        self.stream = stream
        self.pending = bytearray()

    def push(self, data: bytes) -> None:
        if data:
            self.pending = bytearray(data) + self.pending

    def read(self, size: int = 1024 * 1024) -> bytes:
        if self.pending:
            chunk = bytes(self.pending[:size])
            del self.pending[:size]
            return chunk
        # BufferedReader.read(size) may wait for the full requested size on a
        # request body. read1() returns whatever has arrived, which prevents
        # a small upload from waiting for the browser's response round-trip.
        read1 = getattr(self.stream, "read1", None)
        return read1(size) if read1 else self.stream.read(size)

    def readline(self) -> bytes:
        line = bytearray()
        while True:
            index = self.pending.find(b"\n")
            if index >= 0:
                line.extend(self.pending[: index + 1])
                del self.pending[: index + 1]
                return bytes(line)
            if self.pending:
                line.extend(self.pending)
                self.pending.clear()
            chunk = self.stream.readline()
            if not chunk:
                return bytes(line)
            line.extend(chunk)
            if chunk.endswith(b"\n"):
                return bytes(line)


def parse_multipart_stream(content_type: str, stream: Any, staging_dir: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    """Parse form fields while streaming uploaded files into ``staging_dir``.

    Only a small rolling buffer is retained while a file is being written, so
    large RAW uploads do not occupy an equally large Python byte string.
    """
    marker = "boundary="
    if marker not in content_type:
        raise ValueError("request must be multipart/form-data")
    boundary = content_type.split(marker, 1)[1].split(";", 1)[0].strip().strip('"')
    if not boundary:
        raise ValueError("multipart boundary is missing")
    delimiter = b"--" + boundary.encode("utf-8")
    boundary_marker = b"\r\n" + delimiter
    reader = MultipartReader(stream)
    fields: dict[str, str] = {}
    files: list[dict[str, str]] = []
    staging_dir.mkdir(parents=True, exist_ok=True)

    def consume_payload(target: Any = None) -> bytes:
        buffer = bytearray()
        field_payload = io.BytesIO() if target is None else None
        keep = len(boundary_marker) + 4
        while True:
            chunk = reader.read()
            if not chunk:
                raise ValueError("incomplete multipart upload")
            buffer.extend(chunk)
            position = buffer.find(boundary_marker)
            if position >= 0:
                payload = bytes(buffer[:position])
                if target is not None:
                    target.write(payload)
                else:
                    field_payload.write(payload)
                    if field_payload.tell() > 2 * 1024 * 1024:
                        raise ValueError("multipart form field is too large")
                reader.push(bytes(buffer[position + 2 :]))
                return field_payload.getvalue() if target is None else b""
            if len(buffer) > keep:
                output = bytes(buffer[:-keep])
                if target is not None:
                    target.write(output)
                else:
                    field_payload.write(output)
                    if field_payload.tell() > 2 * 1024 * 1024:
                        raise ValueError("multipart form field is too large")
                del buffer[:-keep]

    first = reader.readline().strip()
    if first != delimiter:
        raise ValueError("invalid multipart boundary")
    while True:
        raw_headers: list[str] = []
        while True:
            line = reader.readline()
            if not line:
                raise ValueError("incomplete multipart headers")
            stripped = line.rstrip(b"\r\n")
            if not stripped:
                break
            raw_headers.append(stripped.decode("utf-8", errors="replace"))
        headers: dict[str, str] = {}
        for line in raw_headers:
            if ":" in line:
                key, value = line.split(":", 1)
                headers[key.strip().lower()] = value.strip()
        disposition = parse_content_disposition(headers.get("content-disposition", ""))
        name = disposition.get("name", "")
        filename = disposition.get("filename", "")
        if filename:
            target = unique_path(staging_dir, filename)
            with target.open("wb") as handle:
                consume_payload(handle)
            files.append({
                "field_name": name,
                "original_name": safe_filename(filename),
                "stored_name": target.name,
                "staged_path": str(target),
            })
        else:
            payload = consume_payload()
            if name:
                fields[name] = payload.decode("utf-8", errors="replace").strip()
        boundary_line = reader.readline().strip()
        if boundary_line == delimiter + b"--":
            break
        if boundary_line != delimiter:
            raise ValueError("invalid multipart boundary terminator")
    return fields, files


@dataclass
class PortalConfig:
    host: str
    port: int
    root: Path
    parser_path: Path
    python_executable: str
    jobs_dir_override: Path | None = None
    max_upload_bytes: int = MAX_UPLOAD_BYTES

    @property
    def state_dir(self) -> Path:
        return self.root / "state"

    @property
    def jobs_dir(self) -> Path:
        return (self.jobs_dir_override or (self.root / "jobs")).resolve()

    @property
    def tasks_path(self) -> Path:
        return self.state_dir / "tasks.json"


@dataclass
class TaskStore:
    config: PortalConfig
    lock: threading.Lock = field(default_factory=threading.Lock)

    def _load_unlocked(self) -> list[dict[str, Any]]:
        if not self.config.tasks_path.exists():
            return []
        try:
            payload = json.loads(self.config.tasks_path.read_text(encoding="utf-8"))
            return payload if isinstance(payload, list) else []
        except json.JSONDecodeError:
            return []

    def _save_unlocked(self, tasks: list[dict[str, Any]]) -> None:
        self.config.state_dir.mkdir(parents=True, exist_ok=True)
        tmp = self.config.tasks_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(tasks, ensure_ascii=False, indent=2), encoding="utf-8")
        tmp.replace(self.config.tasks_path)

    def load(self) -> list[dict[str, Any]]:
        with self.lock:
            return self._load_unlocked()

    def save(self, tasks: list[dict[str, Any]]) -> None:
        with self.lock:
            self._save_unlocked(tasks)

    def update(self, task_id: str, **updates: Any) -> dict[str, Any]:
        with self.lock:
            tasks = self._load_unlocked()
            for task in tasks:
                if task.get("task_id") == task_id:
                    task.update(updates)
                    task["updated_at"] = now_iso()
                    self._save_unlocked(tasks)
                    return task
            raise KeyError(task_id)

    def get(self, task_id: str) -> dict[str, Any] | None:
        for task in self.load():
            if task.get("task_id") == task_id:
                return task
        return None

    def create(self, task: dict[str, Any]) -> dict[str, Any]:
        with self.lock:
            tasks = self._load_unlocked()
            tasks.insert(0, task)
            self._save_unlocked(tasks)
            return task

    def recover_interrupted(self) -> int:
        """Recover tasks left behind when the previous server stopped."""
        recovered = 0
        with self.lock:
            tasks = self._load_unlocked()
            for task in tasks:
                if task.get("status") == "canceling" or task.get("cancel_requested"):
                    task.update(
                        status="canceled",
                        stage="canceled",
                        progress=int(task.get("progress") or 0),
                        cancel_requested=False,
                        finished_at=now_iso(),
                        error="Task canceled while the previous server process was stopping.",
                    )
                    task["updated_at"] = now_iso()
                    recovered += 1
                elif task.get("status") == "running":
                    task.update(
                        status="waiting",
                        stage="requeued after restart",
                        progress=0,
                        started_at="",
                        finished_at="",
                        error="Previous server stopped before completion; task was requeued automatically.",
                    )
                    task["updated_at"] = now_iso()
                    recovered += 1
            if recovered:
                self._save_unlocked(tasks)
        return recovered

    def claim_waiting(self, task_id: str) -> dict[str, Any] | None:
        """Atomically move one waiting task to running for the worker."""
        with self.lock:
            tasks = self._load_unlocked()
            for task in tasks:
                if task.get("task_id") == task_id:
                    if task.get("status") != "waiting":
                        return None
                    task.update(
                        status="running",
                        stage="preparing",
                        progress=5,
                        started_at=now_iso(),
                        error="",
                        cancel_requested=False,
                    )
                    task["updated_at"] = now_iso()
                    self._save_unlocked(tasks)
                    return task
            return None

    def request_cancel(self, task_id: str) -> dict[str, Any]:
        """Atomically request cancellation for a waiting or running task."""
        with self.lock:
            tasks = self._load_unlocked()
            for task in tasks:
                if task.get("task_id") != task_id:
                    continue
                status = str(task.get("status") or "")
                if status == "waiting":
                    task.update(
                        status="canceled",
                        stage="canceled",
                        progress=0,
                        cancel_requested=False,
                        finished_at=now_iso(),
                        error="Canceled by user.",
                    )
                elif status == "running":
                    task.update(
                        status="canceling",
                        stage="canceling",
                        cancel_requested=True,
                        error="Cancellation requested by user; stopping the current process.",
                    )
                elif status in {"canceling", "canceled"}:
                    return task
                else:
                    raise ValueError("only waiting or running tasks can be canceled")
                task["updated_at"] = now_iso()
                self._save_unlocked(tasks)
                return task
            raise KeyError(task_id)

    def mark_canceled(self, task_id: str, reason: str) -> dict[str, Any] | None:
        """Persist the terminal canceled state unless a task already finished."""
        with self.lock:
            tasks = self._load_unlocked()
            for task in tasks:
                if task.get("task_id") != task_id:
                    continue
                if task.get("status") == "finished":
                    return task
                task.update(
                    status="canceled",
                    stage="canceled",
                    cancel_requested=False,
                    finished_at=task.get("finished_at") or now_iso(),
                    error=reason,
                )
                task["updated_at"] = now_iso()
                self._save_unlocked(tasks)
                return task
            return None

    def finish_if_not_canceled(self, task_id: str, **updates: Any) -> dict[str, Any] | None:
        """Finish only if cancellation did not win the state transition race."""
        with self.lock:
            tasks = self._load_unlocked()
            for task in tasks:
                if task.get("task_id") != task_id:
                    continue
                if task.get("status") != "running" or task.get("cancel_requested"):
                    return None
                task.update(updates)
                task["updated_at"] = now_iso()
                self._save_unlocked(tasks)
                return task
            return None


def iso_seconds(started_at: object, finished_at: object = "") -> float | None:
    values = [str(started_at or ""), str(finished_at or "")]
    if not values[0]:
        return None
    try:
        start = calendar.timegm(time.strptime(values[0], "%Y-%m-%dT%H:%M:%SZ"))
        end = calendar.timegm(time.strptime(values[1], "%Y-%m-%dT%H:%M:%SZ")) if values[1] else time.time()
        return max(0.0, end - start)
    except (ValueError, OverflowError):
        return None


def task_progress(task: dict[str, Any]) -> int:
    if task.get("status") == "finished":
        return 100
    if task.get("status") in {"failed", "canceled"}:
        return int(task.get("progress") or 0)
    explicit = task.get("progress")
    if explicit is not None:
        try:
            return max(0, min(99, int(explicit)))
        except (TypeError, ValueError):
            pass
    stage = str(task.get("stage") or "")
    if stage == "waiting":
        return 0
    if "converting" in stage or "processing" in stage:
        return 20
    if "building" in stage:
        return 75
    return 5


def task_snapshots(tasks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    waiting = [task for task in tasks if task.get("status") == "waiting"]
    waiting.sort(key=lambda item: str(item.get("created_at") or ""))
    positions = {str(task.get("task_id")): index + 1 for index, task in enumerate(waiting)}
    snapshots: list[dict[str, Any]] = []
    for task in tasks:
        snapshot = dict(task)
        snapshot.pop("job_dir", None)
        snapshot.pop("output_dir", None)
        snapshot["progress"] = task_progress(task)
        snapshot["queue_position"] = positions.get(str(task.get("task_id")))
        snapshot["duration_seconds"] = iso_seconds(task.get("started_at"), task.get("finished_at"))
        snapshot["total_size_bytes"] = sum(int(item.get("size_bytes") or 0) for item in task.get("files") or []) + sum(
            int(item.get("size_bytes") or 0) for item in (task.get("reference_files") or [])
        )
        snapshot["can_retry"] = task.get("status") == "failed"
        snapshot["can_cancel"] = task.get("status") in {"waiting", "running"}
        snapshot["has_result"] = task.get("status") == "finished" and bool(task.get("report_url"))
        snapshots.append(snapshot)
    return snapshots


def task_directory(config: PortalConfig, task_id: str) -> Path:
    """Return a task directory only when it stays inside the configured jobs root."""
    root = config.jobs_dir.resolve()
    candidate = (root / task_id).resolve()
    if not candidate.is_relative_to(root) or candidate == root:
        raise ValueError("invalid task id")
    return candidate


def clear_task_outputs(config: PortalConfig, task_id: str) -> None:
    """Remove only generated retry artifacts for one known task."""
    task_dir = task_directory(config, task_id)
    # Keep converted mzML because RAW upload copies may be removed after the
    # first successful conversion. Retries can reuse the converted inputs.
    for child_name in ("output", "exports"):
        child = (task_dir / child_name).resolve()
        if child.exists() and child.is_relative_to(task_dir) and child != task_dir:
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()


class TaskWorker(threading.Thread):
    daemon = True

    def __init__(self, store: TaskStore):
        super().__init__(name="lcms-task-worker")
        self.store = store
        self.stop_requested = threading.Event()

    def run(self) -> None:
        while not self.stop_requested.is_set():
            task = self.next_waiting_task()
            if task is None:
                self.stop_requested.wait(2.0)
                continue
            self.process_task(task)

    def next_waiting_task(self) -> dict[str, Any] | None:
        tasks = self.store.load()
        waiting = [task for task in tasks if task.get("status") == "waiting"]
        waiting.sort(key=lambda item: str(item.get("created_at") or ""))
        return waiting[0] if waiting else None

    def log(self, task_id: str, message: str) -> None:
        task_dir = self.store.config.jobs_dir / task_id
        task_dir.mkdir(parents=True, exist_ok=True)
        with (task_dir / "worker.log").open("a", encoding="utf-8") as handle:
            handle.write(f"[{now_iso()}] {message}\n")

    def cancellation_requested(self, task_id: str) -> bool:
        task = self.store.get(task_id)
        return bool(
            task
            and (
                task.get("cancel_requested")
                or task.get("status") in {"canceling", "canceled"}
            )
        )

    def check_cancellation(self, task_id: str) -> None:
        if self.cancellation_requested(task_id):
            raise TaskCanceled("Task cancellation requested by user.")

    def terminate_process_tree(self, process: subprocess.Popen[str]) -> None:
        """Stop a converter/search process and its descendants on cancellation."""
        if process.poll() is not None:
            return
        try:
            process.terminate()
        except OSError:
            pass
        try:
            process.wait(timeout=3)
            return
        except subprocess.TimeoutExpired:
            pass
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
                timeout=10,
            )
        else:
            try:
                process.kill()
            except OSError:
                pass
        try:
            process.wait(timeout=3)
        except subprocess.TimeoutExpired:
            pass

    def run_command(self, task_id: str, command: list[str], cwd: Path) -> None:
        self.check_cancellation(task_id)
        self.log(task_id, "RUN " + " ".join(command))
        creation_flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
        process = subprocess.Popen(
            command,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            creationflags=creation_flags,
        )
        while True:
            try:
                stdout, stderr = process.communicate(timeout=0.5)
                break
            except subprocess.TimeoutExpired:
                if self.cancellation_requested(task_id):
                    self.log(task_id, f"CANCEL requested; terminating process tree PID {process.pid}.")
                    self.terminate_process_tree(process)
                    stdout, stderr = process.communicate()
                    if stdout.strip():
                        self.log(task_id, stdout.strip())
                    if stderr.strip():
                        self.log(task_id, stderr.strip())
                    raise TaskCanceled("Task canceled while the process was running.")
        if stdout.strip():
            self.log(task_id, stdout.strip())
        if stderr.strip():
            self.log(task_id, stderr.strip())
        if self.cancellation_requested(task_id):
            raise TaskCanceled("Task canceled after the process completed.")
        if process.returncode != 0:
            raise RuntimeError(f"Command failed with exit code {process.returncode}: {' '.join(command)}")

    def convert_raw(self, task_id: str, raw_path: Path, mzml_dir: Path) -> Path:
        parser = self.store.config.parser_path
        if not parser.exists():
            raise FileNotFoundError(f"ThermoRawFileParser not found: {parser}")
        target = mzml_dir / f"{raw_path.stem}.mzML"
        if target.exists():
            return target
        command = [str(parser), "-i", str(raw_path), "-o", str(mzml_dir), "-f", "2", "-m", "0", "-l", "3"]
        self.run_command(task_id, command, WORKSPACE)
        if target.exists():
            return target
        matches = sorted(mzml_dir.glob(f"{raw_path.stem}*.mzML"), key=lambda item: item.stat().st_mtime, reverse=True)
        if matches:
            return matches[0]
        raise FileNotFoundError(f"Converted mzML was not created for {raw_path.name}")

    def cleanup_task_raw_copies(
        self,
        task_id: str,
        task_dir: Path,
        raw_dir: Path,
        input_files: list[dict[str, Any]],
    ) -> None:
        """Delete only uploaded RAW copies after every conversion succeeded."""
        task_root = task_directory(self.store.config, task_id)
        safe_raw_dir = raw_dir.resolve()
        if not safe_raw_dir.is_relative_to(task_root) or safe_raw_dir == task_root:
            raise ValueError("refusing to clean RAW files outside the task directory")
        deleted = 0
        for file_info in input_files:
            stored_name = str(file_info.get("stored_name") or "")
            if Path(stored_name).suffix.lower() != ".raw":
                continue
            candidate = (safe_raw_dir / stored_name).resolve()
            if not candidate.is_relative_to(safe_raw_dir) or not candidate.is_relative_to(task_root):
                raise ValueError("refusing to clean an unsafe RAW task path")
            if not candidate.is_file():
                continue
            try:
                candidate.unlink()
                deleted += 1
            except OSError as error:
                self.log(task_id, f"WARNING: could not remove task RAW copy {candidate.name}: {error}")
        if deleted:
            self.log(task_id, f"Removed {deleted} uploaded RAW copy/copies after successful mzML conversion; source RAW files were not touched.")

    def process_task(self, task: dict[str, Any]) -> None:
        task_id = str(task["task_id"])
        claimed = self.store.claim_waiting(task_id)
        if claimed is None:
            return
        task = claimed
        task_dir = self.store.config.jobs_dir / task_id
        raw_dir = task_dir / "raw"
        mzml_dir = task_dir / "mzML"
        output_dir = task_dir / "output"
        references_dir = task_dir / "references"
        try:
            self.log(task_id, "Task started.")
            mzml_dir.mkdir(parents=True, exist_ok=True)
            converted: list[str] = []
            input_files = task.get("files", [])
            total_inputs = max(1, len(input_files))
            for index, file_info in enumerate(input_files):
                source = raw_dir / str(file_info["stored_name"])
                suffix = source.suffix.lower()
                file_progress = 10 + int(35 * index / total_inputs)
                self.store.update(task_id, stage=f"processing {source.name}", progress=file_progress)
                if suffix == ".mzml":
                    target = mzml_dir / source.name
                    if not target.exists():
                        shutil.copy2(source, target)
                    converted.append(str(target))
                elif suffix == ".raw":
                    cached_mzml = mzml_dir / f"{source.stem}.mzML"
                    if not source.exists() and cached_mzml.exists():
                        self.log(task_id, f"Reusing existing converted mzML for removed RAW copy: {cached_mzml.name}")
                        converted.append(str(cached_mzml))
                    else:
                        converted.append(str(self.convert_raw(task_id, source, mzml_dir)))
                else:
                    raise ValueError(f"Unsupported file type: {source.name}")
            if len(converted) < 2:
                raise ValueError("At least two RAW/mzML files are required for comparison.")
            self.cleanup_task_raw_copies(task_id, task_dir, raw_dir, input_files)
            self.store.update(task_id, stage="building peak-first comparison", progress=55)
            output_dir.mkdir(parents=True, exist_ok=True)
            params = task.get("params") or {}
            sample_names_path = task_dir / "sample_names.json"
            sample_names_path.write_text(
                json.dumps(
                    {
                        Path(str(file_info["stored_name"])).stem: str(file_info.get("sample_name") or Path(str(file_info["stored_name"])).stem)
                        for file_info in input_files
                    },
                    ensure_ascii=False,
                    indent=2,
                ),
                encoding="utf-8",
            )
            command = [
                self.store.config.python_executable,
                str(WORKSPACE / "lcms_feature_mvp" / "run_peak_first_compare.py"),
                "--input-dir",
                str(mzml_dir),
                "--output-dir",
                str(output_dir),
                "--project-id",
                task_id,
                "--top-n-peaks",
                str(params.get("top_n_peaks", DEFAULT_TOP_N_TIC_PEAKS)),
                "--top-n-mz",
                str(params.get("top_n_mz", DEFAULT_TOP_N_MZ)),
                "--top-n-changed-mz",
                str(params.get("top_n_changed_mz", DEFAULT_TOP_N_CHANGED_MZ)),
                "--mz-tolerance-mode",
                "da",
                "--mz-tolerance-da",
                str(params.get("mz_tolerance_da", 0.16)),
                "--mz-tolerance-ppm",
                str(params.get("mz_tolerance_ppm", 10.0)),
                "--max-spectrum-points-per-scan",
                str(params.get("max_spectrum_points_per_scan", DEFAULT_MAX_SPECTRUM_POINTS_PER_SCAN)),
                "--max-peaks-per-scan",
                str(params.get("max_peaks_per_scan", DEFAULT_MAX_PEAKS_PER_SCAN)),
                "--spectrum-min-intensity",
                "0",
                "--sample-names-json",
                str(sample_names_path),
            ]
            self.run_command(task_id, command, WORKSPACE)
            report = output_dir / "lcms_peak_first_compare.html"
            db = output_dir / "lcms_peak_first_compare.sqlite"
            if not report.exists() or not db.exists():
                raise FileNotFoundError("Peak-first report or SQLite output was not created.")
            reference_meta = dict(task.get("references") or {})
            sequence_meta = dict(reference_meta.get("sequence") or {})
            structure_meta = dict(reference_meta.get("structure") or {})
            analysis_metadata: dict[str, object] = {
                "method": {
                    "centroid_match_tolerance_ppm": 10.0,
                    "xic_half_window_da": 0.16,
                    "isotope_spacing_da": 1.00335483507,
                    "analysis_limits": {
                        "top_tic_peaks": int(params.get("top_n_peaks", DEFAULT_TOP_N_TIC_PEAKS)),
                        "top_mz_per_peak": int(params.get("top_n_mz", DEFAULT_TOP_N_MZ)),
                        "top_changed_mz": int(params.get("top_n_changed_mz", DEFAULT_TOP_N_CHANGED_MZ)),
                        "max_spectrum_points_per_scan": int(params.get("max_spectrum_points_per_scan", DEFAULT_MAX_SPECTRUM_POINTS_PER_SCAN)),
                        "max_peaks_per_scan": int(params.get("max_peaks_per_scan", DEFAULT_MAX_PEAKS_PER_SCAN)),
                    },
                },
                "ms1": {"status": "completed"},
                "ms2": {"status": "skipped_no_fasta", "reason": "未提供 FASTA 序列，未执行后续 MS2 计算。"},
                "structure_mapping": {"status": "skipped_no_fasta", "reason": "未提供 FASTA 序列，未执行蛋白质结构映射。"},
                "references": reference_meta,
            }
            fasta_path = Path(str(sequence_meta.get("path") or "")) if sequence_meta.get("path") else None
            if fasta_path and not fasta_path.is_absolute():
                fasta_path = task_dir / fasta_path
            if fasta_path and fasta_path.exists():
                self.store.update(task_id, stage="building MS/MS identification", progress=78)
                msms_command = [
                    self.store.config.python_executable,
                    str(WORKSPACE / "lcms_feature_mvp" / "run_msms_compare.py"),
                    "--fasta",
                    str(fasta_path),
                    "--output-dir",
                    str(output_dir),
                    "--feature-sqlite",
                    str(db),
                    "--project-id",
                    task_id,
                ]
                for mzml_path in converted:
                    msms_command.extend(["--mzml", mzml_path])
                self.run_command(task_id, msms_command, WORKSPACE)
                msms_report = output_dir / "lcms_msms_report.html"
                if not msms_report.exists():
                    raise FileNotFoundError("FASTA was provided but the MS/MS report was not created.")
                analysis_metadata["ms2"] = {
                    "status": "completed",
                    "report_available": True,
                    "fasta": sequence_meta.get("name", fasta_path.name),
                }
                if structure_meta.get("provided"):
                    structure_ready = bool(structure_meta.get("path"))
                    if structure_meta.get("source") == "pdb_id" and not structure_meta.get("path"):
                        try:
                            structure_text, _ = pdb_structure_text(str(structure_meta.get("pdb_id")))
                            structure_target = references_dir / "structure.cif"
                            structure_target.write_text(structure_text, encoding="utf-8")
                            structure_meta["path"] = "references/structure.cif"
                            reference_meta["structure"] = structure_meta
                            task["references"] = reference_meta
                            self.store.update(task_id, references=reference_meta)
                            structure_ready = True
                        except Exception as structure_error:
                            self.log(task_id, f"PDB structure fetch deferred/failed: {structure_error}")
                    analysis_metadata["structure_mapping"] = {
                        "status": "available" if structure_ready else "skipped_structure_unavailable",
                        "reason": (
                            "已提供 FASTA 与结构文件/PDB 结构；报告中可进行序列到结构的映射。"
                            if structure_ready
                            else "PDB 结构未能下载，未执行蛋白质结构映射。"
                        ),
                    }
                else:
                    analysis_metadata["structure_mapping"] = {
                        "status": "skipped_no_structure",
                        "reason": "已执行 MS2，但未提供结构文件或 PDB ID，未执行结构映射。",
                    }
            elif structure_meta.get("provided"):
                analysis_metadata["structure_mapping"] = {
                    "status": "skipped_no_fasta",
                    "reason": "已提供结构，但缺少 FASTA 序列；按规则不执行结构映射。",
                }
            update_bootstrap_metadata(db, analysis_metadata)
            if structure_meta.get("provided") and structure_meta.get("path"):
                structure_path = (task_dir / str(structure_meta["path"])).resolve()
                task_root = task_dir.resolve()
                if structure_path.is_relative_to(task_root) and structure_path.exists():
                    structure_payload = {
                        "available": True,
                        "format": "cif" if structure_path.suffix.lower() in {".cif", ".mmcif"} else "pdb",
                        "text": structure_path.read_text(encoding="utf-8", errors="replace"),
                        "label": structure_meta.get("name") or "任务结构",
                        "source_url": (
                            f"https://www.rcsb.org/structure/{quote(str(structure_meta.get('pdb_id')))}"
                            if structure_meta.get("pdb_id") else ""
                        ),
                        "source_meta": {"kind": "pdb" if structure_meta.get("pdb_id") else "uploaded"},
                    }
                    upsert_sqlite_artifact(db, "task_reference_structure", structure_payload)
            self.check_cancellation(task_id)
            if self.store.finish_if_not_canceled(
                task_id,
                status="finished",
                stage="finished",
                progress=100,
                finished_at=now_iso(),
                report_url=f"/tasks/{quote(task_id)}/report/",
                output_dir=str(output_dir),
                cancel_requested=False,
            ) is None:
                raise TaskCanceled("Task canceled before finalizing the result.")
            self.log(task_id, "Task finished.")
        except TaskCanceled as exc:
            self.log(task_id, "CANCELED " + str(exc))
            self.store.mark_canceled(task_id, str(exc))
        except Exception as exc:  # pragma: no cover - operational path
            if self.cancellation_requested(task_id):
                self.log(task_id, "CANCELED " + str(exc))
                self.store.mark_canceled(task_id, "Task canceled while processing.")
                return
            self.log(task_id, "FAILED " + str(exc))
            self.log(task_id, traceback.format_exc())
            self.store.update(task_id, status="failed", stage="failed", progress=0, error=str(exc), finished_at=now_iso())


def page_shell(title: str, body: str) -> bytes:
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      --ink:#17253d;
      --muted:#718096;
      --line:#dce6f2;
      --canvas:#f3f7fb;
      --card:#ffffff;
      --navy:#0a1730;
      --blue:#2563eb;
      --cyan:#0ea5b7;
      --green:#15966b;
      --amber:#c77710;
      --red:#c2414d;
      --shadow:0 18px 45px rgba(31,58,92,.08);
    }}
    * {{ box-sizing:border-box; }}
    html {{ min-width:320px; background:var(--canvas); }}
    body {{ margin:0; min-height:100vh; font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI","Microsoft YaHei",sans-serif; color:var(--ink); background:var(--canvas); line-height:1.5; overflow-x:hidden; }}
    body::before {{ content:""; position:fixed; z-index:-1; width:580px; height:580px; right:-220px; top:80px; border-radius:50%; background:radial-gradient(circle,rgba(38,99,235,.08),rgba(38,99,235,0) 68%); pointer-events:none; }}
    header {{ position:sticky; top:0; z-index:20; min-height:72px; padding:14px clamp(18px,4vw,56px); display:flex; gap:24px; align-items:center; justify-content:space-between; color:#fff; background:rgba(10,23,48,.96); border-bottom:1px solid rgba(148,191,255,.18); box-shadow:0 8px 24px rgba(10,23,48,.14); backdrop-filter:blur(14px); }}
    .header-brand {{ display:flex; min-width:0; align-items:center; gap:12px; }}
    .brand-mark {{ position:relative; display:grid; width:38px; height:38px; flex:0 0 auto; place-items:center; color:#dffbff; font-size:18px; font-weight:800; border:1px solid rgba(114,230,240,.62); border-radius:12px; background:linear-gradient(145deg,#1c6ee8,#0ea5b7); box-shadow:0 6px 18px rgba(14,165,183,.26); }}
    .brand-mark::after {{ content:""; position:absolute; right:7px; bottom:7px; width:6px; height:6px; border:2px solid #dffbff; border-radius:50%; }}
    .brand-kicker {{ color:#8ca9ce; font-size:9px; font-weight:800; letter-spacing:.18em; line-height:1.2; }}
    header h1 {{ margin:2px 0 0; color:#f6faff; font-size:17px; font-weight:750; letter-spacing:.01em; line-height:1.2; }}
    .top-nav {{ display:flex; min-width:0; align-items:center; gap:14px; }}
    header a {{ color:#d8e8ff; font-size:13px; font-weight:650; text-decoration:none; }}
    header a:hover {{ color:#fff; }}
    .top-nav a {{ padding:8px 12px; border:1px solid rgba(153,192,237,.2); border-radius:9px; background:rgba(255,255,255,.07); }}
    .nav-divider {{ width:1px; height:22px; background:rgba(155,184,222,.24); }}
    .header-note {{ color:#9fb6d5; font-size:11px; white-space:nowrap; }}
    main {{ position:relative; width:min(1480px,100%); margin:0 auto; padding:30px clamp(16px,3.5vw,52px) 58px; display:grid; gap:20px; }}
    section {{ position:relative; min-width:0; padding:24px; overflow:hidden; border:1px solid var(--line); border-radius:18px; background:var(--card); box-shadow:var(--shadow); }}
    section::before {{ content:""; position:absolute; left:0; top:0; width:100%; height:3px; background:linear-gradient(90deg,rgba(37,99,235,.8),rgba(14,165,183,.32),transparent 80%); opacity:.72; }}
    h2 {{ margin:0 0 16px; color:#182945; font-size:16px; font-weight:760; letter-spacing:-.01em; }}
    h2::first-letter {{ color:var(--blue); }}
    label {{ display:block; margin:13px 0 6px; color:#344965; font-size:12px; font-weight:750; letter-spacing:.01em; }}
    input, select, textarea, button {{ border:1px solid #cfdcea; border-radius:9px; padding:10px 12px; color:var(--ink); font:inherit; font-size:13px; background:#fff; outline:none; transition:border-color .18s ease,box-shadow .18s ease,background .18s ease,transform .18s ease; }}
    input::placeholder, textarea::placeholder {{ color:#a0aec0; }}
    input:focus, select:focus, textarea:focus {{ border-color:#6ca5f5; box-shadow:0 0 0 3px rgba(37,99,235,.12); }}
    input[type=file] {{ width:100%; padding:0; border:0; background:transparent; color:#61748e; font-size:12px; }}
    input[type=file]::file-selector-button {{ margin-right:9px; padding:8px 11px; border:1px solid #bfd1ea; border-radius:8px; color:#245291; font:inherit; font-weight:700; background:#eef5ff; cursor:pointer; }}
    input[type=file]::file-selector-button:hover {{ background:#e2efff; }}
    textarea {{ display:block; width:100%; resize:vertical; line-height:1.55; }}
    button {{ cursor:pointer; font-weight:700; }}
    button:hover:not(:disabled) {{ transform:translateY(-1px); }}
    button:disabled {{ cursor:default; opacity:.54; }}
    #submitTask {{ display:inline-flex; align-items:center; gap:8px; padding:11px 17px; color:#fff; border-color:#1d59cf; background:linear-gradient(135deg,#2563eb,#1f77c9); box-shadow:0 8px 18px rgba(37,99,235,.2); }}
    #submitTask::before {{ content:"↗"; font-size:15px; }}
    table {{ width:100%; min-width:0; table-layout:fixed; border-collapse:separate; border-spacing:0; font-size:12px; }}
    th, td {{ padding:13px 12px; text-align:left; vertical-align:top; border-bottom:1px solid #e9eff6; }}
    th {{ color:#7b8da5; font-size:10px; font-weight:800; letter-spacing:.09em; text-transform:uppercase; background:#f8fafd; }}
    th:first-child {{ border-radius:10px 0 0 10px; }}
    th:last-child {{ border-radius:0 10px 10px 0; }}
    tr[data-task] {{ cursor:pointer; transition:background .16s ease; }}
    tr[data-task]:hover {{ background:#f7fbff; }}
    tr:last-child td {{ border-bottom:0; }}
    .grid {{ display:grid; grid-template-columns:minmax(360px,.62fr) minmax(0,1.38fr); gap:20px; align-items:start; }}
    .grid > section {{ min-width:0; }}
    .page-intro {{ display:flex; align-items:flex-end; justify-content:space-between; gap:24px; padding:6px 4px 2px; }}
    .page-intro-copy {{ max-width:760px; }}
    .eyebrow {{ margin:0 0 8px; color:var(--cyan); font-size:10px; font-weight:850; letter-spacing:.18em; text-transform:uppercase; }}
    .page-intro h2 {{ margin:0; color:#152947; font-size:clamp(25px,3vw,36px); letter-spacing:-.045em; line-height:1.12; }}
    .intro-note {{ max-width:360px; margin:0 0 3px; color:#71839c; font-size:12px; line-height:1.65; text-align:right; }}
    .workflow-strip {{ display:grid; grid-template-columns:repeat(4,1fr); gap:0; padding:14px 18px; border:1px solid #dce8f4; border-radius:14px; background:linear-gradient(100deg,#eef6ff,#f6fbfc); box-shadow:0 10px 28px rgba(31,58,92,.04); }}
    .workflow-step {{ position:relative; display:flex; align-items:center; gap:10px; min-width:0; padding:0 18px; }}
    .workflow-step:first-child {{ padding-left:0; }}
    .workflow-step:not(:last-child)::after {{ content:""; position:absolute; right:0; top:50%; width:22px; height:1px; background:#c4d9ed; }}
    .step-dot {{ display:grid; width:28px; height:28px; flex:0 0 auto; place-items:center; color:#fff; font-size:11px; font-weight:800; border-radius:9px; background:#3673d9; box-shadow:0 5px 12px rgba(54,115,217,.18); }}
    .workflow-step:nth-child(2) .step-dot {{ background:#0ea5b7; }}
    .workflow-step:nth-child(3) .step-dot {{ background:#5366d8; }}
    .workflow-step:nth-child(4) .step-dot {{ background:#15966b; }}
    .step-copy {{ min-width:0; }}
    .step-title {{ display:block; color:#244262; font-size:12px; font-weight:760; white-space:nowrap; }}
    .step-caption {{ display:block; overflow:hidden; color:#8192a8; font-size:10px; text-overflow:ellipsis; white-space:nowrap; }}
    .section-head {{ display:flex; align-items:flex-start; justify-content:space-between; gap:16px; }}
    .section-kicker {{ margin:0 0 3px; color:#6e83a0; font-size:10px; font-weight:800; letter-spacing:.12em; text-transform:uppercase; }}
    .section-head h2 {{ margin-bottom:3px; }}
    .section-meta {{ margin:0; color:#8a9bb0; font-size:11px; white-space:nowrap; }}
    .field-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:0 14px; }}
    .field-grid label {{ margin-top:8px; }}
    .field-grid input {{ width:100%; }}
    .field-span {{ grid-column:1/-1; }}
    .upload-box {{ padding:17px; border:1px dashed #9bb9da; border-radius:14px; background:linear-gradient(145deg,#f8fbff,#f3f9fb); transition:border-color .18s ease,background .18s ease; }}
    .upload-box:focus-within {{ border-color:#4b91ea; background:#f2f8ff; }}
    .upload-icon {{ display:grid; width:34px; height:34px; margin-bottom:9px; place-items:center; color:#2365c8; font-size:18px; font-weight:800; border-radius:10px; background:#e2efff; }}
    .upload-title {{ margin-bottom:3px; color:#294866; font-size:13px; font-weight:780; }}
    .upload-caption {{ margin-bottom:12px; color:#8596aa; font-size:11px; }}
    .upload-box .drop-hint {{ margin-top:11px; }}
    .note {{ color:var(--muted); font-size:11px; line-height:1.55; }}
    .status-waiting {{ color:var(--amber); font-weight:750; }}
    .status-running {{ color:#1d70bc; font-weight:750; }}
    .status-canceling {{ color:#b26b13; font-weight:750; }}
    .status-finished {{ color:var(--green); font-weight:750; }}
    .status-failed {{ color:var(--red); font-weight:750; }}
    .status-canceled {{ color:#78879a; font-weight:750; }}
    .progress {{ min-height:20px; }}
    .upload-summary {{ display:flex; gap:9px; flex-wrap:wrap; margin:12px 0; }}
    .sample-name-list {{ display:grid; gap:9px; margin:12px 0; }}
    .sample-name-row {{ display:grid; grid-template-columns:minmax(0,1fr) minmax(220px,.9fr); gap:12px; align-items:center; padding:11px 12px; border:1px solid #dbe7f3; border-radius:12px; background:#fff; box-shadow:0 4px 12px rgba(31,58,92,.03); }}
    .sample-file-meta {{ display:flex; min-width:0; flex-direction:column; gap:2px; }}
    .sample-file-size {{ color:#1d3554; font-size:16px; font-weight:800; line-height:1.2; }}
    .sample-file-label {{ min-width:0; color:#536a86; font-size:11px; font-weight:650; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }}
    .sample-name-field {{ display:flex; min-width:0; gap:7px; align-items:center; margin:0; color:#536a86; font-size:11px; font-weight:750; white-space:nowrap; }}
    .sample-name-field input {{ min-width:0; flex:1 1 auto; }}
    .sample-name-input {{ width:100%; }}
    @media (max-width:520px) {{ .sample-name-row {{ grid-template-columns:1fr; gap:5px; }} }}
    .metric {{ min-width:78px; padding:9px 12px; color:#73869f; font-size:10px; font-weight:700; border:1px solid #e1eaf4; border-radius:11px; background:#f7fafd; }}
    .metric b {{ display:block; margin-bottom:1px; color:#1d3554; font-size:17px; font-weight:800; line-height:1.25; }}
    #queueSummary .metric:nth-child(1) {{ border-color:#d9e6fa; background:#f2f7ff; }}
    #queueSummary .metric:nth-child(2) {{ border-color:#d2eef0; background:#f0fbfb; }}
    #queueSummary .metric:nth-child(3) {{ border-color:#d6eee4; background:#f2fbf7; }}
    #queueSummary .metric:nth-child(4) {{ border-color:#e0e8f0; background:#f7fafd; }}
    .progress-track {{ height:7px; margin:6px 0 4px; overflow:hidden; min-width:115px; border-radius:10px; background:#e6edf5; }}
    .progress-bar {{ height:100%; border-radius:10px; background:linear-gradient(90deg,#3675e6,#0ea5b7); transition:width .25s ease; }}
    .controls {{ display:flex; gap:10px; align-items:center; margin-top:8px; }}
    .controls span {{ color:#71839a; font-size:12px; white-space:nowrap; }}
    .controls input {{ width:160px; }}
    .action-link, .action-button {{ display:inline-flex; align-items:center; white-space:nowrap; margin:0 6px 4px 0; text-decoration:none; }}
    .action-link {{ padding:7px 10px; color:#205ab6; font-size:11px; font-weight:750; border:1px solid #c9dcf5; border-radius:8px; background:#eff6ff; }}
    .action-link:hover {{ background:#e2efff; }}
    .action-button {{ padding:7px 10px; color:#51647d; font-size:11px; border-color:#d5e0ec; background:#fff; }}
    .action-button:hover {{ color:#1f518f; border-color:#a9c5e5; background:#f5f9ff; }}
    .task-row-selected {{ background:#eef7ff !important; box-shadow:inset 3px 0 0 #3182ce; }}
    .task-list-toolbar {{ display:flex; gap:8px; flex-wrap:wrap; align-items:center; margin:15px 0 11px; }}
    .task-list-toolbar label {{ display:inline; margin:0; color:#667b95; font-size:11px; }}
    .task-search {{ min-width:220px; flex:1 1 260px; }}
    .task-list-toolbar select {{ min-width:110px; }}
    .task-list-toolbar > button {{ padding:9px 11px; color:#60738a; font-size:11px; background:#fff; }}
    .task-pagination {{ display:flex; gap:5px; align-items:center; flex-wrap:wrap; min-height:30px; margin:8px 0 9px; }}
    .task-pagination button {{ min-width:29px; padding:6px 8px; font-size:11px; }}
    .task-pagination button.active {{ color:white; border-color:#2563eb; background:#2563eb; }}
    .task-pagination button:disabled {{ cursor:default; opacity:.45; }}
    .table-wrap {{ width:100%; overflow:hidden; border:1px solid #e5edf5; border-radius:12px; }}
    #taskTable {{ table-layout:fixed; }}
    #taskTable th, #taskTable td {{ white-space:normal; overflow-wrap:anywhere; }}
    #taskTable th:nth-child(1), #taskTable td:nth-child(1) {{ width:27%; }}
    #taskTable th:nth-child(2), #taskTable td:nth-child(2) {{ width:25%; }}
    #taskTable th:nth-child(3), #taskTable td:nth-child(3) {{ width:20%; }}
    #taskTable th:nth-child(4), #taskTable td:nth-child(4) {{ width:13%; }}
    #taskTable th:nth-child(5), #taskTable td:nth-child(5) {{ width:15%; }}
    .task-title {{ color:#1d3554; font-size:13px; font-weight:800; line-height:1.35; }}
    .task-status {{ display:flex; align-items:center; gap:6px; margin-top:5px; font-size:11px; line-height:1.25; }}
    .task-status::before {{ content:""; display:inline-block; width:7px; height:7px; flex:0 0 auto; border-radius:50%; background:currentColor; box-shadow:0 0 0 3px rgba(37,99,235,.10); }}
    .task-queue {{ color:#8a9bb0; font-weight:500; }}
    .sample-files-cell {{ color:#536a86; }}
    .sample-file {{ display:block; max-width:100%; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; line-height:1.5; }}
    .task-actions {{ overflow:visible !important; }}
    .error-text {{ max-width:290px; color:var(--red); font-size:11px; white-space:pre-wrap; }}
    .help-list {{ display:grid; gap:5px; margin:15px 0 0 18px; padding:0; color:#7789a0; font-size:11px; line-height:1.55; }}
    .drop-hint {{ color:#74879f; font-size:11px; line-height:1.5; }}
    #taskDetail {{ min-height:150px; }}
    #detailSummary .upload-summary {{ margin-top:2px; }}
    #detailSummary > div:last-child {{ margin-top:14px; padding:12px 14px; color:#62758d; font-size:11px; line-height:1.8; border:1px solid #e3ebf4; border-radius:11px; background:#f9fbfd; }}
    #detailActions {{ margin-top:14px; }}
    #platformPanel {{ color:#d8e8fb; border-color:#173762; background:linear-gradient(135deg,#0d2040,#122d52); box-shadow:0 18px 42px rgba(10,23,48,.16); }}
    #platformPanel::before {{ background:linear-gradient(90deg,#36b9d0,#66a3ff,transparent); }}
    #platformPanel h2, #platformPanel .section-kicker {{ color:#f3f8ff; }}
    #platformPanel .section-kicker {{ color:#73c9dc; }}
    #platformPanel .note {{ color:#a1b8d3; }}
    #platformStatus .metric {{ min-width:142px; color:#9db5d0; border-color:rgba(145,192,235,.18); background:rgba(255,255,255,.07); }}
    #platformStatus .metric b {{ color:#f1f7ff; font-size:14px; }}
    #platformStatus .metric:first-child b::before {{ content:""; display:inline-block; width:7px; height:7px; margin:0 7px 1px 0; border-radius:50%; background:#47d6a2; box-shadow:0 0 0 4px rgba(71,214,162,.12); }}
    #logPanel {{ padding-bottom:20px; }}
    #taskLog {{ min-height:92px; margin:12px 0 0; white-space:pre-wrap; color:#c8ddfa; font:11px/1.65 "SFMono-Regular",Consolas,"Liberation Mono",monospace; border:1px solid #182f54; border-radius:12px; background:#08172e; box-shadow:inset 0 1px 0 rgba(255,255,255,.04); }}
    @media (max-width:1100px) {{ .grid {{ grid-template-columns:1fr; }} }}
    @media (max-width:760px) {{ header {{ padding:13px 16px; }} .header-note, .nav-divider {{ display:none; }} .top-nav a {{ padding:7px 9px; }} main {{ padding:24px 12px 40px; gap:15px; }} section {{ padding:18px 15px; border-radius:14px; }} .workflow-strip {{ grid-template-columns:1fr 1fr; gap:14px 0; padding:14px; }} .workflow-step {{ padding:0 7px; }} .workflow-step:first-child {{ padding-left:0; }} .workflow-step:not(:last-child)::after {{ display:none; }} .field-grid {{ grid-template-columns:1fr; }} .field-span {{ grid-column:auto; }} .controls {{ align-items:flex-start; flex-direction:column; gap:5px; }} .controls input {{ width:100%; }} .section-head {{ display:block; }} .section-meta {{ display:block; margin-top:5px; }} .task-list-toolbar {{ align-items:stretch; }} .task-list-toolbar label {{ align-self:center; }} .task-search {{ min-width:0; flex-basis:100%; }} #platformStatus .metric {{ min-width:calc(50% - 5px); }} #taskTable th, #taskTable td {{ padding:10px 7px; }} .task-title {{ font-size:12px; }} }}
    @media (max-width:440px) {{ header h1 {{ font-size:15px; }} .brand-kicker {{ display:none; }} .brand-mark {{ width:34px; height:34px; }} .workflow-strip {{ grid-template-columns:1fr; }} .workflow-step {{ padding:0; }} .workflow-step:not(:last-child)::after {{ display:block; right:auto; left:14px; top:30px; width:1px; height:14px; }} .step-caption {{ white-space:normal; }} }}
  </style>
</head>
<body>
  <header>
    <div class="header-brand">
      <span class="brand-mark">∿</span>
      <div>
        <div class="brand-kicker">LC-MS ANALYTICS / INTERNAL</div>
        <h1>LC-MS/MS分析任务平台</h1>
      </div>
    </div>
    <nav class="top-nav" aria-label="主导航">
      <a href="/">任务首页</a>
      <span class="nav-divider"></span>
      <span class="header-note">数据上传&nbsp; → &nbsp;排队运行&nbsp; → &nbsp;结果查看 / 导出</span>
    </nav>
  </header>
  <main>{body}</main>
</body>
</html>""".encode("utf-8")


INDEX_BODY = r"""
<div class="workflow-strip" aria-label="分析流程">
  <div class="workflow-step"><span class="step-dot">01</span><span class="step-copy"><span class="step-title">上传数据</span><span class="step-caption">RAW / mzML / 参考资料</span></span></div>
  <div class="workflow-step"><span class="step-dot">02</span><span class="step-copy"><span class="step-title">进入队列</span><span class="step-caption">按创建时间顺序处理</span></span></div>
  <div class="workflow-step"><span class="step-dot">03</span><span class="step-copy"><span class="step-title">自动分析</span><span class="step-caption">MS1 · MS2 · Feature</span></span></div>
  <div class="workflow-step"><span class="step-dot">04</span><span class="step-copy"><span class="step-title">查看结果</span><span class="step-caption">在线报告与离线导出</span></span></div>
</div>
<div class="grid">
  <section>
    <div class="section-head">
      <div><p class="section-kicker">01 / Analysis setup</p><h2>数据上传与分析参数</h2></div>
    </div>
    <form id="taskForm">
      <div class="field-grid">
        <div class="field-span"><label>任务名称</label><input name="task_name" required placeholder="例如：202407 PTM 两样品比对"></div>
      </div>
      <label>样品数据 <span class="note">至少 2 个 RAW 或 mzML 文件</span></label>
      <div class="upload-box">
        <div class="upload-icon">↑</div>
        <div class="upload-title">选择要比较的原始数据</div>
        <div class="upload-caption">支持多选上传 · 单次上传上限 8 GB</div>
        <input id="files" name="files" type="file" accept=".raw,.RAW,.mzML,.mzml" multiple required>
        <div id="fileSummary" class="sample-name-list" hidden></div>
        <div class="drop-hint">迁移包已自带 ThermoRawFileParser，可直接处理 RAW；mzML 可直接分析。</div>
      </div>
      <div class="field-grid">
        <div class="field-span">
          <label>蛋白质 FASTA 序列 <span class="note">可选 · 文件或文本二选一</span></label>
          <input id="sequence_file" name="sequence_file" type="file" accept=".fa,.fasta,.fas,.faa,.txt">
        </div>
        <div class="field-span">
          <textarea id="sequence_text" name="sequence_text" rows="4" placeholder=">HC\nEVQL...\n>LC\nDIQ..."></textarea>
          <div class="note">未提供序列时仍可完成 MS1 / Feature 分析，但会跳过 MS2 计算。</div>
        </div>
        <div class="field-span">
          <label>蛋白质结构 <span class="note">可选 · 文件或 PDB ID 二选一</span></label>
          <input id="structure_file" name="structure_file" type="file" accept=".pdb,.ent,.cif,.mmcif">
          <div class="controls"><span>或输入 PDB ID</span><input id="pdb_id" name="pdb_id" maxlength="12" placeholder="例如 3V4P"></div>
          <div class="note">只有同时具备 FASTA 与结构时才进行蛋白质结构映射。</div>
        </div>
      </div>
      <div class="form-actions"><button id="submitTask" type="submit">上传并加入队列</button><span class="note">提交后可在右侧队列查看实时进度</span></div>
      <div id="uploadStatus" class="note progress"></div>
      <ul class="help-list">
        <li>上传上限为 8 GB；大文件上传期间请保持页面打开。</li>
        <li>任务按创建时间依次运行，同一台服务器默认一次运行一个任务。</li>
        <li>未提供 FASTA：完成 MS1/Feature 分析，但跳过 MS2 和蛋白质结构映射。</li>
        <li>完成后请通过“查看报告”进入在线报告；MS2 证据保留在主报告内。</li>
      </ul>
    </form>
  </section>
  <section>
    <div class="section-head">
      <div><p class="section-kicker">02 / Queue monitor</p><h2>运行序列与队列</h2></div>
      <p class="section-meta">每 3 秒自动刷新</p>
    </div>
    <div id="queueSummary" class="upload-summary"></div>
    <div class="note">选择任务行查看详情、运行阶段和日志；完成后可直接打开在线报告。</div>
    <div class="task-list-toolbar">
      <label for="taskSearch">搜索任务</label>
       <input id="taskSearch" class="task-search" type="search" placeholder="任务名称、任务 ID、样品名称或文件名">
      <label for="taskStatusFilter">状态</label>
      <select id="taskStatusFilter">
        <option value="">全部状态</option>
        <option value="waiting">排队中</option>
        <option value="running">运行中</option>
        <option value="finished">已完成</option>
        <option value="failed">失败</option>
        <option value="canceled">已取消</option>
      </select>
      <button id="clearTaskFilters" type="button">清除</button>
    </div>
    <div id="taskPagination" class="task-pagination" aria-live="polite"></div>
    <div class="table-wrap"><table id="taskTable"></table></div>
  </section>
</div>
<section id="taskDetail">
  <div class="section-head">
    <div><p class="section-kicker">03 / Selected task</p><h2>选中任务详情与结果</h2></div>
    <p class="section-meta">点击队列中的任意任务</p>
  </div>
  <div id="detailSummary" class="note">请选择一个任务。</div>
  <div id="detailActions"></div>
</section>
<section id="platformPanel">
  <div class="section-head">
    <div><p class="section-kicker">04 / Service health</p><h2>平台状态与运行说明</h2></div>
    <p class="section-meta">Local department service</p>
  </div>
  <div id="platformStatus" class="upload-summary"><span class="metric"><b>读取中…</b>平台状态</span></div>
  <div class="note">本版本不含权限管理；建议把任务目录放在容量充足、定期备份的磁盘上。</div>
</section>
<section id="logPanel">
  <div class="section-head">
    <div><p class="section-kicker">05 / Execution log</p><h2>任务日志</h2></div>
    <p class="section-meta">用于排查转换与分析状态</p>
  </div>
  <div class="note" id="logTitle">点击任务行查看日志。</div>
  <pre id="taskLog"></pre>
</section>
<script>
const $ = id => document.getElementById(id);
let selectedTaskId = "";
let latestTasks = [];
const TASK_PAGE_SIZE = 10;
let taskPage = 1;
const statusText = {waiting:"排队中", running:"运行中", canceling:"正在取消", finished:"已完成", failed:"失败", canceled:"已取消"};
function escapeHtml(text){ return String(text ?? "").replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch])); }
function formatBytes(value){ const n=Number(value||0); if(!n) return "0 B"; const units=["B","KB","MB","GB","TB"]; let i=0, v=n; while(v>=1024&&i<units.length-1){v/=1024;i++;} return `${v.toFixed(i?1:0)} ${units[i]}`; }
 function formatDuration(value){ const n=Math.max(0,Math.round(Number(value||0))); if(!n) return "—"; const h=Math.floor(n/3600), m=Math.floor((n%3600)/60), s=n%60; return h?`${h}h ${m}m ${s}s`:m?`${m}m ${s}s`:`${s}s`; }
 function referenceSummary(task){
   const refs=task.references||{}, seq=refs.sequence||{}, structure=refs.structure||{};
   const sequence=seq.provided?`FASTA：${seq.source||"已提供"}`:"FASTA：未提供";
   const structureText=structure.provided?`结构：${structure.source==="pdb_id"?`PDB ${structure.pdb_id||""}`:(structure.name||"已提供")}`:"结构：未提供";
   return `${sequence}；${structureText}`;
 }
function taskResultLinks(task){
  if(task.status !== "finished") return "";
  return `<a class="action-link" href="${escapeHtml(task.report_url||"")}">查看报告</a>`;
}
function taskActionButton(task){
  if(task.status === "waiting") return `<button class="action-button" data-task-action="cancel" data-task-id="${escapeHtml(task.task_id)}">取消排队</button>`;
  if(task.status === "running") return `<button class="action-button" data-task-action="cancel" data-task-id="${escapeHtml(task.task_id)}">取消任务</button>`;
  if(task.status === "failed") return `<button class="action-button" data-task-action="retry" data-task-id="${escapeHtml(task.task_id)}">失败重试</button>`;
  return "";
}
async function fetchJson(url, options){
  const response = await fetch(url, {cache:"no-store", ...(options||{})});
  if(!response.ok) throw new Error(await response.text());
  return await response.json();
}
function renderTaskDetail(){
  const task=latestTasks.find(item=>item.task_id===selectedTaskId);
  if(!task){ $("detailSummary").textContent="请选择一个任务。"; $("detailActions").innerHTML=""; return; }
  const files=(task.files||[]).map(file=>`${escapeHtml(file.sample_name||file.original_name)} ← ${escapeHtml(file.original_name)}（${formatBytes(file.size_bytes)}）`).join("；");
  const queue=task.status==="waiting"&&task.queue_position?`队列第 ${task.queue_position} 位`:(statusText[task.status]||task.status);
  const error=task.error?`<div class="error-text">错误：${escapeHtml(task.error)}</div>`:"";
  $("detailSummary").innerHTML=`<div class="upload-summary"><span class="metric"><b>${escapeHtml(statusText[task.status]||task.status)}</b>状态</span><span class="metric"><b>${Number(task.progress||0)}%</b>进度</span><span class="metric"><b>${escapeHtml(queue)}</b>队列</span><span class="metric"><b>${escapeHtml(formatDuration(task.duration_seconds))}</b>耗时</span><span class="metric"><b>${formatBytes(task.total_size_bytes)}</b>上传数据</span></div><div>任务 ID：${escapeHtml(task.task_id)}<br>样品文件：${files||"—"}<br>参考资料：${escapeHtml(referenceSummary(task))}<br>阶段：${escapeHtml(task.stage||"—")}<br>创建时间：${escapeHtml(task.created_at||"—")}</div>${error}`;
  $("detailActions").innerHTML=`<p>${taskResultLinks(task)} ${taskActionButton(task)}</p>`;
  bindTaskActions();
}
function bindTaskActions(){
  document.querySelectorAll("[data-task-action]").forEach(button=>button.onclick=async event=>{ event.stopPropagation(); await taskAction(button.dataset.taskAction, button.dataset.taskId); });
}
async function loadPlatformStatus(){
  try {
    const payload=await fetchJson("/api/status");
    const parser=payload.parser_available?"可用":"未配置（mzML 仍可用）";
    $("platformStatus").innerHTML=`<span class="metric"><b>${escapeHtml(payload.service||"在线")}</b>平台状态</span><span class="metric"><b>${escapeHtml(payload.python_version||"—")}</b>Python</span><span class="metric"><b>${escapeHtml(parser)}</b>RAW 转换器</span><span class="metric"><b>${formatBytes(payload.max_upload_bytes)}</b>上传上限</span><span class="metric"><b>${escapeHtml(payload.jobs_dir||"—")}</b>任务目录</span>`;
  } catch(err){ $("platformStatus").innerHTML=`<span class="metric"><b>读取失败</b>${escapeHtml(err.message||err)}</span>`; }
}
async function taskAction(action, taskId){
  try { await fetchJson(`/api/tasks/${encodeURIComponent(taskId)}/${action}`, {method:"POST"}); selectedTaskId=taskId; await loadTasks(); await loadLog(); }
  catch(err){ $("uploadStatus").textContent=`操作失败：${err.message||err}`; }
}
function taskSearchText(task){
  const refs=task.references||{}, sequence=refs.sequence||{}, structure=refs.structure||{};
  return [task.task_name, task.task_id, task.status, statusText[task.status], task.stage,
     ...(task.files||[]).flatMap(file=>[file.sample_name, file.original_name, file.name]),
    sequence.source, structure.source, structure.pdb_id]
    .filter(Boolean).join(" ").toLocaleLowerCase();
}
function filteredTasks(){
  const query=$("taskSearch").value.trim().toLocaleLowerCase();
  const status=$("taskStatusFilter").value;
  return latestTasks.filter(task=>(!status||task.status===status)&&(!query||taskSearchText(task).includes(query)));
}
function renderTaskPagination(total){
  const container=$("taskPagination");
  const totalPages=Math.max(1, Math.ceil(total/TASK_PAGE_SIZE));
  taskPage=Math.min(Math.max(1, taskPage), totalPages);
  if(!total){ container.innerHTML='<span class="note">没有匹配的任务</span>'; return; }
  const pageSet=new Set([1,totalPages,taskPage,taskPage-1,taskPage+1]);
  const pages=[...pageSet].filter(page=>page>=1&&page<=totalPages).sort((a,b)=>a-b);
  let pageButtons="", previous=0;
  pages.forEach(page=>{
    if(previous&&page-previous>1) pageButtons+='<span class="note">…</span>';
    pageButtons+=`<button type="button" class="${page===taskPage?"active":""}" data-task-page="${page}">${page}</button>`;
    previous=page;
  });
  container.innerHTML=`<button type="button" data-task-page="${taskPage-1}" ${taskPage<=1?"disabled":""}>上一页</button>${pageButtons}<button type="button" data-task-page="${taskPage+1}" ${taskPage>=totalPages?"disabled":""}>下一页</button><span class="note">第 ${taskPage}/${totalPages} 页，共 ${total} 条</span>`;
  container.querySelectorAll("[data-task-page]").forEach(button=>button.onclick=()=>{ if(button.disabled) return; taskPage=Number(button.dataset.taskPage); renderTaskRows(); });
}
function renderTaskRows(){
  const matchingTasks=filteredTasks();
  const totalPages=Math.max(1, Math.ceil(matchingTasks.length/TASK_PAGE_SIZE));
  taskPage=Math.min(Math.max(1, taskPage), totalPages);
  renderTaskPagination(matchingTasks.length);
  const start=(taskPage-1)*TASK_PAGE_SIZE;
  const rows = matchingTasks.slice(start,start+TASK_PAGE_SIZE).map(task => {
    const progress=Math.max(0,Math.min(100,Number(task.progress||0)));
    const queue=task.status==="waiting"&&task.queue_position?`第 ${task.queue_position} 位`:"—";
     const files=(task.files||[]).map(file=>`<span class="sample-file" title="${escapeHtml(file.sample_name||file.original_name)}"><b>${escapeHtml(file.sample_name||file.original_name)}</b> <span class="note">← ${escapeHtml(file.original_name)} · ${formatBytes(file.size_bytes)}</span></span>`).join("");
    const error=task.error?`<div class="error-text">${escapeHtml(task.error)}</div>`:"";
    return `<tr data-task="${escapeHtml(task.task_id)}" class="${task.task_id===selectedTaskId?"task-row-selected":""}"><td><div class="task-title">${escapeHtml(task.task_name)}</div><div class="task-status status-${escapeHtml(task.status)}">${escapeHtml(statusText[task.status]||task.status)}${queue!=="—"?`<span class="task-queue">· ${escapeHtml(queue)}</span>`:""}</div></td><td>${escapeHtml(task.stage||"—")}<div class="progress-track"><div class="progress-bar" style="width:${progress}%"></div></div><span class="note">${progress}% · ${escapeHtml(formatDuration(task.duration_seconds))}</span>${error}</td><td class="sample-files-cell">${files||"—"}</td><td>${escapeHtml(task.created_at||"—")}</td><td class="task-actions">${taskResultLinks(task)} ${taskActionButton(task)}</td></tr>`;
  }).join("");
   $("taskTable").innerHTML=`<tr><th>任务</th><th>阶段 / 进度</th><th>样品名称 / 文件</th><th>创建时间</th><th>结果 / 操作</th></tr>${rows||'<tr><td colspan="5" class="note">暂无匹配任务。</td></tr>'}`;
  document.querySelectorAll("tr[data-task]").forEach(row=>row.onclick=()=>{ selectedTaskId=row.dataset.task; renderTaskDetail(); loadLog(); renderTaskRows(); });
  bindTaskActions();
  renderTaskDetail();
}
async function loadTasks(){
  const payload = await fetchJson("/api/tasks");
  latestTasks = payload.tasks || [];
  const waiting=latestTasks.filter(task=>task.status==="waiting").length;
  const running=latestTasks.filter(task=>task.status==="running").length;
  const finished=latestTasks.filter(task=>task.status==="finished").length;
  $("queueSummary").innerHTML=`<span class="metric"><b>${waiting}</b>排队</span><span class="metric"><b>${running}</b>运行中</span><span class="metric"><b>${finished}</b>已完成</span><span class="metric"><b>${latestTasks.length}</b>全部任务</span>`;
  renderTaskRows();
}
async function loadLog(){
  if(!selectedTaskId) return;
  $("logTitle").textContent = `任务 ${selectedTaskId} 日志`;
  try { const payload=await fetchJson(`/api/tasks/${encodeURIComponent(selectedTaskId)}/log`); $("taskLog").textContent=payload.log||""; }
  catch(err){ $("taskLog").textContent=String(err.message||err); }
}
function renderSampleNameInputs(){
  const files=[...($("files").files||[])];
  const container=$("fileSummary");
  container.hidden=!files.length;
  container.innerHTML=files.map((file,index)=>`<div class="sample-name-row"><div class="sample-file-meta"><span class="sample-file-size">${formatBytes(file.size)}</span><span class="sample-file-label" title="${escapeHtml(file.name)}">${escapeHtml(file.name)}</span></div><label class="sample-name-field">填写样品名<input class="sample-name-input" name="sample_name_${index}" value="${escapeHtml(file.name)}" maxlength="160" required aria-label="${escapeHtml(file.name)} 的样品名称" placeholder="请输入样品名"></label></div>`).join("");
}
$("files").onchange=renderSampleNameInputs;
$("taskForm").onsubmit = async event => {
  event.preventDefault();
  const files=[...($("files").files||[])];
  if(files.length<2){ $("uploadStatus").textContent="请至少选择两个 RAW 或 mzML 文件。"; return; }
  const invalid=files.find(file=>!(/\.(raw|mzml)$/i).test(file.name));
  if(invalid){ $("uploadStatus").textContent=`文件类型不支持：${invalid.name}`; return; }
  const sequenceFile=$("sequence_file").files?.[0], sequenceText=$("sequence_text").value.trim();
  if(sequenceFile && sequenceText){ $("uploadStatus").textContent="FASTA 文件和粘贴序列只能选择一种。"; return; }
  const structureFile=$("structure_file").files?.[0], pdbId=$("pdb_id").value.trim();
  if(structureFile && pdbId){ $("uploadStatus").textContent="结构文件和 PDB ID 只能选择一种。"; return; }
  if(pdbId && !/^[A-Za-z0-9]{4,12}$/.test(pdbId)){ $("uploadStatus").textContent="PDB ID 必须为 4-12 位字母或数字。"; return; }
   const data = new FormData(event.target);
   const sampleNameInputs=[...document.querySelectorAll(".sample-name-input")];
   const sampleNames=sampleNameInputs.map(input=>input.value.trim());
   if(sampleNames.length!==files.length||sampleNames.some(name=>!name)){ $("uploadStatus").textContent="请为每个文件填写样品名称。"; return; }
   if(new Set(sampleNames.map(name=>name.toLocaleLowerCase())).size!==sampleNames.length){ $("uploadStatus").textContent="样品名称不能重复，请分别填写唯一名称。"; return; }
   sampleNameInputs.forEach((_,index)=>data.delete(`sample_name_${index}`));
   data.append("sample_names",JSON.stringify(sampleNames));
  $("submitTask").disabled=true;
  $("uploadStatus").textContent = "正在上传，RAW 文件较大时请耐心等待...";
  try {
    const payload = await fetchJson("/api/tasks", {method:"POST", body:data});
    $("uploadStatus").textContent = `已加入队列：${payload.task.task_id}`;
     event.target.reset(); $("fileSummary").innerHTML=""; $("fileSummary").hidden=true;
    selectedTaskId = payload.task.task_id;
    await loadTasks(); await loadLog();
 } catch (err) { $("uploadStatus").textContent = `上传失败：${err.message || err}`; }
 finally { $("submitTask").disabled=false; }
 };
 $("sequence_file").onchange=event=>{ if(event.target.files?.length) $("sequence_text").value=""; };
 $("sequence_text").oninput=event=>{ if(event.target.value.trim()) $("sequence_file").value=""; };
 $("structure_file").onchange=event=>{ if(event.target.files?.length) $("pdb_id").value=""; };
 $("pdb_id").oninput=event=>{ if(event.target.value.trim()) $("structure_file").value=""; };
$("taskSearch").oninput=()=>{ taskPage=1; renderTaskRows(); };
$("taskStatusFilter").onchange=()=>{ taskPage=1; renderTaskRows(); };
$("clearTaskFilters").onclick=()=>{ $("taskSearch").value=""; $("taskStatusFilter").value=""; taskPage=1; renderTaskRows(); };
loadTasks().catch(err => $("taskTable").innerHTML = `<tr><td>${escapeHtml(err.message || err)}</td></tr>`);
loadPlatformStatus();
setInterval(() => { loadTasks().catch(()=>{}); if(selectedTaskId) loadLog(); }, 3000);
</script>
"""


def task_output_dir(config: PortalConfig, task_id: str) -> Path:
    return config.jobs_dir / task_id / "output"


def task_db_path(config: PortalConfig, task_id: str) -> Path:
    return task_output_dir(config, task_id) / "lcms_peak_first_compare.sqlite"


def task_html_path(config: PortalConfig, task_id: str) -> Path:
    return task_output_dir(config, task_id) / "lcms_peak_first_compare.html"


def rewrite_report_html(task_id: str, content: str) -> str:
    prefix = f"/tasks/{quote(task_id)}/report"
        # The report shell is generated from the current template so previously
        # completed tasks also receive the current unified MS1/MS2 interface.
    content = PEAK_FIRST_TEMPLATE
    content = content.replace('"/api/', f'"{prefix}/api/')
    content = content.replace("'/api/", f"'{prefix}/api/")
    content = content.replace("`/api/", f"`{prefix}/api/")
    content = content.replace('src="/assets/', f'src="{prefix}/assets/')
    return content


def _read_artifact_raw(db_path: Path, key: str) -> object:
    return read_artifact(db_path, key)


def _read_optional_artifact(db_path: Path, key: str) -> object | None:
    try:
        return read_artifact(db_path, key)
    except KeyError:
        return None


def standalone_offline_script(task_id: str, db_path: Path) -> str:
    bootstrap = _read_artifact_raw(db_path, "bootstrap")
    if not isinstance(bootstrap, dict):
        raise ValueError("bootstrap artifact is not a JSON object")
    sample_ids = [str(item) for item in bootstrap.get("sample_ids", [])]
    spectra = {sample_id: _read_artifact_raw(db_path, f"spectra:{sample_id}") for sample_id in sample_ids}
    payload = {
        "task_id": task_id,
        "bootstrap": bootstrap,
        "features": read_features(db_path),
        "spectra": spectra,
        "msms": _read_optional_artifact(db_path, "msms_identifications"),
        "task_reference_structure": _read_optional_artifact(db_path, "task_reference_structure"),
    }
    compressed = gzip.compress(json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode("utf-8"), compresslevel=6)
    encoded = base64.b64encode(compressed).decode("ascii")
    return f"""
<script>
window.__LCMS_OFFLINE_REPORT__ = {{
  taskId: {json.dumps(task_id, ensure_ascii=False)},
  compressedBase64: "{encoded}"
}};
(function(){{
  const offline = window.__LCMS_OFFLINE_REPORT__;
  let dataPromise = null;
  function b64ToBytes(text){{
    const binary = atob(text);
    const bytes = new Uint8Array(binary.length);
    for(let i=0;i<binary.length;i++) bytes[i] = binary.charCodeAt(i);
    return bytes;
  }}
  async function loadData(){{
    if(dataPromise) return dataPromise;
    dataPromise = (async()=>{{
      const bytes = b64ToBytes(offline.compressedBase64);
      if(!("DecompressionStream" in window)){{
        throw new Error("当前浏览器不支持离线压缩数据解码，请使用新版 Chrome 或 Edge 打开此 HTML。");
      }}
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
      const text = await new Response(stream).text();
      return JSON.parse(text);
    }})();
    return dataPromise;
  }}
  function jsonResponse(payload, status=200){{
    return new Response(JSON.stringify(payload), {{status, headers:{{"Content-Type":"application/json; charset=utf-8"}}}});
  }}
  function mzTolerance(mz, params){{
    if(String(params.mz_tolerance_mode || "da") === "ppm") return Math.max(mz * Number(params.mz_tolerance_ppm || 20) / 1000000, 1e-9);
    return Number(params.mz_tolerance_da || 0.5);
  }}
  function integrate(rt, intensity){{
    let area = 0;
    for(let i=1;i<rt.length;i++) area += (intensity[i-1] + intensity[i]) * 0.5 * Math.max(rt[i] - rt[i-1], 0);
    return area;
  }}
  function xicPayload(store, peakId, targetMz, fullRun){{
    const payload = store.bootstrap;
    const params = payload.params || {{}};
    const peak = (payload.peak_results || []).find(item => item && item.tic_peak_id === peakId);
    if(!peak) throw new Error(`peak not found: ${{peakId}}`);
    const tolerance = mzTolerance(targetMz, params);
    let featureGroup = null;
    const groups = Array.isArray(peak.feature_groups) ? peak.feature_groups : [];
    if(groups.length){{
      featureGroup = groups.reduce((best, group) => {{
        const delta = Math.abs(Number(group.representative_mz || targetMz) - targetMz);
        return !best || delta < best.delta ? {{group, delta}} : best;
      }}, null);
      featureGroup = featureGroup && featureGroup.delta <= Math.max(tolerance * 2, 1e-9) ? featureGroup.group : null;
    }}
    const featureBoundsBySample = {{}};
    if(featureGroup && featureGroup.features_by_sample){{
      for(const [sample, feature] of Object.entries(featureGroup.features_by_sample)){{
        if(feature && feature.rt_start != null && feature.rt_end != null){{
          featureBoundsBySample[sample] = {{
            rt_start: Number(feature.rt_start),
            rt_end: Number(feature.rt_end),
            rt_apex: Number(feature.rt_apex || feature.aligned_rt_apex || feature.rt_start),
            area: Number(feature.area || 0),
            height: Number(feature.height || 0),
            match_status: String(feature.match_status || "")
          }};
        }}
      }}
    }}
    const featureBounds = Object.values(featureBoundsBySample);
    const integrationRtStart = featureBounds.length ? Math.min(...featureBounds.map(item => item.rt_start)) : Number(peak.rt_start || 0);
    const integrationRtEnd = featureBounds.length ? Math.max(...featureBounds.map(item => item.rt_end)) : Number(peak.rt_end || 0);
    const contextMin = Math.max(Number(params.xic_context_min || 3), 0.5);
    const contextRtStart = Math.max(0, integrationRtStart - contextMin);
    const contextRtEnd = integrationRtEnd + contextMin;
    const sampleIds = (payload.sample_ids || []).map(String);
    const sampleCutBounds = peak.sample_cut_bounds || {{}};
    const xicBySample = {{}};
    const integrationBySample = {{}};
    for(const sample of sampleIds){{
      const scans = store.spectra[sample] || [];
      const sampleBounds = sampleCutBounds[sample] || {{}};
      const localShift = Number(sampleBounds.local_rt_shift || 0);
      const fb = featureBoundsBySample[sample];
      const sampleStart = fb ? fb.rt_start : integrationRtStart;
      const sampleEnd = fb ? fb.rt_end : integrationRtEnd;
      const rtValues = [];
      const intensities = [];
      for(const scan of scans){{
        const rt = Number(scan.aligned_rt || scan.rt || 0) + localShift;
        if(!fullRun && (rt < contextRtStart || rt > contextRtEnd)) continue;
        let total = 0;
        const mz = scan.mz || [];
        const intensity = scan.intensity || [];
        for(let i=0;i<mz.length;i++) if(Math.abs(Number(mz[i]) - targetMz) <= tolerance) total += Number(intensity[i] || 0);
        rtValues.push(rt);
        intensities.push(total);
      }}
      const intRt = [];
      const intY = [];
      for(let i=0;i<rtValues.length;i++){{
        if(rtValues[i] >= sampleStart && rtValues[i] <= sampleEnd){{
          intRt.push(rtValues[i]);
          intY.push(intensities[i]);
        }}
      }}
      const calcArea = integrate(intRt, intY);
      const height = intY.length ? Math.max(...intY) : 0;
      const apexIndex = intY.indexOf(height);
      const apex = apexIndex >= 0 && height > 0 ? intRt[apexIndex] : null;
      xicBySample[sample] = {{rt: rtValues, intensity: intensities}};
      integrationBySample[sample] = {{
        area: fb ? fb.area : calcArea,
        height: fb ? fb.height : height,
        rt_apex: fb ? fb.rt_apex : apex,
        rt_start: sampleStart,
        rt_end: sampleEnd,
        local_rt_shift: localShift,
        match_status: fb ? fb.match_status : ""
      }};
    }}
    const allRt = Object.values(xicBySample).flatMap(series => series.rt || []);
    return {{
      peak_id: peakId,
      feature_group_id: featureGroup ? featureGroup.feature_group_id : null,
      full_run: fullRun,
      target_mz: targetMz,
      mz_tolerance: tolerance,
      rt_start: fullRun && allRt.length ? Math.min(...allRt) : integrationRtStart,
      rt_end: fullRun && allRt.length ? Math.max(...allRt) : integrationRtEnd,
      integration_rt_start: integrationRtStart,
      integration_rt_end: integrationRtEnd,
      context_rt_start: contextRtStart,
      context_rt_end: contextRtEnd,
      sample_cut_bounds: sampleCutBounds,
      feature_bounds_by_sample: featureBoundsBySample,
      xic_by_sample: xicBySample,
      integration_by_sample: integrationBySample
    }};
  }}
  const originalFetch = window.fetch.bind(window);
  window.fetch = async function(input, init={{}}){{
    const urlText = typeof input === "string" ? input : input.url;
    const url = new URL(urlText, window.location.href);
    const path = url.pathname;
    if(!path.includes("/api/")) return originalFetch(input, init);
    const store = await loadData();
    if(path.endsWith("/api/comparisons")) return jsonResponse({{default_comparison: offline.taskId, comparisons:[{{id: offline.taskId, label: offline.taskId}}]}});
    if(path.endsWith("/api/bootstrap")) return jsonResponse(store.bootstrap);
    if(path.endsWith("/api/msms")) return store.msms ? jsonResponse(store.msms) : jsonResponse({{error:"MS/MS identification has not been run"}}, 404);
    if(path.endsWith("/api/task-reference")) return store.task_reference_structure ? jsonResponse(store.task_reference_structure) : jsonResponse({{available:false,reason:"未提供结构文件或 PDB ID。"}});
    if(path.endsWith("/api/pdb-structure")){{
      const pdbId = String(url.searchParams.get("pdb_id") || "").toUpperCase();
      if(!/^[A-Z0-9]{{4,12}}$/.test(pdbId)) return jsonResponse({{error:"invalid PDB ID"}}, 400);
      return originalFetch(`https://files.rcsb.org/download/${{pdbId}}.cif`, {{cache:"force-cache"}});
    }}
    if(path.endsWith("/api/features")){{
      if(String(init.method || "GET").toUpperCase() === "POST"){{
        const body = JSON.parse(init.body || "{{}}");
        store.features = Array.isArray(body.features) ? body.features : [];
        try {{ localStorage.setItem(`lcms_features_${{offline.taskId}}`, JSON.stringify(store.features)); }} catch(_err) {{}}
        return jsonResponse({{saved_count: store.features.length, features: store.features}});
      }}
      try {{
        const saved = localStorage.getItem(`lcms_features_${{offline.taskId}}`);
        if(saved) store.features = JSON.parse(saved);
      }} catch(_err) {{}}
      return jsonResponse({{features: store.features || []}});
    }}
    if(path.endsWith("/api/xic")){{
      const peakId = url.searchParams.get("peak_id") || "";
      const mz = Number(url.searchParams.get("mz"));
      const fullRun = ["1","true","yes","full"].includes(String(url.searchParams.get("full") || "0").toLowerCase());
      return jsonResponse(xicPayload(store, peakId, mz, fullRun));
    }}
    return jsonResponse({{error:`offline API not found: ${{path}}`}}, 404);
  }};
}})();
</script>
"""


def build_standalone_report(config: PortalConfig, task_id: str) -> Path:
    task_dir = config.jobs_dir / task_id
    output_dir = task_output_dir(config, task_id)
    html_path = task_html_path(config, task_id)
    db_path = task_db_path(config, task_id)
    if not html_path.exists() or not db_path.exists():
        raise FileNotFoundError("report HTML or SQLite output is missing")
    server_path = WORKSPACE / "lcms_feature_mvp" / "serve_peak_first_compare.py"
    export_dir = task_dir / "exports"
    export_dir.mkdir(parents=True, exist_ok=True)
    standalone_path = export_dir / f"{task_id}_standalone_interactive_report.html"
    viewer_asset = VENDOR_DIR / "3Dmol-min.js"
    source_mtime = max(html_path.stat().st_mtime_ns, db_path.stat().st_mtime_ns, viewer_asset.stat().st_mtime_ns, server_path.stat().st_mtime_ns)
    if standalone_path.exists() and standalone_path.stat().st_mtime_ns > source_mtime:
        return standalone_path
    content = PEAK_FIRST_TEMPLATE
    if viewer_asset.exists():
        viewer_javascript = viewer_asset.read_text(encoding="utf-8").replace("</script", "<\\/script")
        content = content.replace(
            '<script src="/assets/3Dmol-min.js"></script>',
            f"<script>{viewer_javascript}</script>",
            1,
        )
    offline_script = standalone_offline_script(task_id, db_path)
    if "<script>\nlet DATA" in content:
        content = content.replace("<script>\nlet DATA", offline_script + "\n<script>\nlet DATA", 1)
    else:
        content = content.replace("</head>", offline_script + "\n</head>", 1)
    badge = '<span style="color:#bfdbfe">单文件离线报告，无需启动端口服务</span>'
    content = content.replace("</header>", f"  {badge}\n</header>", 1)
    standalone_path.write_text(content, encoding="utf-8")
    return standalone_path


def build_static_offline_package(config: PortalConfig, task_id: str) -> Path:
    """Build a no-server offline report package.

    The package must be extracted first, then `index.html` can be opened
    directly. Large data lives in `data/lcms_offline_data.js` instead of being
    embedded in the HTML file.
    """
    task_dir = config.jobs_dir / task_id
    html_path = task_html_path(config, task_id)
    db_path = task_db_path(config, task_id)
    if not html_path.exists() or not db_path.exists():
        raise FileNotFoundError("report HTML or SQLite output is missing")
    server_path = WORKSPACE / "lcms_feature_mvp" / "serve_peak_first_compare.py"
    export_dir = task_dir / "exports"
    export_dir.mkdir(parents=True, exist_ok=True)
    zip_path = export_dir / f"{task_id}_static_offline_report.zip"
    viewer_asset = VENDOR_DIR / "3Dmol-min.js"
    source_mtime = max(html_path.stat().st_mtime_ns, db_path.stat().st_mtime_ns, viewer_asset.stat().st_mtime_ns)
    if zip_path.exists() and zip_path.stat().st_mtime_ns > source_mtime and zip_path.stat().st_mtime_ns > server_path.stat().st_mtime_ns:
        return zip_path

    # Reuse the tested standalone generator, then move its large base64 payload
    # into a sibling JS asset so index.html stays lightweight.
    standalone_path = build_standalone_report(config, task_id)
    content = standalone_path.read_text(encoding="utf-8")
    match = re.search(r'compressedBase64:\s*"([A-Za-z0-9+/=]+)"', content)
    if not match:
        raise ValueError("standalone report payload was not found")
    encoded = match.group(1)
    content = content[: match.start(1) - 1] + "window.__LCMS_OFFLINE_DATA_COMPRESSED__" + content[match.end(1) + 1 :]
    data_script_tag = '<script src="data/lcms_offline_data.js"></script>\n'
    content = content.replace("<script>\nwindow.__LCMS_OFFLINE_REPORT__", data_script_tag + "<script>\nwindow.__LCMS_OFFLINE_REPORT__", 1)
    content = content.replace("单文件离线报告，无需启动端口服务", "静态离线报告包，无需启动端口服务")

    data_js = (
        "// LC-MS offline report data. Keep this file next to index.html in the extracted package.\n"
        f"window.__LCMS_OFFLINE_DATA_COMPRESSED__ = \"{encoded}\";\n"
    )
    readme = f"""# LC-MS 静态离线报告包

任务 ID: {task_id}

## 使用方法

1. 先完整解压这个 zip。
2. 双击打开 `index.html`。
3. 不需要启动端口服务。

请不要只在压缩包预览窗口中双击 HTML；需要先解压，否则浏览器可能无法加载旁边的 `data/lcms_offline_data.js`。

## 文件说明

- `index.html`: 轻量报告前端
- `data/lcms_offline_data.js`: 压缩后的报告数据

这个包不包含 RAW/mzML 原始数据。
"""
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        archive.writestr("index.html", content)
        archive.writestr("data/lcms_offline_data.js", data_js)
        archive.writestr("README.md", readme)
    return zip_path


def package_readme(task_id: str) -> str:
    return f"""# LC-MS 可交互报告包

任务 ID: {task_id}

## 使用方法

Windows:

```powershell
powershell -ExecutionPolicy Bypass -File .\\start_report_server.ps1
```

或者双击 `start_report_server.bat`。

启动后打开：

```text
http://127.0.0.1:8780/
```

## 包含内容

- `output/lcms_peak_first_compare.html`: Peak-first 交互报告页面
- `output/lcms_peak_first_compare.sqlite`: 报告所需的本地结果库
- `serve_peak_first_compare.py`: 轻量本地查看服务器

这个包不包含 RAW/mzML 原始数据，因此体积更小。XIC、Feature 表、热图和保存 Feature 等交互功能由 SQLite 结果库支持。
"""


def build_report_package(config: PortalConfig, task_id: str) -> Path:
    task_dir = config.jobs_dir / task_id
    output_dir = task_output_dir(config, task_id)
    html_path = task_html_path(config, task_id)
    db_path = task_db_path(config, task_id)
    server_path = WORKSPACE / "lcms_feature_mvp" / "serve_peak_first_compare.py"
    viewer_asset = VENDOR_DIR / "3Dmol-min.js"
    if not html_path.exists() or not db_path.exists():
        raise FileNotFoundError("report HTML or SQLite output is missing")
    if not server_path.exists():
        raise FileNotFoundError(f"report server script is missing: {server_path}")
    export_dir = task_dir / "exports"
    export_dir.mkdir(parents=True, exist_ok=True)
    zip_path = export_dir / f"{task_id}_interactive_report.zip"
    references_dir = task_dir / "references"
    reference_files = [path for path in references_dir.glob("*") if path.is_file()] if references_dir.exists() else []
    sources = [html_path, db_path, server_path, viewer_asset, *reference_files]
    source_mtime = max(item.stat().st_mtime_ns for item in sources)
    if zip_path.exists() and zip_path.stat().st_mtime_ns > source_mtime:
        return zip_path
    start_ps1 = """param([int]$Port = 8780)
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
python (Join-Path $Root "serve_peak_first_compare.py") --output-dir (Join-Path $Root "output") --port $Port
"""
    start_bat = """@echo off
cd /d "%~dp0"
powershell -ExecutionPolicy Bypass -File "%~dp0start_report_server.ps1"
pause
"""
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        archive.writestr("output/lcms_peak_first_compare.html", PEAK_FIRST_TEMPLATE)
        archive.write(db_path, "output/lcms_peak_first_compare.sqlite")
        archive.write(server_path, "serve_peak_first_compare.py")
        archive.write(viewer_asset, "ui/vendor/3Dmol-min.js")
        for reference_path in reference_files:
            archive.write(reference_path, f"references/{reference_path.name}")
        archive.writestr("start_report_server.ps1", start_ps1)
        archive.writestr("start_report_server.bat", start_bat)
        archive.writestr("README.md", package_readme(task_id))
        metadata_path = task_dir / "worker.log"
        if metadata_path.exists():
            archive.write(metadata_path, "worker.log")
    return zip_path


def make_handler(config: PortalConfig, store: TaskStore) -> type[BaseHTTPRequestHandler]:
    mcp = MCPApplication(config, store, audit=lambda message: print(f"MCP {message}"))

    class Handler(BaseHTTPRequestHandler):
        server_version = "LCMSDepartmentPlatform/0.1"

        def log_message(self, fmt: str, *args: object) -> None:
            print(f"{self.address_string()} - {fmt % args}")

        def send_bytes(self, body: bytes, content_type: str, status: HTTPStatus = HTTPStatus.OK) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def send_json(self, payload: object, status: HTTPStatus = HTTPStatus.OK) -> None:
            self.send_bytes(json.dumps(payload, ensure_ascii=False).encode("utf-8"), "application/json; charset=utf-8", status)

        def send_error_json(self, message: str, status: HTTPStatus = HTTPStatus.BAD_REQUEST) -> None:
            self.send_json({"error": message}, status)

        def send_file(self, path: Path, content_type: str, download_name: str | None = None) -> None:
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(path.stat().st_size))
            if download_name:
                self.send_header("Content-Disposition", f'attachment; filename="{download_name}"')
            self.end_headers()
            with path.open("rb") as handle:
                shutil.copyfileobj(handle, self.wfile)

        def task_from_report_path(self, path: str) -> tuple[str, str] | None:
            parts = path.strip("/").split("/")
            if len(parts) >= 3 and parts[0] == "tasks" and parts[2] == "report":
                task_id = unquote(parts[1])
                rest = "/" + "/".join(parts[3:]) if len(parts) > 3 else "/"
                return task_id, rest
            return None

        def do_GET(self) -> None:  # noqa: N802
            parsed = urlparse(self.path)
            path = parsed.path
            try:
                if path == "/":
                    self.send_bytes(page_shell("LC-MS/MS分析任务平台", INDEX_BODY), "text/html; charset=utf-8")
                    return
                if path == "/api/tasks":
                    self.send_json({"tasks": task_snapshots(store.load())})
                    return
                if path == "/api/status":
                    self.send_json(
                        {
                            "service": "online",
                            "python_version": sys.version.split()[0],
                            "parser_available": config.parser_path.exists(),
                            "parser_path": str(config.parser_path),
                            "max_upload_bytes": config.max_upload_bytes,
                            "jobs_dir": str(config.jobs_dir),
                            "mcp_endpoint": "/mcp",
                            "mcp_read_only": True,
                        }
                    )
                    return
                if path == "/mcp":
                    body = json.dumps(
                        {"error": "MCP endpoint accepts POST JSON-RPC requests"},
                        ensure_ascii=False,
                    ).encode("utf-8")
                    self.send_response(HTTPStatus.METHOD_NOT_ALLOWED)
                    self.send_header("Allow", "POST")
                    self.send_header("Content-Type", "application/json; charset=utf-8")
                    self.send_header("Content-Length", str(len(body)))
                    self.end_headers()
                    self.wfile.write(body)
                    return
                if path.startswith("/api/tasks/") and path.count("/") == 3:
                    task_id = unquote(path.split("/")[3])
                    task = store.get(task_id)
                    if not task:
                        self.send_error_json("task not found", HTTPStatus.NOT_FOUND)
                        return
                    self.send_json({"task": task_snapshots([task])[0]})
                    return
                if path.startswith("/api/tasks/") and path.endswith("/log"):
                    task_id = unquote(path.split("/")[3])
                    log_path = config.jobs_dir / task_id / "worker.log"
                    self.send_json({"task_id": task_id, "log": log_path.read_text(encoding="utf-8") if log_path.exists() else ""})
                    return
                if path.startswith("/tasks/") and path.endswith("/export.zip"):
                    task_id = unquote(path.strip("/").split("/")[1])
                    task = store.get(task_id)
                    if not task:
                        self.send_error_json("task not found", HTTPStatus.NOT_FOUND)
                        return
                    if task.get("status") != "finished":
                        self.send_error_json("task is not finished", HTTPStatus.CONFLICT)
                        return
                    package = build_report_package(config, task_id)
                    self.send_file(package, "application/zip", package.name)
                    return
                if path.startswith("/tasks/") and path.endswith("/export.html"):
                    task_id = unquote(path.strip("/").split("/")[1])
                    task = store.get(task_id)
                    if not task:
                        self.send_error_json("task not found", HTTPStatus.NOT_FOUND)
                        return
                    if task.get("status") != "finished":
                        self.send_error_json("task is not finished", HTTPStatus.CONFLICT)
                        return
                    report = build_standalone_report(config, task_id)
                    self.send_file(report, "text/html; charset=utf-8", report.name)
                    return
                if path.startswith("/tasks/") and path.endswith("/export-static.zip"):
                    task_id = unquote(path.strip("/").split("/")[1])
                    task = store.get(task_id)
                    if not task:
                        self.send_error_json("task not found", HTTPStatus.NOT_FOUND)
                        return
                    if task.get("status") != "finished":
                        self.send_error_json("task is not finished", HTTPStatus.CONFLICT)
                        return
                    package = build_static_offline_package(config, task_id)
                    self.send_file(package, "application/zip", package.name)
                    return
                report = self.task_from_report_path(path)
                if report:
                    task_id, rest = report
                    task = store.get(task_id)
                    if not task:
                        self.send_error_json("task not found", HTTPStatus.NOT_FOUND)
                        return
                    if task.get("status") != "finished":
                        self.send_error_json("task is not finished", HTTPStatus.CONFLICT)
                        return
                    self.handle_report_get(task_id, rest, parse_qs(parsed.query))
                    return
                self.send_error_json(f"Not found: {path}", HTTPStatus.NOT_FOUND)
            except Exception as exc:  # pragma: no cover - browser-facing error path
                self.send_error_json(str(exc), HTTPStatus.INTERNAL_SERVER_ERROR)

        def handle_report_get(self, task_id: str, rest: str, query: dict[str, list[str]]) -> None:
            db_path = task_db_path(config, task_id)
            if rest in {"/", "/lcms_peak_first_compare.html"}:
                content = task_html_path(config, task_id).read_text(encoding="utf-8")
                self.send_bytes(rewrite_report_html(task_id, content).encode("utf-8"), "text/html; charset=utf-8")
                return
            if rest == "/api/comparisons":
                task = store.get(task_id) or {}
                self.send_json({"default_comparison": task_id, "comparisons": [{"id": task_id, "label": task.get("task_name") or task_id}]})
                return
            if rest == "/api/bootstrap":
                self.send_json(read_bootstrap(db_path))
                return
            if rest == "/api/msms":
                try:
                    self.send_json(read_msms_identifications(db_path))
                except KeyError:
                    self.send_error_json("MS/MS identification has not been run", HTTPStatus.NOT_FOUND)
                return
            if rest == "/api/task-reference":
                stored_structure = _read_optional_artifact(db_path, "task_reference_structure")
                if isinstance(stored_structure, dict) and stored_structure.get("available"):
                    self.send_json(stored_structure)
                    return
                task = store.get(task_id) or {}
                references = dict(task.get("references") or {})
                structure = dict(references.get("structure") or {})
                if not structure.get("provided"):
                    self.send_json({"available": False, "reason": "未提供结构文件或 PDB ID。"})
                    return
                structure_path = str(structure.get("path") or "")
                structure_text = ""
                format_name = "pdb"
                if structure_path:
                    candidate = (task_directory(config, task_id) / structure_path).resolve()
                    task_root = task_directory(config, task_id)
                    if not candidate.is_relative_to(task_root) or not candidate.exists():
                        raise FileNotFoundError("task structure file is missing")
                    structure_text = candidate.read_text(encoding="utf-8", errors="replace")
                    format_name = "cif" if candidate.suffix.lower() in {".cif", ".mmcif"} else "pdb"
                elif structure.get("pdb_id"):
                    structure_text, _ = pdb_structure_text(str(structure.get("pdb_id")))
                    format_name = "cif"
                if not structure_text.strip():
                    self.send_json({"available": False, "reason": "结构内容为空。"})
                    return
                self.send_json({
                    "available": True,
                    "format": format_name,
                    "text": structure_text,
                    "label": structure.get("name") or "任务结构",
                    "source_url": (
                        f"https://www.rcsb.org/structure/{quote(str(structure.get('pdb_id')))}"
                        if structure.get("pdb_id") else ""
                    ),
                    "source_meta": {"kind": "pdb" if structure.get("pdb_id") else "uploaded"},
                })
                return
            if rest == "/api/pdb-structure":
                pdb_id = str(query.get("pdb_id", [""])[0])
                structure_text, _ = pdb_structure_text(pdb_id)
                self.send_bytes(structure_text.encode("utf-8"), "chemical/x-mmcif; charset=utf-8")
                return
            if rest == "/assets/3Dmol-min.js":
                viewer_asset = VENDOR_DIR / "3Dmol-min.js"
                if viewer_asset.exists():
                    self.send_file(viewer_asset, "text/javascript; charset=utf-8")
                else:
                    self.send_error_json("3D structure viewer asset is not installed", HTTPStatus.NOT_FOUND)
                return
            if rest == "/api/features":
                self.send_json({"features": read_features(db_path)})
                return
            if rest == "/api/xic":
                peak_id = str(query.get("peak_id", [""])[0])
                mz = float(query.get("mz", ["nan"])[0])
                full_run = str(query.get("full", ["0"])[0]).lower() in {"1", "true", "yes", "full"}
                self.send_json(xic_payload(db_path, peak_id, mz, full_run))
                return
            self.send_error_json(f"Report path not found: {rest}", HTTPStatus.NOT_FOUND)

        def do_POST(self) -> None:  # noqa: N802
            parsed = urlparse(self.path)
            path = parsed.path
            try:
                if path == "/mcp":
                    self.handle_mcp_post()
                    return
                if path == "/api/tasks":
                    self.handle_create_task()
                    return
                action_match = re.fullmatch(r"/api/tasks/([^/]+)/(cancel|retry)", path)
                if action_match:
                    self.handle_task_action(unquote(action_match.group(1)), action_match.group(2))
                    return
                report = self.task_from_report_path(path)
                if report and report[1] == "/api/task-reference/pdb":
                    self.handle_report_pdb_save(report[0])
                    return
                if report and report[1] == "/api/task-reference/upload":
                    self.handle_report_structure_upload(report[0])
                    return
                if report and report[1] == "/api/features":
                    task_id = report[0]
                    db_path = task_db_path(config, task_id)
                    length = int(self.headers.get("Content-Length") or "0")
                    body = self.rfile.read(length)
                    payload = json.loads(body.decode("utf-8")) if body else {}
                    features = payload.get("features", [])
                    if not isinstance(features, list):
                        raise ValueError("features must be a list")
                    count = replace_features(db_path, features)
                    self.send_json({"saved_count": count, "features": read_features(db_path)})
                    return
                self.send_error_json(f"Not found: {path}", HTTPStatus.NOT_FOUND)
            except Exception as exc:
                self.send_error_json(str(exc), HTTPStatus.BAD_REQUEST)

        def handle_mcp_post(self) -> None:
            length = int(self.headers.get("Content-Length") or "0")
            if length <= 0 or length > 2 * 1024 * 1024:
                raise ValueError("MCP request is empty or exceeds the 2 MB limit")
            body = self.rfile.read(length)
            payload = mcp.handle_json(body)
            if payload is None or payload == []:
                self.send_response(HTTPStatus.ACCEPTED)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            self.send_bytes(
                json.dumps(payload, ensure_ascii=False).encode("utf-8"),
                "application/json; charset=utf-8",
            )

        def handle_report_structure_upload(self, task_id: str) -> None:
            task = store.get(task_id)
            if not task:
                self.send_error_json("task not found", HTTPStatus.NOT_FOUND)
                return
            if task.get("status") != "finished":
                self.send_error_json("task is not finished", HTTPStatus.CONFLICT)
                return
            db_path = task_db_path(config, task_id)
            if not db_path.exists():
                self.send_error_json("task report database is missing", HTTPStatus.NOT_FOUND)
                return
            length = int(self.headers.get("Content-Length") or "0")
            if length <= 0 or length > MAX_STRUCTURE_BYTES + 2 * 1024 * 1024:
                raise ValueError("structure upload is empty or exceeds the 40 MB limit")
            content_type = self.headers.get("Content-Type", "")
            staging_dir = Path(tempfile.mkdtemp(prefix="structure_", dir=str(config.root)))
            try:
                fields, files = parse_multipart_stream(content_type, self.rfile, staging_dir)
                uploads = [item for item in files if str(item.get("field_name") or "structure_file") == "structure_file"]
                if len(uploads) != 1:
                    raise ValueError("exactly one PDB or mmCIF structure file is required")
                source = Path(str(uploads[0]["staged_path"]))
                suffix = source.suffix.lower()
                if suffix not in STRUCTURE_SUFFIXES:
                    raise ValueError(f"unsupported structure file type: {source.name}")
                validate_structure_file(source)
                task_dir = task_directory(config, task_id)
                references_dir = task_dir / "references"
                references_dir.mkdir(parents=True, exist_ok=True)
                target = references_dir / f"structure{suffix}"
                shutil.copy2(source, target)
                validate_structure_file(target)
                structure_name = str(uploads[0].get("original_name") or target.name)
                structure_path = str(target.relative_to(task_dir))
                references = dict(task.get("references") or {})
                sequence = dict(references.get("sequence") or {})
                if not sequence.get("provided"):
                    raise ValueError("该任务没有 FASTA 序列，不能进行蛋白质结构映射")
                structure = {
                    "provided": True,
                    "source": "report_upload",
                    "name": structure_name,
                    "pdb_id": "",
                    "path": structure_path.replace("\\", "/"),
                }
                references["structure"] = structure
                reference_files = [item for item in task.get("reference_files") or [] if item.get("kind") != "structure"]
                reference_files.append({
                    "kind": "structure",
                    "original_name": structure_name,
                    "stored_name": target.name,
                    "size_bytes": target.stat().st_size,
                })
                (references_dir / "manifest.json").write_text(
                    json.dumps({"sequence": sequence, "structure": structure}, ensure_ascii=False, indent=2),
                    encoding="utf-8",
                )
                structure_payload = structure_payload_from_file(target, structure_name, {"kind": "uploaded"})
                structure_payload["path"] = structure_path
                upsert_sqlite_artifact(db_path, "task_reference_structure", structure_payload)
                update_bootstrap_metadata(db_path, {
                    "references": references,
                    "structure_mapping": {
                        "status": "available",
                        "reason": "已从报告页面上传结构文件，可进行蛋白质结构映射。",
                        "source": "report_upload",
                        "path": structure_path,
                    },
                })
                updated = store.update(task_id, references=references, reference_files=reference_files)
                self.send_json({
                    "task": task_snapshots([updated])[0],
                    "structure": structure_payload,
                })
            finally:
                if staging_dir.exists():
                    shutil.rmtree(staging_dir, ignore_errors=True)

        def handle_report_pdb_save(self, task_id: str) -> None:
            """Persist a PDB ID entered from an already-generated report."""
            task = store.get(task_id)
            if not task:
                self.send_error_json("task not found", HTTPStatus.NOT_FOUND)
                return
            if task.get("status") != "finished":
                self.send_error_json("task is not finished", HTTPStatus.CONFLICT)
                return
            db_path = task_db_path(config, task_id)
            if not db_path.exists():
                self.send_error_json("task report database is missing", HTTPStatus.NOT_FOUND)
                return

            length = int(self.headers.get("Content-Length") or "0")
            if length <= 0 or length > 64 * 1024:
                raise ValueError("PDB save request is empty or too large")
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            pdb_id = str(payload.get("pdb_id") or "").strip().upper()
            if not re.fullmatch(r"[A-Z0-9]{4,12}", pdb_id):
                raise ValueError("PDB ID must contain 4-12 letters or digits")

            references = dict(task.get("references") or {})
            sequence = dict(references.get("sequence") or {})
            if not sequence.get("provided"):
                raise ValueError("该任务没有 FASTA 序列，不能进行蛋白质结构映射")

            structure_text, _ = pdb_structure_text(pdb_id)
            task_dir = task_directory(config, task_id)
            references_dir = task_dir / "references"
            references_dir.mkdir(parents=True, exist_ok=True)
            target = references_dir / "structure.cif"
            target.write_text(structure_text, encoding="utf-8")
            validate_structure_file(target)

            structure_name = str(payload.get("name") or f"PDB {pdb_id}")[:240]
            structure_path = str(target.relative_to(task_dir)).replace("\\", "/")
            structure = {
                "provided": True,
                "source": "pdb_id",
                "name": structure_name,
                "pdb_id": pdb_id,
                "path": structure_path,
            }
            references["structure"] = structure
            reference_files = [
                item for item in task.get("reference_files") or []
                if item.get("kind") != "structure"
            ]
            reference_files.append({
                "kind": "structure",
                "original_name": structure_name,
                "stored_name": target.name,
                "size_bytes": target.stat().st_size,
            })
            (references_dir / "manifest.json").write_text(
                json.dumps({"sequence": sequence, "structure": structure}, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )

            structure_payload = structure_payload_from_file(
                target,
                structure_name,
                {"kind": "pdb", "pdb_id": pdb_id},
            )
            structure_payload["path"] = structure_path
            structure_payload["source_url"] = f"https://www.rcsb.org/structure/{quote(pdb_id)}"
            upsert_sqlite_artifact(db_path, "task_reference_structure", structure_payload)
            update_bootstrap_metadata(db_path, {
                "references": references,
                "structure_mapping": {
                    "status": "available",
                    "reason": "已从报告页面输入 PDB ID 并保存结构，可进行蛋白质结构映射。",
                    "source": "pdb_id",
                    "pdb_id": pdb_id,
                    "path": structure_path,
                },
            })
            updated = store.update(task_id, references=references, reference_files=reference_files)
            self.send_json({
                "task": task_snapshots([updated])[0],
                "references": references,
                "structure": structure_payload,
            })

        def handle_create_task(self) -> None:
            length = int(self.headers.get("Content-Length") or "0")
            if length <= 0:
                raise ValueError("empty request")
            if length > config.max_upload_bytes:
                raise ValueError(f"upload is too large; maximum is {config.max_upload_bytes // (1024 * 1024 * 1024)} GB")
            content_type = self.headers.get("Content-Type", "")
            staging_dir = Path(tempfile.mkdtemp(prefix="upload_", dir=str(config.root)))
            task_dir: Path | None = None
            task_persisted = False
            try:
                fields, files = parse_multipart_stream(content_type, self.rfile, staging_dir)
                spectrum_uploads = [item for item in files if str(item.get("field_name") or "files") == "files"]
                if len(spectrum_uploads) < 2:
                    raise ValueError("at least two RAW/mzML files are required")
                sequence_uploads = [item for item in files if str(item.get("field_name") or "") == "sequence_file"]
                structure_uploads = [item for item in files if str(item.get("field_name") or "") == "structure_file"]
                if len(sequence_uploads) > 1 or len(structure_uploads) > 1:
                    raise ValueError("only one FASTA and one structure file can be supplied")
                sequence_text = str(fields.get("sequence_text") or "").strip()
                pdb_id = str(fields.get("pdb_id") or "").strip().upper()
                if sequence_uploads and sequence_text:
                    raise ValueError("请只使用 FASTA 文件或序列粘贴框中的一种方式")
                if structure_uploads and pdb_id:
                    raise ValueError("请只使用结构文件或 PDB ID 中的一种方式")
                if pdb_id and not re.fullmatch(r"[A-Z0-9]{4,12}", pdb_id):
                    raise ValueError("PDB ID 必须为 4-12 位字母或数字")
                task_name = fields.get("task_name") or f"LCMS task {int(time.time())}"
                try:
                    submitted_sample_names = json.loads(str(fields.get("sample_names") or "[]"))
                except json.JSONDecodeError as error:
                    raise ValueError("样品名称数据格式无效") from error
                if not isinstance(submitted_sample_names, list) or len(submitted_sample_names) != len(spectrum_uploads):
                    raise ValueError("请为每个上传文件填写样品名称")
                sample_names = [sample_name_value(value, str(item["original_name"])) for value, item in zip(submitted_sample_names, spectrum_uploads)]
                if len({name.casefold() for name in sample_names}) != len(sample_names):
                    raise ValueError("样品名称不能重复，请分别填写唯一名称")
                base_id = safe_id(task_name)
                task_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{base_id}_{uuid.uuid4().hex[:8]}"
                task_dir = config.jobs_dir / task_id
                raw_dir = task_dir / "raw"
                raw_dir.mkdir(parents=True, exist_ok=False)
                references_dir = task_dir / "references"
                references_dir.mkdir(parents=True, exist_ok=True)
                saved_files: list[dict[str, object]] = []
                reference_files: list[dict[str, object]] = []
                for file_index, file_info in enumerate(spectrum_uploads):
                    original_name = str(file_info["original_name"])
                    suffix = Path(original_name).suffix.lower()
                    if suffix not in SPECTRUM_SUFFIXES:
                        raise ValueError(f"unsupported file type: {original_name}")
                    staged_path = Path(str(file_info["staged_path"]))
                    size_bytes = staged_path.stat().st_size
                    if size_bytes <= 0:
                        raise ValueError(f"empty upload is not allowed: {original_name}")
                    target = unique_path(raw_dir, original_name)
                    shutil.move(str(staged_path), str(target))
                    saved_files.append(
                        {
                            "original_name": original_name,
                            "sample_name": sample_names[file_index],
                            "stored_name": target.name,
                            "size_bytes": size_bytes,
                        }
                    )
                references: dict[str, object] = {
                    "sequence": {"provided": False, "source": "none"},
                    "structure": {"provided": False, "source": "none"},
                }
                if sequence_uploads or sequence_text:
                    sequence_target = references_dir / "sequence.fasta"
                    if sequence_uploads:
                        source = Path(str(sequence_uploads[0]["staged_path"]))
                        suffix = source.suffix.lower()
                        if suffix not in FASTA_SUFFIXES:
                            raise ValueError(f"unsupported FASTA file type: {source.name}")
                        shutil.copy2(source, sequence_target)
                        sequence_source = "uploaded_file"
                        sequence_name = str(sequence_uploads[0]["original_name"])
                    else:
                        sequence_target.write_text(sequence_text + "\n", encoding="utf-8")
                        sequence_source = "pasted_text"
                        sequence_name = "pasted_sequence.fasta"
                    sequence_lengths = validate_fasta_file(sequence_target)
                    sequence_size = sequence_target.stat().st_size
                    reference_files.append({
                        "kind": "sequence",
                        "original_name": sequence_name,
                        "stored_name": sequence_target.name,
                        "size_bytes": sequence_size,
                    })
                    references["sequence"] = {
                        "provided": True,
                        "source": sequence_source,
                        "name": sequence_name,
                        "path": "references/sequence.fasta",
                        "chain_lengths": sequence_lengths,
                    }
                if structure_uploads or pdb_id:
                    if structure_uploads:
                        source = Path(str(structure_uploads[0]["staged_path"]))
                        suffix = source.suffix.lower()
                        if suffix not in STRUCTURE_SUFFIXES:
                            raise ValueError(f"unsupported structure file type: {source.name}")
                        structure_target = references_dir / f"structure{suffix if suffix else '.pdb'}"
                        shutil.copy2(source, structure_target)
                        validate_structure_file(structure_target)
                        structure_source = "uploaded_file"
                        structure_name = str(structure_uploads[0]["original_name"])
                        structure_path = str(structure_target.relative_to(task_dir))
                        structure_size = structure_target.stat().st_size
                    else:
                        structure_source = "pdb_id"
                        structure_name = f"PDB {pdb_id}"
                        structure_path = ""
                        structure_size = 0
                    references["structure"] = {
                        "provided": True,
                        "source": structure_source,
                        "name": structure_name,
                        "pdb_id": pdb_id if pdb_id else "",
                        "path": structure_path,
                    }
                    if structure_size:
                        reference_files.append({
                            "kind": "structure",
                            "original_name": structure_name,
                            "stored_name": Path(structure_path).name,
                            "size_bytes": structure_size,
                        })
                manifest = {
                    "sequence": references["sequence"],
                    "structure": references["structure"],
                }
                (references_dir / "manifest.json").write_text(
                    json.dumps(manifest, ensure_ascii=False, indent=2),
                    encoding="utf-8",
                )
                task = {
                    "task_id": task_id,
                    "task_name": task_name,
                    "status": "waiting",
                    "stage": "waiting",
                    "created_at": now_iso(),
                    "updated_at": now_iso(),
                    "started_at": "",
                    "finished_at": "",
                    "progress": 0,
                    "retry_count": 0,
                    "error": "",
                    "files": saved_files,
                    "reference_files": reference_files,
                    "references": references,
                    "params": {
                        "reference_sample": "",
                        "top_n_peaks": DEFAULT_TOP_N_TIC_PEAKS,
                        "top_n_mz": DEFAULT_TOP_N_MZ,
                        "top_n_changed_mz": DEFAULT_TOP_N_CHANGED_MZ,
                        "mz_tolerance_da": 0.16,
                        "mz_tolerance_ppm": 10.0,
                        "mz_tolerance_mode": "da",
                        "max_spectrum_points_per_scan": DEFAULT_MAX_SPECTRUM_POINTS_PER_SCAN,
                        "max_peaks_per_scan": DEFAULT_MAX_PEAKS_PER_SCAN,
                    },
                    "job_dir": str(task_dir),
                    "report_url": "",
                }
                store.create(task)
                task_persisted = True
                self.send_json({"task": task_snapshots([task])[0]}, HTTPStatus.CREATED)
            finally:
                if staging_dir.exists():
                    shutil.rmtree(staging_dir, ignore_errors=True)
                if task_dir is not None and not task_persisted and task_dir.exists():
                    jobs_root = config.jobs_dir.resolve()
                    candidate = task_dir.resolve()
                    if candidate.is_relative_to(jobs_root) and candidate != jobs_root:
                        shutil.rmtree(candidate, ignore_errors=True)

        def handle_task_action(self, task_id: str, action: str) -> None:
            task = store.get(task_id)
            if not task:
                self.send_error_json("task not found", HTTPStatus.NOT_FOUND)
                return
            if action == "cancel":
                updated = store.request_cancel(task_id)
                self.send_json({"task": task_snapshots([updated])[0]})
                return
            if action == "retry":
                if task.get("status") != "failed":
                    raise ValueError("only failed tasks can be retried")
                clear_task_outputs(config, task_id)
                updated = store.update(
                    task_id,
                    status="waiting",
                    stage="waiting",
                    progress=0,
                    started_at="",
                    finished_at="",
                    error="",
                    report_url="",
                    output_dir="",
                    retry_count=int(task.get("retry_count") or 0) + 1,
                )
                self.send_json({"task": task_snapshots([updated])[0]})
                return
            raise ValueError(f"unsupported action: {action}")

    return Handler


def default_parser_path() -> Path:
    """Prefer the self-contained converter shipped with the migration package."""
    bundled = WORKSPACE / "lcms_department_platform" / "tools" / "ThermoRawFileParser" / "ThermoRawFileParser.exe"
    legacy = WORKSPACE / ".local-tools" / "ThermoRawFileParser" / "current" / "ThermoRawFileParser.exe"
    for candidate in (bundled, legacy):
        if candidate.exists():
            return candidate
    return bundled


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Serve the LC-MS department upload and queue portal.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8770)
    parser.add_argument("--root", default=str(WORKSPACE / "lcms_department_platform"))
    parser.add_argument(
        "--jobs-dir",
        default="",
        help="Optional task directory root; defaults to <root>/jobs.",
    )
    parser.add_argument("--parser-path", default=str(default_parser_path()))
    parser.add_argument("--python", default=sys.executable)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = PortalConfig(
        host=args.host,
        port=args.port,
        root=Path(args.root).resolve(),
        parser_path=Path(args.parser_path).resolve(),
        python_executable=str(args.python),
        jobs_dir_override=Path(args.jobs_dir).resolve() if args.jobs_dir else None,
    )
    config.state_dir.mkdir(parents=True, exist_ok=True)
    config.jobs_dir.mkdir(parents=True, exist_ok=True)
    store = TaskStore(config)
    recovered = store.recover_interrupted()
    if recovered:
        print(f"Requeued {recovered} task(s) left running by a previous server process.")
    worker = TaskWorker(store)
    worker.start()
    server = ThreadingHTTPServer((config.host, config.port), make_handler(config, store))
    print(f"Serving LC-MS department platform at http://{config.host}:{config.port}/")
    print(f"Task root: {config.root}")
    print(f"ThermoRawFileParser: {config.parser_path}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server.")
        worker.stop_requested.set()


if __name__ == "__main__":
    main()
