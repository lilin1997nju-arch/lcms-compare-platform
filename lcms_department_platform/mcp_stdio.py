"""STDIO entry point for the local LC-MS MCP server.

The packaged executable lives beside the Electron executable and therefore
shares its default ``data`` directory and optional data-root configuration.
Only JSON-RPC messages are written to stdout; diagnostics go to stderr.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

if not getattr(sys, "frozen", False):
    source_root = Path(__file__).resolve().parents[1]
    for source_path in (source_root, source_root / "lcms_feature_mvp"):
        if str(source_path) not in sys.path:
            sys.path.insert(0, str(source_path))

try:
    from mcp_server import MCPApplication
except ModuleNotFoundError:
    from lcms_department_platform.mcp_server import MCPApplication


CONFIG_FILE_NAME = "lcms-desktop-config.json"


def installation_root() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parents[1]


def configured_data_root(override: str = "") -> Path:
    explicit = str(override or os.environ.get("LCMS_DATA_ROOT") or "").strip()
    if explicit:
        return Path(explicit).expanduser().resolve()
    install_root = installation_root()
    config_path = install_root / CONFIG_FILE_NAME
    try:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
        configured = str(payload.get("dataRoot") or "").strip() if isinstance(payload, dict) else ""
        if configured:
            return Path(configured).expanduser().resolve()
    except (OSError, ValueError, TypeError):
        pass
    return (install_root / "data").resolve()


@dataclass(frozen=True)
class MCPConfig:
    root: Path

    @property
    def jobs_dir(self) -> Path:
        return (self.root / "jobs").resolve()

    @property
    def tasks_path(self) -> Path:
        return self.root / "state" / "tasks.json"


class ReadOnlyTaskStore:
    def __init__(self, config: MCPConfig) -> None:
        self.config = config

    def load(self) -> list[dict[str, Any]]:
        try:
            payload = json.loads(self.config.tasks_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return []
        return [item for item in payload if isinstance(item, dict)] if isinstance(payload, list) else []

    def get(self, task_id: str) -> dict[str, Any] | None:
        return next((task for task in self.load() if task.get("task_id") == task_id), None)


def _error_response(code: int, message: str, request_id: object = None) -> dict[str, Any]:
    return {"jsonrpc": "2.0", "id": request_id, "error": {"code": code, "message": message}}


def serve_stdio(app: MCPApplication, input_stream: BinaryIO, output_stream: BinaryIO) -> None:
    for raw_line in input_stream:
        if not raw_line.strip():
            continue
        try:
            response = app.handle_json(raw_line)
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            response = _error_response(-32700, f"parse error: {exc}")
        except Exception as exc:  # Keep the process alive after one malformed request.
            response = _error_response(-32600, f"invalid request: {exc}")
        if response is None or response == []:
            continue
        output_stream.write(json.dumps(response, ensure_ascii=False, separators=(",", ":")).encode("utf-8") + b"\n")
        output_stream.flush()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LC-MS department platform MCP server (STDIO).")
    parser.add_argument("--data-root", default="", help="Override the shared LC-MS data directory.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = configured_data_root(args.data_root)
    config = MCPConfig(root=root)
    app = MCPApplication(config, ReadOnlyTaskStore(config), audit=lambda message: print(message, file=sys.stderr, flush=True))
    serve_stdio(app, sys.stdin.buffer, sys.stdout.buffer)


if __name__ == "__main__":
    main()
