#!/usr/bin/env python3
"""Serve the standalone Peak-first LC-MS compare workbench."""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
import time
from contextlib import closing
from functools import lru_cache
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Serve LC-MS Peak-first Compare V2.")
    parser.add_argument("--output-dir", default="lcms_feature_mvp/outputs_peak_first")
    parser.add_argument(
        "--comparison",
        action="append",
        default=[],
        help="Named comparison as id|label|output_dir. Can be repeated.",
    )
    parser.add_argument("--sqlite-name", default="lcms_peak_first_compare.sqlite")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8768)
    return parser.parse_args()


@lru_cache(maxsize=64)
def read_artifact_cached(db_path_text: str, mtime_ns: int, key: str) -> object:
    del mtime_ns  # part of the cache key so rebuilt SQLite files invalidate cached JSON
    with closing(sqlite3.connect(db_path_text)) as connection:
        row = connection.execute(
            "SELECT payload_json FROM peak_first_artifacts WHERE artifact_key = ?",
            (key,),
        ).fetchone()
    if not row:
        raise KeyError(key)
    return json.loads(str(row[0]))


def read_artifact(db_path: Path, key: str) -> object:
    resolved = db_path.resolve()
    return read_artifact_cached(str(resolved), resolved.stat().st_mtime_ns, key)


def ensure_feature_table(connection: sqlite3.Connection) -> None:
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS saved_lcms_features (
            feature_id TEXT PRIMARY KEY,
            feature_json TEXT NOT NULL,
            created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
            updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
        )
        """
    )


def read_features(db_path: Path) -> list[dict[str, object]]:
    with closing(sqlite3.connect(db_path)) as connection:
        ensure_feature_table(connection)
        rows = connection.execute(
            "SELECT feature_json FROM saved_lcms_features ORDER BY created_at, feature_id"
        ).fetchall()
    features = []
    for row in rows:
        try:
            item = json.loads(str(row[0]))
            if isinstance(item, dict):
                features.append(item)
        except json.JSONDecodeError:
            continue
    return features


def replace_features(db_path: Path, features: list[object]) -> int:
    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with closing(sqlite3.connect(db_path)) as connection:
        ensure_feature_table(connection)
        connection.execute("DELETE FROM saved_lcms_features")
        for index, feature in enumerate(features, start=1):
            if not isinstance(feature, dict):
                continue
            feature_id = str(feature.get("feature_id") or feature.get("region_id") or f"PF_{index:05d}")
            feature["feature_id"] = feature_id
            created = str(feature.get("created_at") or now)
            connection.execute(
                """
                INSERT INTO saved_lcms_features (feature_id, feature_json, created_at, updated_at)
                VALUES (?, ?, ?, ?)
                """,
                (feature_id, json.dumps(feature, ensure_ascii=False, separators=(",", ":")), created, now),
            )
        connection.commit()
    return len([item for item in features if isinstance(item, dict)])


def mz_tolerance(mz: float, params: dict[str, object]) -> float:
    if str(params.get("mz_tolerance_mode") or "da") == "ppm":
        return max(mz * float(params.get("mz_tolerance_ppm") or 20.0) / 1_000_000.0, 1e-9)
    return float(params.get("mz_tolerance_da") or 0.5)


def integrate(rt: list[float], intensity: list[float]) -> float:
    if len(rt) < 2:
        return 0.0
    area = 0.0
    for index in range(1, len(rt)):
        area += (intensity[index - 1] + intensity[index]) * 0.5 * max(rt[index] - rt[index - 1], 0.0)
    return area


def xic_payload(db_path: Path, peak_id: str, target_mz: float, full_run: bool = False) -> dict[str, object]:
    payload = read_artifact(db_path, "bootstrap")
    assert isinstance(payload, dict)
    params = payload.get("params") or {}
    assert isinstance(params, dict)
    peak = next((item for item in payload.get("peak_results", []) if isinstance(item, dict) and item.get("tic_peak_id") == peak_id), None)
    if peak is None:
        raise KeyError(peak_id)
    tolerance = mz_tolerance(target_mz, params)
    feature_group = None
    feature_groups = peak.get("feature_groups") if isinstance(peak.get("feature_groups"), list) else []
    if feature_groups:
        feature_group = min(
            feature_groups,
            key=lambda group: abs(float(group.get("representative_mz") or target_mz) - target_mz) if isinstance(group, dict) else float("inf"),
        )
        if not isinstance(feature_group, dict) or abs(float(feature_group.get("representative_mz") or target_mz) - target_mz) > max(tolerance * 2.0, 1e-9):
            feature_group = None
    feature_bounds_by_sample: dict[str, dict[str, float]] = {}
    if isinstance(feature_group, dict):
        features = feature_group.get("features_by_sample")
        if isinstance(features, dict):
            for sample_id, feature in features.items():
                if not isinstance(feature, dict):
                    continue
                rt_start = feature.get("rt_start")
                rt_end = feature.get("rt_end")
                if rt_start is None or rt_end is None:
                    continue
                feature_bounds_by_sample[str(sample_id)] = {
                    "rt_start": float(rt_start),
                    "rt_end": float(rt_end),
                    "rt_apex": float(feature.get("rt_apex") or feature.get("aligned_rt_apex") or rt_start),
                    "area": float(feature.get("area") or 0.0),
                    "height": float(feature.get("height") or 0.0),
                    "match_status": str(feature.get("match_status") or ""),
                }
    if feature_bounds_by_sample:
        integration_rt_start = min(bounds["rt_start"] for bounds in feature_bounds_by_sample.values())
        integration_rt_end = max(bounds["rt_end"] for bounds in feature_bounds_by_sample.values())
    else:
        integration_rt_start = float(peak.get("rt_start") or 0.0)
        integration_rt_end = float(peak.get("rt_end") or 0.0)
    display_rt_start = integration_rt_start
    display_rt_end = integration_rt_end
    context_min = max(float(params.get("xic_context_min") or 3.0), 0.5)
    context_rt_start = max(0.0, integration_rt_start - context_min)
    context_rt_end = integration_rt_end + context_min
    sample_ids = [str(item) for item in payload.get("sample_ids", [])]
    sample_cut_bounds = peak.get("sample_cut_bounds") if isinstance(peak.get("sample_cut_bounds"), dict) else {}
    xic_by_sample: dict[str, dict[str, list[float]]] = {}
    integration_by_sample: dict[str, dict[str, float]] = {}
    for sample_id in sample_ids:
        spectra = read_artifact(db_path, f"spectra:{sample_id}")
        assert isinstance(spectra, list)
        sample_bounds = sample_cut_bounds.get(sample_id, {}) if isinstance(sample_cut_bounds, dict) else {}
        local_shift = float(sample_bounds.get("local_rt_shift") or 0.0) if isinstance(sample_bounds, dict) else 0.0
        sample_feature_bounds = feature_bounds_by_sample.get(sample_id)
        sample_integration_start = sample_feature_bounds["rt_start"] if sample_feature_bounds else integration_rt_start
        sample_integration_end = sample_feature_bounds["rt_end"] if sample_feature_bounds else integration_rt_end
        rt_values: list[float] = []
        intensities: list[float] = []
        for scan in spectra:
            if not isinstance(scan, dict):
                continue
            rt = float(scan.get("aligned_rt") or scan.get("rt") or 0.0) + local_shift
            if not full_run and (rt < context_rt_start or rt > context_rt_end):
                continue
            total = 0.0
            for mz, intensity in zip(scan.get("mz", []), scan.get("intensity", [])):
                if abs(float(mz) - target_mz) <= tolerance:
                    total += float(intensity)
            rt_values.append(rt)
            intensities.append(total)
        integration_pairs = [
            (rt, intensity)
            for rt, intensity in zip(rt_values, intensities)
            if sample_integration_start <= rt <= sample_integration_end
        ]
        area = integrate([rt for rt, _ in integration_pairs], [intensity for _, intensity in integration_pairs])
        integration_intensities = [intensity for _, intensity in integration_pairs]
        integration_rt_values = [rt for rt, _ in integration_pairs]
        height = max(integration_intensities, default=0.0)
        apex = integration_rt_values[integration_intensities.index(height)] if integration_intensities and height > 0 else None
        xic_by_sample[sample_id] = {"rt": rt_values, "intensity": intensities}
        integration_by_sample[sample_id] = {
            "area": sample_feature_bounds.get("area", area) if sample_feature_bounds else area,
            "height": sample_feature_bounds.get("height", height) if sample_feature_bounds else height,
            "rt_apex": sample_feature_bounds.get("rt_apex", apex) if sample_feature_bounds else apex,
            "rt_start": sample_integration_start,
            "rt_end": sample_integration_end,
            "local_rt_shift": local_shift,
            "match_status": sample_feature_bounds.get("match_status", "") if sample_feature_bounds else "",
        }
    all_rt_values = [
        rt
        for series in xic_by_sample.values()
        for rt in series.get("rt", [])
    ]
    display_rt_start = min(all_rt_values) if full_run and all_rt_values else display_rt_start
    display_rt_end = max(all_rt_values) if full_run and all_rt_values else display_rt_end
    return {
        "peak_id": peak_id,
        "feature_group_id": feature_group.get("feature_group_id") if isinstance(feature_group, dict) else None,
        "full_run": full_run,
        "target_mz": target_mz,
        "mz_tolerance": tolerance,
        "rt_start": display_rt_start,
        "rt_end": display_rt_end,
        "integration_rt_start": integration_rt_start,
        "integration_rt_end": integration_rt_end,
        "context_rt_start": context_rt_start,
        "context_rt_end": context_rt_end,
        "sample_cut_bounds": sample_cut_bounds,
        "feature_bounds_by_sample": feature_bounds_by_sample,
        "xic_by_sample": xic_by_sample,
        "integration_by_sample": integration_by_sample,
    }


def _safe_comparison_id(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value.strip())
    return cleaned or "comparison"


def parse_comparison_spec(spec: str, sqlite_name: str) -> dict[str, object]:
    parts = spec.split("|", 2)
    if len(parts) != 3:
        raise ValueError(f"Comparison must be id|label|output_dir: {spec}")
    comparison_id, label, output = parts
    output_dir = Path(output)
    db_path = output_dir / sqlite_name
    html_path = output_dir / "lcms_peak_first_compare.html"
    if not html_path.exists():
        raise FileNotFoundError(f"Peak-first HTML not found: {html_path}")
    if not db_path.exists():
        raise FileNotFoundError(f"Peak-first SQLite not found: {db_path}")
    return {
        "id": _safe_comparison_id(comparison_id),
        "label": label,
        "output_dir": output_dir,
        "db_path": db_path,
        "html_path": html_path,
    }


def make_handler(comparisons: dict[str, dict[str, object]], default_id: str) -> type[BaseHTTPRequestHandler]:
    def comparison_for(query: dict[str, list[str]]) -> dict[str, object]:
        requested = str(query.get("comparison", [default_id])[0] or default_id)
        return comparisons.get(requested) or comparisons[default_id]

    def db_for(query: dict[str, list[str]]) -> Path:
        return Path(comparison_for(query)["db_path"])

    class Handler(BaseHTTPRequestHandler):
        server_version = "LCMSPeakFirstCompare/0.1"

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

        def do_GET(self) -> None:  # noqa: N802
            parsed = urlparse(self.path)
            path = parsed.path
            query = parse_qs(parsed.query)
            try:
                if path in {"/", "/lcms_peak_first_compare.html"}:
                    html_path = Path(comparisons[default_id]["html_path"])
                    self.send_bytes(html_path.read_bytes(), "text/html; charset=utf-8")
                    return
                if path == "/api/comparisons":
                    self.send_json(
                        {
                            "default_comparison": default_id,
                            "comparisons": [
                                {"id": item["id"], "label": item["label"]}
                                for item in comparisons.values()
                            ],
                        }
                    )
                    return
                if path == "/lcms_xcalibur_workbench.html":
                    sibling = Path(comparison_for(query)["output_dir"]) / "lcms_xcalibur_workbench.html"
                    if sibling.exists():
                        self.send_bytes(sibling.read_bytes(), "text/html; charset=utf-8")
                    else:
                        self.send_json({"error": "Original workbench HTML is not in this output directory"}, HTTPStatus.NOT_FOUND)
                    return
                if path == "/api/bootstrap":
                    self.send_json(read_artifact(db_for(query), "bootstrap"))
                    return
                if path == "/api/features":
                    self.send_json({"features": read_features(db_for(query))})
                    return
                if path == "/api/xic":
                    peak_id = str(query.get("peak_id", [""])[0])
                    mz = float(query.get("mz", ["nan"])[0])
                    if not peak_id or not math.isfinite(mz):
                        self.send_json({"error": "peak_id and numeric mz are required"}, HTTPStatus.BAD_REQUEST)
                        return
                    full_run = str(query.get("full", ["0"])[0]).lower() in {"1", "true", "yes", "full"}
                    self.send_json(xic_payload(db_for(query), peak_id, mz, full_run))
                    return
                self.send_json({"error": f"Not found: {path}"}, HTTPStatus.NOT_FOUND)
            except Exception as exc:  # pragma: no cover - browser-facing error path
                self.send_json({"error": str(exc)}, HTTPStatus.INTERNAL_SERVER_ERROR)

        def do_POST(self) -> None:  # noqa: N802
            parsed = urlparse(self.path)
            if parsed.path != "/api/features":
                self.send_json({"error": "Not found"}, HTTPStatus.NOT_FOUND)
                return
            length = int(self.headers.get("Content-Length") or "0")
            body = self.rfile.read(length)
            try:
                payload = json.loads(body.decode("utf-8")) if body else {}
                features = payload.get("features", [])
                if not isinstance(features, list):
                    raise ValueError("features must be a list")
                query = parse_qs(parsed.query)
                db_path = db_for(query)
                count = replace_features(db_path, features)
                self.send_json({"saved_count": count, "features": read_features(db_path)})
            except Exception as exc:
                self.send_json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)

    return Handler


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    db_path = output_dir / args.sqlite_name
    html_path = output_dir / "lcms_peak_first_compare.html"
    if not html_path.exists():
        raise FileNotFoundError(f"Peak-first HTML not found: {html_path}")
    if not db_path.exists():
        raise FileNotFoundError(f"Peak-first SQLite not found: {db_path}")
    comparison_specs = args.comparison or [f"default|Default|{output_dir}"]
    comparison_list = [parse_comparison_spec(spec, args.sqlite_name) for spec in comparison_specs]
    comparisons = {str(item["id"]): item for item in comparison_list}
    default_id = str(comparison_list[0]["id"])
    server = ThreadingHTTPServer((args.host, args.port), make_handler(comparisons, default_id))
    print(f"Serving LC-MS Peak-first Compare at http://{args.host}:{args.port}/")
    print("Comparisons:")
    for item in comparison_list:
        print(f"  {item['id']}: {item['label']} -> {item['output_dir']}")
    print("Press Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server.")


if __name__ == "__main__":
    main()
