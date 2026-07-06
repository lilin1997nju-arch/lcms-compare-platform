#!/usr/bin/env python3
"""Serve the LC-MS workbench HTML and SQLite-backed API locally."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import mimetypes
import sqlite3
import statistics
import time
from contextlib import closing
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from core.lcms_workbench import (
    cohort_difference_region,
    difference_region,
    heatmap_signal_floor,
    local_max,
    neighborhood_radius,
    normalize_grid,
    similarity_score,
    sort_difference_regions,
    top_spatially_distinct_regions,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Serve the LC-MS workbench from local SQLite storage.")
    parser.add_argument("--output-dir", default="lcms_feature_mvp/outputs_workbench")
    parser.add_argument(
        "--comparison",
        action="append",
        default=[],
        help="Named comparison as id|label|output_dir. Can be repeated.",
    )
    parser.add_argument("--sqlite-name", default="lcms_workbench.sqlite")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    return parser.parse_args()


def _safe_comparison_id(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value.strip())
    return cleaned or "comparison"


def parse_comparison_spec(spec: str, sqlite_name: str) -> dict[str, object]:
    parts = spec.split("|", 2)
    if len(parts) != 3:
        raise ValueError(f"Comparison must be id|label|output_dir: {spec}")
    comparison_id, label, output = parts
    output_dir = Path(output).resolve()
    db_path = output_dir / sqlite_name
    html_path = output_dir / "lcms_xcalibur_workbench.html"
    if not html_path.exists():
        raise FileNotFoundError(f"Workbench HTML not found: {html_path}")
    if not db_path.exists():
        raise FileNotFoundError(f"Workbench SQLite database not found: {db_path}")
    return {
        "id": _safe_comparison_id(comparison_id),
        "label": label,
        "output_dir": output_dir,
        "db_path": db_path,
        "html_path": html_path,
    }


def read_payload_json(db_path: Path) -> str:
    with closing(sqlite3.connect(db_path)) as connection:
        row = connection.execute(
            "SELECT payload_json FROM workbench_artifacts WHERE artifact_key = ?",
            ("payload",),
        ).fetchone()
        if row is not None:
            return str(row[0])
        bootstrap_row = connection.execute(
            "SELECT payload_json FROM workbench_artifacts WHERE artifact_key = ?",
            ("bootstrap",),
        ).fetchone()
        if bootstrap_row is None:
            raise FileNotFoundError(f"bootstrap artifact not found in {db_path}")
        payload = json.loads(str(bootstrap_row[0]))
        spectra_rows = connection.execute(
            "SELECT artifact_key, payload_json FROM workbench_artifacts WHERE artifact_key LIKE 'spectra:%'"
        ).fetchall()
        heatmap_rows = connection.execute(
            "SELECT artifact_key, payload_json FROM workbench_artifacts WHERE artifact_key LIKE 'heatmap:%'"
        ).fetchall()
    spectra_by_key = {key.removeprefix("spectra:"): json.loads(value) for key, value in spectra_rows}
    payload["spectra"] = {
        sample_id: spectra_by_key[sample_id]
        for sample_id in payload.get("sample_ids", [])
        if sample_id in spectra_by_key
    }
    heatmaps_by_key = {key: json.loads(value) for key, value in heatmap_rows}
    payload["heatmaps"] = []
    for item in payload.get("heatmap_catalog", []):
        key = item.get("comparison_key") or item.get("sample_id")
        method = item.get("normalization_method") or "none"
        heatmap = heatmaps_by_key.get(f"heatmap:{key}:{method}")
        if heatmap is not None:
            payload["heatmaps"].append(heatmap)
    return json.dumps(payload, ensure_ascii=False)


def read_artifact_json(db_path: Path, artifact_key: str) -> str:
    with closing(sqlite3.connect(db_path)) as connection:
        row = connection.execute(
            "SELECT payload_json FROM workbench_artifacts WHERE artifact_key = ?",
            (artifact_key,),
        ).fetchone()
    if row is None:
        raise FileNotFoundError(f"artifact not found: {artifact_key}")
    return str(row[0])


def read_heatmap(db_path: Path, key: str, method: str) -> dict[str, object]:
    return json.loads(read_artifact_json(db_path, f"heatmap:{key}:{method}"))


def read_bootstrap(db_path: Path) -> dict[str, object]:
    return json.loads(read_artifact_json(db_path, "bootstrap"))


def ensure_feature_table(connection: sqlite3.Connection) -> None:
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS saved_lcms_features (
            feature_id TEXT PRIMARY KEY,
            feature_json TEXT NOT NULL,
            created_at REAL NOT NULL,
            updated_at REAL NOT NULL
        )
        """
    )


def read_saved_features(db_path: Path) -> list[object]:
    with closing(sqlite3.connect(db_path)) as connection:
        ensure_feature_table(connection)
        rows = connection.execute(
            "SELECT feature_json FROM saved_lcms_features ORDER BY created_at, feature_id"
        ).fetchall()
    return [json.loads(str(row[0])) for row in rows]


def replace_saved_features(db_path: Path, features: list[object]) -> int:
    now = time.time()
    with closing(sqlite3.connect(db_path)) as connection:
        ensure_feature_table(connection)
        existing_created_at = {
            str(feature_id): float(created_at)
            for feature_id, created_at in connection.execute(
                "SELECT feature_id, created_at FROM saved_lcms_features"
            ).fetchall()
        }
        connection.execute("DELETE FROM saved_lcms_features")
        for index, feature in enumerate(features):
            if not isinstance(feature, dict):
                continue
            feature_id = str(feature.get("region_id") or feature.get("feature_id") or f"LCMSF_{index + 1}_{int(now * 1000)}")
            feature["region_id"] = feature_id
            created_at = existing_created_at.get(feature_id, now)
            connection.execute(
                """
                INSERT INTO saved_lcms_features (feature_id, feature_json, created_at, updated_at)
                VALUES (?, ?, ?, ?)
                """,
                (feature_id, json.dumps(feature, ensure_ascii=False, separators=(",", ":")), created_at, now),
            )
        connection.commit()
    return len([feature for feature in features if isinstance(feature, dict)])


def feature_matrix_columns(sample_ids: list[str]) -> list[str]:
    fixed = [
        "feature_id",
        "comparison",
        "comparison_type",
        "normalization",
        "aligned_rt",
        "rt_start",
        "rt_end",
        "mz",
        "mz_start",
        "mz_end",
        "similarity_score",
        "difference_score",
        "fold_change",
        "difference_type",
        "sample_presence",
    ]
    sample_columns: list[str] = []
    for sample_id in sample_ids:
        sample_columns.extend(
            [
                f"{sample_id} raw",
                f"{sample_id} normalized",
                f"{sample_id} present",
                f"{sample_id} scan_id",
                f"{sample_id} raw_rt",
                f"{sample_id} aligned_rt",
            ]
        )
    return fixed + sample_columns


def build_feature_matrix(features: list[object], sample_ids: list[str]) -> dict[str, object]:
    columns = feature_matrix_columns(sample_ids)
    rows: list[dict[str, object]] = []
    for feature in features:
        if not isinstance(feature, dict):
            continue
        sample_intensities = feature.get("sample_intensities") or {}
        source_scans = feature.get("source_scans") or {}
        if not isinstance(sample_intensities, dict):
            sample_intensities = {}
        if not isinstance(source_scans, dict):
            source_scans = {}
        def field(name: str) -> object:
            value = feature.get(name)
            return "" if value is None else value

        row: dict[str, object] = {
            "feature_id": feature.get("region_id") or feature.get("feature_id") or "",
            "comparison": field("comparison_label"),
            "comparison_type": field("comparison_type"),
            "normalization": field("normalization_method"),
            "aligned_rt": field("aligned_rt"),
            "rt_start": field("rt_start"),
            "rt_end": field("rt_end"),
            "mz": field("mz"),
            "mz_start": field("mz_start"),
            "mz_end": field("mz_end"),
            "similarity_score": field("similarity_score"),
            "difference_score": field("difference_score"),
            "fold_change": field("fold_change"),
            "difference_type": field("difference_type"),
            "sample_presence": "|".join(str(item) for item in feature.get("sample_presence") or []),
        }
        for sample_id in sample_ids:
            intensity = sample_intensities.get(sample_id) or {}
            trace = source_scans.get(sample_id) or {}
            if not isinstance(intensity, dict):
                intensity = {}
            if not isinstance(trace, dict):
                trace = {}
            raw = float(intensity.get("raw") or 0.0)
            normalized = float(intensity.get("normalized") or 0.0)
            row[f"{sample_id} raw"] = raw
            row[f"{sample_id} normalized"] = normalized
            row[f"{sample_id} present"] = "present" if raw > 0 else "missing"
            row[f"{sample_id} scan_id"] = trace.get("scan_id") or ""
            row[f"{sample_id} raw_rt"] = "" if trace.get("raw_rt") is None else trace.get("raw_rt")
            row[f"{sample_id} aligned_rt"] = "" if trace.get("aligned_rt") is None else trace.get("aligned_rt")
        rows.append(row)
    return {
        "sample_ids": sample_ids,
        "feature_count": len(rows),
        "columns": columns,
        "rows": rows,
    }


def feature_matrix_csv(matrix: dict[str, object]) -> str:
    columns = [str(item) for item in matrix.get("columns", [])]
    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(columns)
    for row in matrix.get("rows", []):
        if isinstance(row, dict):
            writer.writerow([row.get(column, "") for column in columns])
    return output.getvalue()


def read_feature_matrix(db_path: Path) -> dict[str, object]:
    bootstrap = read_bootstrap(db_path)
    sample_ids = [str(item) for item in bootstrap.get("sample_ids", [])]
    return build_feature_matrix(read_saved_features(db_path), sample_ids)


def optional_float(query: dict[str, list[str]], name: str) -> float | None:
    values = query.get(name) or []
    if not values or values[0] == "":
        return None
    return float(values[0])


def _clamped_index_range(
    minimum: float,
    maximum: float,
    count: int,
    requested_min: float | None,
    requested_max: float | None,
) -> tuple[int, int]:
    if count <= 0 or maximum <= minimum:
        return 0, 0
    step = (maximum - minimum) / count
    low = minimum if requested_min is None else max(minimum, min(maximum, requested_min))
    high = maximum if requested_max is None else max(minimum, min(maximum, requested_max))
    if high < low:
        low, high = high, low
    start = max(0, min(count - 1, int((low - minimum) // step)))
    end = max(start + 1, min(count, int((high - minimum + step - 1e-12) // step) + 1))
    return start, end


def _slice_grid(grid: object, mz_start: int, mz_end: int, rt_start: int, rt_end: int) -> object:
    if not isinstance(grid, list):
        return grid
    return [
        row[rt_start:rt_end] if isinstance(row, list) else row
        for row in grid[mz_start:mz_end]
    ]


def _slice_grid_map(grids: object, mz_start: int, mz_end: int, rt_start: int, rt_end: int) -> object:
    if not isinstance(grids, dict):
        return grids
    return {
        sample_id: _slice_grid(grid, mz_start, mz_end, rt_start, rt_end)
        for sample_id, grid in grids.items()
    }


def heatmap_window(
    heatmap: dict[str, object],
    rt_min: float | None = None,
    rt_max: float | None = None,
    mz_min: float | None = None,
    mz_max: float | None = None,
) -> dict[str, object]:
    original_rt_min = float(heatmap.get("rt_min", 0.0))
    original_rt_max = float(heatmap.get("rt_max", original_rt_min))
    original_mz_min = float(heatmap.get("mz_min", 0.0))
    original_mz_max = float(heatmap.get("mz_max", original_mz_min))
    original_rt_count = int(heatmap.get("rt_bin_count", 0))
    original_mz_count = int(heatmap.get("mz_bin_count", 0))
    rt_start, rt_end = _clamped_index_range(original_rt_min, original_rt_max, original_rt_count, rt_min, rt_max)
    mz_start, mz_end = _clamped_index_range(original_mz_min, original_mz_max, original_mz_count, mz_min, mz_max)
    rt_step = (original_rt_max - original_rt_min) / original_rt_count if original_rt_count else 0.0
    mz_step = (original_mz_max - original_mz_min) / original_mz_count if original_mz_count else 0.0
    window = dict(heatmap)
    window["rt_min"] = original_rt_min + rt_start * rt_step
    window["rt_max"] = original_rt_min + rt_end * rt_step
    window["mz_min"] = original_mz_min + mz_start * mz_step
    window["mz_max"] = original_mz_min + mz_end * mz_step
    window["rt_bin_count"] = rt_end - rt_start
    window["mz_bin_count"] = mz_end - mz_start
    window["source_rt_min"] = original_rt_min
    window["source_rt_max"] = original_rt_max
    window["source_mz_min"] = original_mz_min
    window["source_mz_max"] = original_mz_max
    window["source_rt_bin_offset"] = rt_start
    window["source_mz_bin_offset"] = mz_start
    window["scores"] = _slice_grid(heatmap.get("scores"), mz_start, mz_end, rt_start, rt_end)
    for key in [
        "reference_intensity_grid",
        "sample_intensity_grid",
        "reference_raw_intensity_grid",
        "sample_raw_intensity_grid",
    ]:
        if key in heatmap:
            window[key] = _slice_grid(heatmap.get(key), mz_start, mz_end, rt_start, rt_end)
    for key in ["sample_intensity_grids", "sample_raw_intensity_grids"]:
        if key in heatmap:
            window[key] = _slice_grid_map(heatmap.get(key), mz_start, mz_end, rt_start, rt_end)
    low_points = []
    for point in heatmap.get("low_similarity_points", []) or []:
        if not isinstance(point, dict):
            continue
        rt = float(point.get("rt", 0.0))
        mz = float(point.get("mz", 0.0))
        if window["rt_min"] <= rt <= window["rt_max"] and window["mz_min"] <= mz <= window["mz_max"]:
            low_points.append(point)
    window["low_similarity_points"] = low_points
    return window


def _window_bin_count(span: float, target_width: float, minimum: int, maximum: int) -> int:
    if span <= 0:
        return 1
    return max(minimum, min(maximum, int(math.ceil(span / max(target_width, 1e-9)))))


def _spectra_intensity_grid(
    spectra: list[dict[str, object]],
    rt_bins: list[float],
    mz_bins: list[float],
) -> list[list[float]]:
    rt_count = max(0, len(rt_bins) - 1)
    mz_count = max(0, len(mz_bins) - 1)
    grid = [[0.0 for _ in range(rt_count)] for _ in range(mz_count)]
    if rt_count == 0 or mz_count == 0:
        return grid
    rt_min = rt_bins[0]
    mz_min = mz_bins[0]
    rt_step = rt_bins[1] - rt_bins[0]
    mz_step = mz_bins[1] - mz_bins[0]
    for scan in spectra:
        aligned_rt = float(scan.get("aligned_rt") or scan.get("rt") or 0.0)
        rt_index = int((aligned_rt - rt_min) / rt_step)
        if rt_index < 0 or rt_index >= rt_count:
            continue
        mz_values = scan.get("mz") or []
        intensity_values = scan.get("intensity") or []
        if not isinstance(mz_values, list) or not isinstance(intensity_values, list):
            continue
        for mz, intensity in zip(mz_values, intensity_values):
            mz_value = float(mz)
            mz_index = int((mz_value - mz_min) / mz_step)
            if 0 <= mz_index < mz_count:
                grid[mz_index][rt_index] += float(intensity)
    return grid


def _region_bounds(
    heatmap: dict[str, object],
    rt_min: float | None,
    rt_max: float | None,
    mz_min: float | None,
    mz_max: float | None,
) -> tuple[float, float, float, float]:
    original_rt_min = float(heatmap.get("rt_min", 0.0))
    original_rt_max = float(heatmap.get("rt_max", original_rt_min))
    original_mz_min = float(heatmap.get("mz_min", 0.0))
    original_mz_max = float(heatmap.get("mz_max", original_mz_min))
    low_rt = original_rt_min if rt_min is None else max(original_rt_min, min(original_rt_max, rt_min))
    high_rt = original_rt_max if rt_max is None else max(original_rt_min, min(original_rt_max, rt_max))
    low_mz = original_mz_min if mz_min is None else max(original_mz_min, min(original_mz_max, mz_min))
    high_mz = original_mz_max if mz_max is None else max(original_mz_min, min(original_mz_max, mz_max))
    if high_rt < low_rt:
        low_rt, high_rt = high_rt, low_rt
    if high_mz < low_mz:
        low_mz, high_mz = high_mz, low_mz
    if high_rt <= low_rt:
        high_rt = min(original_rt_max, low_rt + max((original_rt_max - original_rt_min) / max(int(heatmap.get("rt_bin_count", 1)), 1), 1e-6))
    if high_mz <= low_mz:
        high_mz = min(original_mz_max, low_mz + max((original_mz_max - original_mz_min) / max(int(heatmap.get("mz_bin_count", 1)), 1), 1e-6))
    return low_rt, high_rt, low_mz, high_mz


def refined_heatmap_window(
    db_path: Path,
    heatmap: dict[str, object],
    rt_min: float | None = None,
    rt_max: float | None = None,
    mz_min: float | None = None,
    mz_max: float | None = None,
) -> dict[str, object]:
    low_rt, high_rt, low_mz, high_mz = _region_bounds(heatmap, rt_min, rt_max, mz_min, mz_max)
    rt_span = high_rt - low_rt
    mz_span = high_mz - low_mz
    rt_bin_count = _window_bin_count(rt_span, 0.05, 20, 600)
    mz_bin_count = _window_bin_count(mz_span, 1.0, 30, 400)
    rt_step = rt_span / rt_bin_count
    mz_step = mz_span / mz_bin_count
    rt_bins = [low_rt + index * rt_step for index in range(rt_bin_count + 1)]
    mz_bins = [low_mz + index * mz_step for index in range(mz_bin_count + 1)]
    sample_ids = [str(item) for item in heatmap.get("sample_ids", [])]
    if not sample_ids:
        sample_ids = [str(heatmap.get("reference_sample")), str(heatmap.get("sample_id"))]
    spectra_by_sample = {
        sample_id: json.loads(read_artifact_json(db_path, f"spectra:{sample_id}"))
        for sample_id in sample_ids
        if sample_id and sample_id != "__cohort__"
    }
    if not spectra_by_sample:
        return heatmap_window(heatmap, rt_min, rt_max, mz_min, mz_max)

    raw_grids = {
        sample_id: _spectra_intensity_grid(spectra, rt_bins, mz_bins)
        for sample_id, spectra in spectra_by_sample.items()
    }
    signal_floors = {
        sample_id: heatmap_signal_floor(grid)
        for sample_id, grid in raw_grids.items()
    }
    method = str(heatmap.get("normalization_method") or "max")
    normalized_by_sample = {
        sample_id: normalize_grid(grid, method)
        for sample_id, grid in raw_grids.items()
    }
    norm_grids = {sample_id: normalized_by_sample[sample_id][0] for sample_id in raw_grids}
    norm_factors = {sample_id: normalized_by_sample[sample_id][1] for sample_id in raw_grids}
    rt_neighborhood_bins, mz_neighborhood_bins = neighborhood_radius(rt_step, mz_step)

    if heatmap.get("comparison_type") == "pair":
        reference_sample = str(heatmap.get("reference_sample") or sample_ids[0])
        sample_id = str(heatmap.get("sample_id") or sample_ids[-1])
        ref_raw_grid = raw_grids[reference_sample]
        sample_raw_grid = raw_grids[sample_id]
        ref_grid = norm_grids[reference_sample]
        sample_grid = norm_grids[sample_id]
        scores: list[list[float | None]] = []
        differences: list[dict[str, object]] = []
        for mz_index in range(mz_bin_count):
            row: list[float | None] = []
            for rt_index in range(rt_bin_count):
                ref_raw = local_max(ref_raw_grid, mz_index, rt_index, mz_neighborhood_bins, rt_neighborhood_bins)
                sample_raw = local_max(sample_raw_grid, mz_index, rt_index, mz_neighborhood_bins, rt_neighborhood_bins)
                if ref_raw < signal_floors.get(reference_sample, 0.0):
                    ref_raw = 0.0
                if sample_raw < signal_floors.get(sample_id, 0.0):
                    sample_raw = 0.0
                if ref_raw <= 0 and sample_raw <= 0:
                    row.append(None)
                    continue
                ref_norm = 0.0 if ref_raw <= 0 else local_max(ref_grid, mz_index, rt_index, mz_neighborhood_bins, rt_neighborhood_bins)
                sample_norm = 0.0 if sample_raw <= 0 else local_max(sample_grid, mz_index, rt_index, mz_neighborhood_bins, rt_neighborhood_bins)
                score = similarity_score(ref_norm, sample_norm)
                row.append(score)
                if score < 0.75:
                    region = difference_region(reference_sample, sample_id, rt_bins, mz_bins, rt_index, mz_index, score, ref_raw, sample_raw, ref_norm, sample_norm)
                    region["rt_index"] = rt_index
                    region["mz_index"] = mz_index
                    differences.append(region)
            scores.append(row)
        sort_difference_regions(differences)
        window = dict(heatmap)
        window.update({
            "rt_min": low_rt,
            "rt_max": high_rt,
            "mz_min": low_mz,
            "mz_max": high_mz,
            "rt_bin_count": rt_bin_count,
            "mz_bin_count": mz_bin_count,
            "rt_neighborhood_bins": rt_neighborhood_bins,
            "mz_neighborhood_bins": mz_neighborhood_bins,
            "similarity_method": "local_neighborhood_max_refined_window",
            "window_source": "api/heatmap-window-refined",
            "window_bin_target": {"rt_min": 0.05, "mz": 1.0},
            "signal_floor_by_sample": {reference_sample: signal_floors[reference_sample], sample_id: signal_floors[sample_id]},
            "reference_normalization_factor": norm_factors[reference_sample],
            "sample_normalization_factor": norm_factors[sample_id],
            "scores": scores,
            "reference_intensity_grid": ref_grid,
            "sample_intensity_grid": sample_grid,
            "reference_raw_intensity_grid": ref_raw_grid,
            "sample_raw_intensity_grid": sample_raw_grid,
            "low_similarity_points": top_spatially_distinct_regions(differences, rt_neighborhood_bins, mz_neighborhood_bins, 80),
        })
        return window

    normalized_sample_ids = list(raw_grids)
    scores = []
    differences = []
    for mz_index in range(mz_bin_count):
        row = []
        for rt_index in range(rt_bin_count):
            raw_values = {}
            for sample_id in normalized_sample_ids:
                raw_value = local_max(raw_grids[sample_id], mz_index, rt_index, mz_neighborhood_bins, rt_neighborhood_bins)
                raw_values[sample_id] = 0.0 if raw_value < signal_floors.get(sample_id, 0.0) else raw_value
            if max(raw_values.values(), default=0.0) <= 0:
                row.append(None)
                continue
            values = [
                0.0 if raw_values[sample_id] <= 0 else local_max(norm_grids[sample_id], mz_index, rt_index, mz_neighborhood_bins, rt_neighborhood_bins)
                for sample_id in normalized_sample_ids
            ]
            mean_value = statistics.mean(values) if values else 0.0
            cv = statistics.pstdev(values) / mean_value if mean_value > 0 else float("inf")
            score = 0.0 if not math.isfinite(cv) else 1.0 / (1.0 + cv)
            row.append(score)
            if score < 0.75:
                norm_values = {sample_id: values[index] for index, sample_id in enumerate(normalized_sample_ids)}
                region = cohort_difference_region(normalized_sample_ids, rt_bins, mz_bins, rt_index, mz_index, score, raw_values, norm_values)
                region["rt_index"] = rt_index
                region["mz_index"] = mz_index
                differences.append(region)
        scores.append(row)
    sort_difference_regions(differences)
    window = dict(heatmap)
    window.update({
        "rt_min": low_rt,
        "rt_max": high_rt,
        "mz_min": low_mz,
        "mz_max": high_mz,
        "rt_bin_count": rt_bin_count,
        "mz_bin_count": mz_bin_count,
        "rt_neighborhood_bins": rt_neighborhood_bins,
        "mz_neighborhood_bins": mz_neighborhood_bins,
        "similarity_method": "local_neighborhood_max_refined_window",
        "window_source": "api/heatmap-window-refined",
        "window_bin_target": {"rt_min": 0.05, "mz": 1.0},
        "normalization_factors_by_sample": norm_factors,
        "signal_floor_by_sample": signal_floors,
        "scores": scores,
        "sample_intensity_grids": norm_grids,
        "sample_raw_intensity_grids": raw_grids,
        "low_similarity_points": top_spatially_distinct_regions(differences, rt_neighborhood_bins, mz_neighborhood_bins, 80),
    })
    return window


def make_handler(comparisons: dict[str, dict[str, object]], default_id: str) -> type[BaseHTTPRequestHandler]:
    def comparison_for(query: dict[str, list[str]]) -> dict[str, object]:
        requested = str(query.get("comparison", [default_id])[0] or default_id)
        return comparisons.get(requested) or comparisons[default_id]

    def db_for(query: dict[str, list[str]]) -> Path:
        return Path(comparison_for(query)["db_path"])

    def output_for(query: dict[str, list[str]]) -> Path:
        return Path(comparison_for(query)["output_dir"])

    class WorkbenchHandler(BaseHTTPRequestHandler):
        server_version = "LCMSWorkbench/0.1"

        def log_message(self, format: str, *args: object) -> None:
            print(f"{self.address_string()} - {format % args}")

        def send_bytes(self, body: bytes, content_type: str, status: HTTPStatus = HTTPStatus.OK) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
            self.end_headers()
            self.wfile.write(body)

        def send_json(self, payload: object, status: HTTPStatus = HTTPStatus.OK) -> None:
            self.send_bytes(
                json.dumps(payload, ensure_ascii=False).encode("utf-8"),
                "application/json; charset=utf-8",
                status,
            )

        def do_GET(self) -> None:
            parsed = urlparse(self.path)
            path = parsed.path
            query = parse_qs(parsed.query)
            try:
                if path in {"/", "/lcms_xcalibur_workbench.html"}:
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
                db_path = db_for(query)
                output_dir = output_for(query)
                if path == "/api/payload":
                    self.send_bytes(read_payload_json(db_path).encode("utf-8"), "application/json; charset=utf-8")
                    return
                if path == "/api/bootstrap":
                    self.send_bytes(read_artifact_json(db_path, "bootstrap").encode("utf-8"), "application/json; charset=utf-8")
                    return
                if path == "/api/spectra":
                    sample = (query.get("sample") or [""])[0]
                    if not sample:
                        self.send_json({"error": "missing sample query parameter"}, HTTPStatus.BAD_REQUEST)
                        return
                    self.send_bytes(
                        read_artifact_json(db_path, f"spectra:{sample}").encode("utf-8"),
                        "application/json; charset=utf-8",
                    )
                    return
                if path == "/api/heatmap":
                    key = (query.get("key") or [""])[0]
                    method = (query.get("method") or [""])[0]
                    if not key or not method:
                        self.send_json({"error": "missing key or method query parameter"}, HTTPStatus.BAD_REQUEST)
                        return
                    self.send_bytes(
                        read_artifact_json(db_path, f"heatmap:{key}:{method}").encode("utf-8"),
                        "application/json; charset=utf-8",
                    )
                    return
                if path == "/api/heatmap-window":
                    key = (query.get("key") or [""])[0]
                    method = (query.get("method") or [""])[0]
                    if not key or not method:
                        self.send_json({"error": "missing key or method query parameter"}, HTTPStatus.BAD_REQUEST)
                        return
                    window = refined_heatmap_window(
                        db_path,
                        read_heatmap(db_path, key, method),
                        rt_min=optional_float(query, "rt_min"),
                        rt_max=optional_float(query, "rt_max"),
                        mz_min=optional_float(query, "mz_min"),
                        mz_max=optional_float(query, "mz_max"),
                    )
                    self.send_json(window)
                    return
                if path == "/api/metadata":
                    with closing(sqlite3.connect(db_path)) as connection:
                        rows = connection.execute(
                            "SELECT metadata_key, metadata_value FROM workbench_metadata ORDER BY metadata_key"
                        ).fetchall()
                    self.send_json({key: value for key, value in rows})
                    return
                if path == "/api/features":
                    self.send_json({"features": read_saved_features(db_path)})
                    return
                if path == "/api/feature-matrix":
                    self.send_json(read_feature_matrix(db_path))
                    return
                if path == "/api/feature-matrix.csv":
                    self.send_bytes(
                        feature_matrix_csv(read_feature_matrix(db_path)).encode("utf-8-sig"),
                        "text/csv; charset=utf-8",
                    )
                    return
                requested = (output_dir / path.lstrip("/")).resolve()
                if output_dir.resolve() in requested.parents and requested.is_file():
                    mime = mimetypes.guess_type(requested.name)[0] or "application/octet-stream"
                    self.send_bytes(requested.read_bytes(), mime)
                    return
                self.send_json({"error": "not found", "path": path}, HTTPStatus.NOT_FOUND)
            except Exception as exc:
                self.send_json({"error": str(exc), "path": path}, HTTPStatus.INTERNAL_SERVER_ERROR)

        def do_POST(self) -> None:
            parsed = urlparse(self.path)
            path = parsed.path
            try:
                if path != "/api/features":
                    self.send_json({"error": "not found", "path": path}, HTTPStatus.NOT_FOUND)
                    return
                content_length = int(self.headers.get("Content-Length") or "0")
                raw_body = self.rfile.read(content_length) if content_length > 0 else b"{}"
                body = json.loads(raw_body.decode("utf-8"))
                features = body.get("features", body if isinstance(body, list) else [])
                if not isinstance(features, list):
                    self.send_json({"error": "features must be a list"}, HTTPStatus.BAD_REQUEST)
                    return
                query = parse_qs(parsed.query)
                db_path = db_for(query)
                count = replace_saved_features(db_path, features)
                self.send_json({"saved_count": count, "features": read_saved_features(db_path)})
            except Exception as exc:
                self.send_json({"error": str(exc), "path": path}, HTTPStatus.INTERNAL_SERVER_ERROR)

        def do_OPTIONS(self) -> None:
            self.send_bytes(b"", "text/plain; charset=utf-8")

    return WorkbenchHandler


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    db_path = output_dir / args.sqlite_name
    html_path = output_dir / "lcms_xcalibur_workbench.html"
    if not html_path.exists():
        raise FileNotFoundError(f"Workbench HTML not found: {html_path}")
    if not db_path.exists():
        raise FileNotFoundError(f"Workbench SQLite database not found: {db_path}")
    comparison_specs = args.comparison or [f"default|Default|{output_dir}"]
    comparison_list = [parse_comparison_spec(spec, args.sqlite_name) for spec in comparison_specs]
    comparisons = {str(item["id"]): item for item in comparison_list}
    default_id = str(comparison_list[0]["id"])
    server = ThreadingHTTPServer((args.host, args.port), make_handler(comparisons, default_id))
    url = f"http://{args.host}:{args.port}/"
    print(f"LC-MS workbench: {url}")
    print("Comparisons:")
    for item in comparison_list:
        print(f"  {item['id']}: {item['label']} -> {item['output_dir']}")
    print("Press Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
