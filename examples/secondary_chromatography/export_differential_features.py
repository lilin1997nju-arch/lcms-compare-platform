#!/usr/bin/env python3
"""Export ranked Peak-first feature groups for a downstream analysis project."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "FeatureHandoffRecord/v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sqlite", required=True, help="Peak-first SQLite result")
    parser.add_argument("--output", required=True, help="Output JSON path")
    parser.add_argument("--max-similarity", type=float, default=0.85)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--source-result-id", default="")
    return parser.parse_args()


def read_bootstrap(path: Path) -> dict[str, Any]:
    with sqlite3.connect(path) as connection:
        row = connection.execute(
            "SELECT payload_json FROM peak_first_artifacts WHERE artifact_key = ?",
            ("bootstrap",),
        ).fetchone()
    if row is None:
        raise ValueError(f"bootstrap artifact is missing: {path}")
    payload = json.loads(row[0])
    if not isinstance(payload, dict):
        raise ValueError("bootstrap artifact must be a JSON object")
    return payload


def index_feature_groups(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for peak in payload.get("peak_results") or []:
        for group in peak.get("feature_groups") or []:
            group_id = str(group.get("feature_group_id") or "")
            if group_id:
                result[group_id] = group
    return result


def handoff_records(
    payload: dict[str, Any],
    source_result_id: str,
    source_sqlite: Path,
    max_similarity: float,
    limit: int,
) -> list[dict[str, Any]]:
    groups_by_id = index_feature_groups(payload)
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    rows = list(payload.get("global_feature_groups") or [])
    rows.sort(key=lambda row: float(row.get("ranking_score") or 0.0), reverse=True)
    records: list[dict[str, Any]] = []
    for row in rows:
        similarity = float(row.get("similarity_score") or 0.0)
        if similarity > max_similarity:
            continue
        source_group_ids = list(row.get("merged_feature_group_ids") or [])
        if not source_group_ids and row.get("feature_group_id"):
            source_group_ids = [str(row["feature_group_id"])]
        source_group = next(
            (groups_by_id[group_id] for group_id in source_group_ids if group_id in groups_by_id),
            {},
        )
        features_by_sample = dict(source_group.get("features_by_sample") or {})
        area_by_sample = dict(row.get("area_by_sample") or source_group.get("area_by_sample") or {})
        local_shifts = dict(row.get("rt_correction_by_sample") or source_group.get("rt_correction_by_sample") or {})
        sample_measurements: dict[str, dict[str, Any]] = {}
        for sample_id in payload.get("sample_ids") or []:
            feature = dict(features_by_sample.get(sample_id) or {})
            sample_measurements[str(sample_id)] = {
                "area": float(area_by_sample.get(sample_id) or feature.get("area") or 0.0),
                "height": float(feature.get("height") or 0.0),
                "rt_start": feature.get("rt_start"),
                "rt_apex": feature.get("aligned_rt_apex", feature.get("rt_apex")),
                "rt_end": feature.get("rt_end"),
                "global_rt_shift": float(shifts.get(sample_id) or 0.0),
                "local_rt_shift": float(local_shifts.get(sample_id) or 0.0),
                "match_status": feature.get("match_status", "unknown"),
                "gap_filled": bool(feature.get("gap_filled", False)),
            }
        records.append(
            {
                "schema_version": SCHEMA_VERSION,
                "source_result_id": source_result_id,
                "source_sqlite": str(source_sqlite),
                "project_id": payload.get("project_id"),
                "feature_group_id": row.get("feature_group_id"),
                "source_feature_group_ids": source_group_ids,
                "parent_tic_peak_ids": list(
                    row.get("merged_parent_tic_peak_ids")
                    or row.get("source_parent_tic_peak_ids")
                    or [row.get("parent_tic_peak_id")]
                ),
                "representative_mz": row.get("representative_mz"),
                "representative_aligned_rt": row.get("representative_rt"),
                "similarity_score": similarity,
                "difference_score": float(row.get("difference_score") or 0.0),
                "ranking_score": float(row.get("ranking_score") or 0.0),
                "difference_type": row.get("difference_type"),
                "max_fold_change": row.get("max_fold_change"),
                "sample_measurements": sample_measurements,
                "review": {
                    "status": "candidate",
                    "reviewer": None,
                    "reviewed_at": None,
                    "notes": "",
                },
            }
        )
        if len(records) >= max(0, limit):
            break
    return records


def main() -> None:
    args = parse_args()
    sqlite_path = Path(args.sqlite).resolve()
    output_path = Path(args.output).resolve()
    payload = read_bootstrap(sqlite_path)
    result_id = args.source_result_id or sqlite_path.parent.name
    records = handoff_records(
        payload,
        result_id,
        sqlite_path,
        max_similarity=max(0.0, min(1.0, args.max_similarity)),
        limit=args.limit,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps({"schema_version": SCHEMA_VERSION, "features": records}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"Exported {len(records)} feature records to {output_path}")


if __name__ == "__main__":
    main()
