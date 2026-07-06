#!/usr/bin/env python3
"""Export helpers for LC-MS MVP outputs."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from .lcms_models import CandidateMz, LCMSFeature, LCMSFeatureGroup, LCMSRawFile, XICTrace


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def write_candidates_csv(path: Path, candidates: list[CandidateMz]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank", "candidate_mz", "total_intensity", "max_intensity", "scan_count", "rt_start", "rt_end"])
        for candidate in candidates:
            writer.writerow([
                candidate.rank,
                f"{candidate.candidate_mz:.6f}",
                f"{candidate.total_intensity:.6f}",
                f"{candidate.max_intensity:.6f}",
                candidate.scan_count,
                f"{candidate.rt_range[0]:.5f}",
                f"{candidate.rt_range[1]:.5f}",
            ])


def write_feature_matrix_csv(path: Path, groups: list[LCMSFeatureGroup], sample_ids: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "feature_group_id",
            "representative_mz",
            "representative_rt",
            "sample_count",
            "missing_sample_count",
            "cv",
            "max_fold_change",
            "difference_type",
            "confidence_score",
            *sample_ids,
        ])
        for group in groups:
            writer.writerow([
                group.feature_group_id,
                f"{group.representative_mz:.6f}",
                f"{group.representative_rt:.5f}",
                group.sample_count,
                group.missing_sample_count,
                "" if group.cv is None else f"{group.cv:.6f}",
                "" if group.max_fold_change is None else f"{group.max_fold_change:.6f}",
                group.difference_type,
                f"{group.confidence_score:.4f}",
                *[f"{group.area_by_sample.get(sample_id, 0.0):.6f}" for sample_id in sample_ids],
            ])


def payload_for_report(
    raw_files: list[LCMSRawFile],
    candidates_by_sample: dict[str, list[CandidateMz]],
    features: list[LCMSFeature],
    groups: list[LCMSFeatureGroup],
    xics_by_group: dict[str, list[XICTrace]],
    sample_ids: list[str],
) -> dict[str, object]:
    return {
        "raw_files": [asdict(item) for item in raw_files],
        "sample_ids": sample_ids,
        "candidates_by_sample": {
            sample_id: [asdict(candidate) for candidate in candidates]
            for sample_id, candidates in candidates_by_sample.items()
        },
        "features": [asdict(item) for item in features],
        "feature_groups": [asdict(item) for item in groups],
        "xics_by_group": {
            group_id: [asdict(trace) for trace in traces]
            for group_id, traces in xics_by_group.items()
        },
    }

