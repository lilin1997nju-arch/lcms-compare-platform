#!/usr/bin/env python3
"""Shared data objects for the LC-MS feature comparison MVP."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


@dataclass
class LCMSRawFile:
    raw_file_id: str
    sample_id: str
    project_id: str
    file_name: str
    file_path: str
    data_format: str
    mz_min: float
    mz_max: float
    rt_min: float
    rt_max: float
    scan_count: int
    created_at: str = field(default_factory=utc_now)
    parser_status: str = "parsed"


@dataclass
class LCMSSpectrumScan:
    scan_id: str
    raw_file_id: str
    sample_id: str
    rt: float
    ms_level: int
    mz_array: list[float]
    intensity_array: list[float]
    tic: float
    base_peak_mz: float | None
    base_peak_intensity: float


@dataclass
class CandidateMz:
    candidate_mz: float
    total_intensity: float
    max_intensity: float
    scan_count: int
    rt_range: tuple[float, float]
    rank: int = 0


@dataclass
class XICTrace:
    sample_id: str
    raw_file_id: str
    target_mz: float
    mz_tolerance: float
    rt_array: list[float]
    intensity_array: list[float]
    max_intensity: float
    tic_ratio: float | None = None


@dataclass
class LCMSFeature:
    feature_id: str
    sample_id: str
    raw_file_id: str
    mz: float
    mz_tolerance: float
    rt_start: float
    rt_apex: float
    rt_end: float
    aligned_rt_apex: float | None
    area: float
    height: float
    signal_to_noise: float
    feature_group_id: str | None = None
    match_status: str = "unmatched"
    annotation_status: str = "unannotated"


@dataclass
class LCMSFeatureGroup:
    feature_group_id: str
    representative_mz: float
    representative_rt: float
    sample_count: int
    missing_sample_count: int
    cv: float | None
    max_fold_change: float | None
    difference_type: str
    confidence_score: float
    sample_presence: dict[str, bool]
    area_by_sample: dict[str, float]

