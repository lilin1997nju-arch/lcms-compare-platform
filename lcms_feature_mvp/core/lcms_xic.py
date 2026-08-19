#!/usr/bin/env python3
"""Candidate m/z screening and XIC extraction."""

from __future__ import annotations

from .lcms_models import CandidateMz, LCMSSpectrumScan, XICTrace


def ppm_window(target_mz: float, mz_tolerance_ppm: float) -> float:
    return abs(target_mz) * mz_tolerance_ppm / 1_000_000.0


def within_ppm(observed_mz: float, target_mz: float, tolerance_ppm: float) -> bool:
    return abs(observed_mz - target_mz) <= ppm_window(target_mz, tolerance_ppm)


def scans_in_rt_range(scans: list[LCMSSpectrumScan], rt_start: float, rt_end: float) -> list[LCMSSpectrumScan]:
    return [scan for scan in scans if rt_start <= scan.rt <= rt_end and scan.ms_level == 1]


def screen_candidate_mz(
    scans: list[LCMSSpectrumScan],
    rt_start: float,
    rt_end: float,
    mz_min: float = 200.0,
    mz_max: float = 2000.0,
    mz_tolerance_ppm: float = 20.0,
    intensity_threshold: float = 1000.0,
    min_scan_count: int = 3,
    top_n_mz: int = 50,
) -> list[CandidateMz]:
    """Cluster centroid points in an RT window into candidate m/z values."""
    window_scans = scans_in_rt_range(scans, rt_start, rt_end)
    points: list[tuple[float, float, str]] = []
    for scan in window_scans:
        for mz, intensity in zip(scan.mz_array, scan.intensity_array):
            if mz < mz_min or mz > mz_max or intensity < intensity_threshold:
                continue
            points.append((mz, intensity, scan.scan_id))
    points.sort(key=lambda item: item[0])

    clusters: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for mz, intensity, scan_id in points:
        if current is None:
            current = {
                "mz_sum": mz * intensity,
                "weight": intensity,
                "total_intensity": intensity,
                "max_intensity": intensity,
                "scan_ids": {scan_id},
            }
            continue
        center = float(current["mz_sum"]) / float(current["weight"])
        if abs(mz - center) <= ppm_window(center, mz_tolerance_ppm):
            current["mz_sum"] = float(current["mz_sum"]) + mz * intensity
            current["weight"] = float(current["weight"]) + intensity
            current["total_intensity"] = float(current["total_intensity"]) + intensity
            current["max_intensity"] = max(float(current["max_intensity"]), intensity)
            scan_ids = current["scan_ids"]
            assert isinstance(scan_ids, set)
            scan_ids.add(scan_id)
        else:
            clusters.append(current)
            current = {
                "mz_sum": mz * intensity,
                "weight": intensity,
                "total_intensity": intensity,
                "max_intensity": intensity,
                "scan_ids": {scan_id},
            }
    if current is not None:
        clusters.append(current)

    candidates: list[CandidateMz] = []
    for cluster in clusters:
        scan_ids = cluster["scan_ids"]
        assert isinstance(scan_ids, set)
        if len(scan_ids) < min_scan_count:
            continue
        candidates.append(
            CandidateMz(
                candidate_mz=float(cluster["mz_sum"]) / float(cluster["weight"]),
                total_intensity=float(cluster["total_intensity"]),
                max_intensity=float(cluster["max_intensity"]),
                scan_count=len(scan_ids),
                rt_range=(rt_start, rt_end),
            )
        )
    candidates.sort(key=lambda item: item.total_intensity, reverse=True)
    for rank, candidate in enumerate(candidates[:top_n_mz], start=1):
        candidate.rank = rank
    return candidates[:top_n_mz]


def extract_xic(
    scans: list[LCMSSpectrumScan],
    target_mz: float,
    mz_tolerance_ppm: float = 20.0,
    mode: str = "sum",
) -> XICTrace:
    """Extract an RT-intensity curve for target_mz across scans."""
    if mode not in {"sum", "max"}:
        raise ValueError("XIC extraction mode must be 'sum' or 'max'")
    tolerance = ppm_window(target_mz, mz_tolerance_ppm)
    rt_array: list[float] = []
    intensity_array: list[float] = []
    tic_sum = 0.0
    xic_sum = 0.0
    sample_id = scans[0].sample_id if scans else ""
    raw_file_id = scans[0].raw_file_id if scans else ""
    for scan in sorted(scans, key=lambda item: item.rt):
        values = [
            intensity
            for mz, intensity in zip(scan.mz_array, scan.intensity_array)
            if abs(mz - target_mz) <= tolerance
        ]
        intensity = sum(values) if mode == "sum" else (max(values) if values else 0.0)
        rt_array.append(scan.rt)
        intensity_array.append(intensity)
        tic_sum += scan.tic
        xic_sum += intensity
    return XICTrace(
        sample_id=sample_id,
        raw_file_id=raw_file_id,
        target_mz=target_mz,
        mz_tolerance=tolerance,
        rt_array=rt_array,
        intensity_array=intensity_array,
        max_intensity=max(intensity_array, default=0.0),
        tic_ratio=(xic_sum / tic_sum) if tic_sum > 0 else None,
    )
