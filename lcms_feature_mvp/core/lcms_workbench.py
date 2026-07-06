#!/usr/bin/env python3
"""Xcalibur-style LC-MS workbench data preparation."""

from __future__ import annotations

import math
import statistics
from dataclasses import asdict

from .lcms_models import LCMSRawFile, LCMSSpectrumScan


def chromatogram_points(scans: list[LCMSSpectrumScan], rt_shift: float = 0.0) -> list[list[float]]:
    return [
        [scan.rt + rt_shift, scan.tic, scan.base_peak_intensity]
        for scan in sorted(scans, key=lambda item: item.rt)
    ]


def quantile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    pos = max(0.0, min(1.0, q)) * (len(ordered) - 1)
    low = int(math.floor(pos))
    high = int(math.ceil(pos))
    if low == high:
        return ordered[low]
    weight = pos - low
    return ordered[low] * (1.0 - weight) + ordered[high] * weight


def smooth(values: list[float], window: int = 7) -> list[float]:
    if window <= 1 or len(values) < 3:
        return list(values)
    if window % 2 == 0:
        window += 1
    radius = window // 2
    padded = [values[0]] * radius + values + [values[-1]] * radius
    return [
        sum(padded[index - radius:index + radius + 1]) / window
        for index in range(radius, radius + len(values))
    ]


def integrate_segment(times: list[float], values: list[float], start: int, end: int, baseline: float) -> float:
    area = 0.0
    for index in range(start, end):
        y0 = max(0.0, values[index] - baseline)
        y1 = max(0.0, values[index + 1] - baseline)
        area += max(0.0, times[index + 1] - times[index]) * (y0 + y1) / 2.0
    return area


def detect_chromatogram_peaks(
    scans: list[LCMSSpectrumScan],
    rt_start: float | None = None,
    rt_end: float | None = None,
    signal: str = "bpc",
) -> list[dict[str, float]]:
    selected = [
        scan for scan in sorted(scans, key=lambda item: item.rt)
        if (rt_start is None or scan.rt >= rt_start) and (rt_end is None or scan.rt <= rt_end)
    ]
    if len(selected) < 3:
        return []
    times = [scan.rt for scan in selected]
    raw_values = [
        scan.base_peak_intensity if signal == "bpc" else scan.tic
        for scan in selected
    ]
    values = smooth(raw_values, 7)
    baseline = quantile(values, 0.10)
    corrected = [max(0.0, value - baseline) for value in values]
    max_signal = max(corrected, default=0.0)
    if max_signal <= 0:
        return []
    noise = max(quantile(corrected, 0.25), max_signal * 0.01, 1.0)
    min_height = max(noise * 4.0, max_signal * 0.08)
    candidates: list[dict[str, float]] = []
    for index in range(1, len(corrected) - 1):
        if corrected[index] < min_height:
            continue
        if corrected[index] < corrected[index - 1] or corrected[index] <= corrected[index + 1]:
            continue
        foot = max(noise * 1.8, corrected[index] * 0.025)
        left = index
        while left > 0 and corrected[left] > foot:
            left -= 1
        right = index
        while right < len(corrected) - 1 and corrected[right] > foot:
            right += 1
        width = times[right] - times[left]
        if width <= 0.02:
            continue
        area = integrate_segment(times, values, left, right, baseline)
        if area <= 0:
            continue
        raw_apex_index = max(range(left, right + 1), key=lambda idx: raw_values[idx])
        candidates.append(
            {
                "start_time_min": times[left],
                "apex_time_min": times[raw_apex_index],
                "end_time_min": times[right],
                "height": raw_values[raw_apex_index],
                "area": area,
                "width_min": width,
                "baseline": baseline,
                "signal": signal,
            }
        )
    return sorted(candidates, key=lambda item: item["area"], reverse=True)


def detect_main_peak(
    scans: list[LCMSSpectrumScan],
    rt_start: float | None = None,
    rt_end: float | None = None,
    signal: str = "bpc",
) -> dict[str, float] | None:
    peaks = detect_chromatogram_peaks(scans, rt_start, rt_end, signal)
    if peaks:
        peak = dict(peaks[0])
        peak["match_method"] = "largest_area"
        return peak
    selected = [
        scan for scan in sorted(scans, key=lambda item: item.rt)
        if (rt_start is None or scan.rt >= rt_start) and (rt_end is None or scan.rt <= rt_end)
    ]
    if len(selected) < 3:
        return None
    times = [scan.rt for scan in selected]
    raw_values = [
        scan.base_peak_intensity if signal == "bpc" else scan.tic
        for scan in selected
    ]
    apex = max(range(len(raw_values)), key=lambda idx: raw_values[idx])
    peak = {
        "start_time_min": times[max(0, apex - 1)],
        "apex_time_min": times[apex],
        "end_time_min": times[min(len(times) - 1, apex + 1)],
        "height": raw_values[apex],
        "area": raw_values[apex],
        "width_min": 0.0,
        "baseline": quantile(raw_values, 0.10),
        "signal": signal,
        "match_method": "fallback_max_point",
    }
    return peak


def match_peak_to_reference(
    peaks: list[dict[str, float]],
    reference_apex: float,
    match_window_min: float,
) -> dict[str, float] | None:
    if not peaks:
        return None
    in_window = [
        peak for peak in peaks
        if abs(float(peak["apex_time_min"]) - reference_apex) <= match_window_min
    ]
    if not in_window:
        return None
    largest_area = max(float(peak.get("area", 0.0)) for peak in in_window) or 1.0
    largest_height = max(float(peak.get("height", 0.0)) for peak in in_window) or 1.0

    def score(peak: dict[str, float]) -> float:
        distance_ratio = abs(float(peak["apex_time_min"]) - reference_apex) / max(match_window_min, 1e-9)
        area_ratio = float(peak.get("area", 0.0)) / largest_area
        height_ratio = float(peak.get("height", 0.0)) / largest_height
        return area_ratio * 0.55 + height_ratio * 0.25 + (1.0 - distance_ratio) * 0.20

    matched = dict(max(in_window, key=score))
    matched["match_method"] = "reference_window"
    matched["reference_apex_time_min"] = reference_apex
    matched["reference_delta_min"] = float(matched["apex_time_min"]) - reference_apex
    return matched


def match_peak_landmarks(
    reference_peaks: list[dict[str, float]],
    sample_peaks: list[dict[str, float]],
    match_window_min: float,
    max_landmarks: int = 40,
) -> list[dict[str, float]]:
    if not reference_peaks or not sample_peaks:
        return []
    references = sorted(reference_peaks, key=lambda item: item.get("area", 0.0), reverse=True)[:max_landmarks]
    samples = sorted(sample_peaks, key=lambda item: item.get("area", 0.0), reverse=True)[:max_landmarks * 2]
    used_sample_indexes: set[int] = set()
    landmarks: list[dict[str, float]] = []
    for reference_peak in sorted(references, key=lambda item: item["apex_time_min"]):
        reference_rt = float(reference_peak["apex_time_min"])
        candidates = [
            (index, peak)
            for index, peak in enumerate(samples)
            if index not in used_sample_indexes
            and abs(float(peak["apex_time_min"]) - reference_rt) <= match_window_min
        ]
        if not candidates:
            continue
        index, matched = min(candidates, key=lambda item: abs(float(item[1]["apex_time_min"]) - reference_rt))
        used_sample_indexes.add(index)
        sample_rt = float(matched["apex_time_min"])
        landmarks.append(
            {
                "reference_rt": reference_rt,
                "sample_rt": sample_rt,
                "rt_shift": reference_rt - sample_rt,
                "reference_area": float(reference_peak.get("area", 0.0)),
                "sample_area": float(matched.get("area", 0.0)),
            }
        )
    return landmarks


def robust_shift_from_landmarks(landmarks: list[dict[str, float]]) -> float | None:
    if not landmarks:
        return None
    shifts = sorted(float(item["rt_shift"]) for item in landmarks)
    median_shift = statistics.median(shifts)
    if len(shifts) < 4:
        return median_shift
    deviations = [abs(value - median_shift) for value in shifts]
    mad = statistics.median(deviations)
    if mad <= 1e-9:
        return median_shift
    filtered = [
        value
        for value in shifts
        if abs(value - median_shift) <= max(0.05, 3.5 * mad)
    ]
    return statistics.median(filtered or shifts)


def main_peak_alignment(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    reference_sample: str,
    rt_start: float | None = None,
    rt_end: float | None = None,
    signal: str = "bpc",
    match_window_min: float = 3.0,
) -> dict[str, object]:
    ref_peaks = detect_chromatogram_peaks(scans_by_sample[reference_sample], rt_start, rt_end, signal)
    ref_peak = dict(ref_peaks[0]) if ref_peaks else detect_main_peak(scans_by_sample[reference_sample], rt_start, rt_end, signal)
    if ref_peak is not None and "match_method" not in ref_peak:
        ref_peak["match_method"] = "largest_area"
    ref_apex = ref_peak["apex_time_min"] if ref_peak else 0.0
    apex_by_sample: dict[str, float] = {}
    shift_by_sample: dict[str, float] = {}
    main_peak_by_sample: dict[str, dict[str, float] | None] = {}
    landmarks_by_sample: dict[str, list[dict[str, float]]] = {}
    for sample_id, scans in scans_by_sample.items():
        if sample_id == reference_sample:
            peak = ref_peak
            landmarks = [
                {
                    "reference_rt": float(item["apex_time_min"]),
                    "sample_rt": float(item["apex_time_min"]),
                    "rt_shift": 0.0,
                    "reference_area": float(item.get("area", 0.0)),
                    "sample_area": float(item.get("area", 0.0)),
                }
                for item in ref_peaks
            ]
            shift = 0.0
        else:
            peaks = detect_chromatogram_peaks(scans, rt_start, rt_end, signal)
            landmarks = match_peak_landmarks(ref_peaks, peaks, match_window_min)
            peak = match_peak_to_reference(peaks, float(ref_apex), match_window_min)
            if peak is None:
                peak = detect_main_peak(scans, rt_start, rt_end, signal)
                if peak is not None:
                    peak = dict(peak)
                    peak["match_method"] = "fallback_largest_area"
            multi_peak_shift = robust_shift_from_landmarks(landmarks)
            if multi_peak_shift is not None and len(landmarks) >= 2:
                shift = multi_peak_shift
                if peak is not None:
                    peak = dict(peak)
                    peak["match_method"] = "multi_peak_landmark_median"
                    peak["landmark_count"] = len(landmarks)
            else:
                apex = peak["apex_time_min"] if peak else 0.0
                shift = ref_apex - apex
        main_peak_by_sample[sample_id] = peak
        apex = peak["apex_time_min"] if peak else 0.0
        apex_by_sample[sample_id] = apex
        shift_by_sample[sample_id] = shift
        landmarks_by_sample[sample_id] = landmarks
    return {
        "reference_sample": reference_sample,
        "reference_apex_rt": ref_apex,
        "apex_by_sample": apex_by_sample,
        "main_peak_by_sample": main_peak_by_sample,
        "alignment_landmarks_by_sample": landmarks_by_sample,
        "alignment_landmark_count_by_sample": {
            sample_id: len(landmarks)
            for sample_id, landmarks in landmarks_by_sample.items()
        },
        "rt_shift_by_sample": shift_by_sample,
        "alignment_signal": signal,
        "alignment_rt_start": rt_start,
        "alignment_rt_end": rt_end,
        "alignment_match_window_min": match_window_min,
        "alignment_method": "multi_peak_landmark_median_shift",
    }


def top_spectrum_points(
    scan: LCMSSpectrumScan,
    mz_min: float,
    mz_max: float,
    max_peaks: int,
    min_intensity: float,
) -> tuple[list[float], list[float]]:
    pairs = [
        (mz, intensity)
        for mz, intensity in zip(scan.mz_array, scan.intensity_array)
        if mz_min <= mz <= mz_max and intensity >= min_intensity
    ]
    if len(pairs) > max_peaks:
        pairs = sorted(pairs, key=lambda item: item[1], reverse=True)[:max_peaks]
    pairs.sort(key=lambda item: item[0])
    return [mz for mz, _ in pairs], [intensity for _, intensity in pairs]


def spectrum_payload(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    shifts: dict[str, float],
    mz_min: float,
    mz_max: float,
    max_peaks_per_scan: int,
    min_intensity: float,
) -> dict[str, list[dict[str, object]]]:
    payload: dict[str, list[dict[str, object]]] = {}
    for sample_id, scans in scans_by_sample.items():
        shift = shifts.get(sample_id, 0.0)
        sample_scans: list[dict[str, object]] = []
        for scan in sorted(scans, key=lambda item: item.rt):
            mz_array, intensity_array = top_spectrum_points(
                scan,
                mz_min=mz_min,
                mz_max=mz_max,
                max_peaks=max_peaks_per_scan,
                min_intensity=min_intensity,
            )
            sample_scans.append(
                {
                    "scan_id": scan.scan_id,
                    "raw_file_id": scan.raw_file_id,
                    "rt": scan.rt,
                    "aligned_rt": scan.rt + shift,
                    "tic": scan.tic,
                    "base_peak_mz": scan.base_peak_mz,
                    "base_peak_intensity": scan.base_peak_intensity,
                    "mz": mz_array,
                    "intensity": intensity_array,
                }
            )
        payload[sample_id] = sample_scans
    return payload


def nearest_scan(scans: list[LCMSSpectrumScan], aligned_rt: float, shift: float) -> LCMSSpectrumScan:
    return min(scans, key=lambda scan: abs((scan.rt + shift) - aligned_rt))


def build_intensity_grid(
    scans: list[LCMSSpectrumScan],
    rt_bins: list[float],
    mz_bins: list[float],
    rt_shift: float,
) -> list[list[float]]:
    rt_count = max(0, len(rt_bins) - 1)
    mz_count = max(0, len(mz_bins) - 1)
    grid = [[0.0 for _ in range(rt_count)] for _ in range(mz_count)]
    if rt_count == 0 or mz_count == 0:
        return grid
    rt_min = rt_bins[0]
    rt_step = rt_bins[1] - rt_bins[0]
    mz_min = mz_bins[0]
    mz_step = mz_bins[1] - mz_bins[0]
    for scan in scans:
        aligned_rt = scan.rt + rt_shift
        rt_index = int((aligned_rt - rt_min) / rt_step)
        if rt_index < 0 or rt_index >= rt_count:
            continue
        for mz, intensity in zip(scan.mz_array, scan.intensity_array):
            mz_index = int((mz - mz_min) / mz_step)
            if 0 <= mz_index < mz_count:
                grid[mz_index][rt_index] += intensity
    return grid


def similarity_score(reference_intensity: float, sample_intensity: float) -> float:
    ref = max(0.0, reference_intensity)
    sample = max(0.0, sample_intensity)
    score = 1.0 - abs(ref - sample) / (ref + sample + 1e-9)
    return max(0.0, min(1.0, score))


def neighborhood_radius(rt_step: float, mz_step: float) -> tuple[int, int]:
    """Choose a small local window so minor RT/mz bin drift is not overcalled."""
    rt_radius = max(1, min(5, int(math.ceil(0.15 / max(rt_step, 1e-9)))))
    mz_radius = max(1, min(3, int(math.ceil(1.0 / max(mz_step, 1e-9)))))
    return rt_radius, mz_radius


def local_max(
    grid: list[list[float]],
    mz_index: int,
    rt_index: int,
    mz_radius: int,
    rt_radius: int,
) -> float:
    if not grid:
        return 0.0
    mz_start = max(0, mz_index - mz_radius)
    mz_end = min(len(grid), mz_index + mz_radius + 1)
    best = 0.0
    for mz_i in range(mz_start, mz_end):
        row = grid[mz_i]
        rt_start = max(0, rt_index - rt_radius)
        rt_end = min(len(row), rt_index + rt_radius + 1)
        if rt_start < rt_end:
            best = max(best, max(row[rt_start:rt_end], default=0.0))
    return best


def grid_values(grid: list[list[float]]) -> list[float]:
    return [value for row in grid for value in row if value > 0]


def normalization_factor(grid: list[list[float]], method: str) -> float:
    values = grid_values(grid)
    if method == "none":
        return 1.0
    if method in {"tic", "tic_log1p"}:
        return sum(values) or 1.0
    if method == "max":
        return max(values, default=1.0) or 1.0
    if method == "median":
        return statistics.median(values) if values else 1.0
    raise ValueError(f"Unsupported heatmap normalization method: {method}")


def heatmap_signal_floor(grid: list[list[float]]) -> float:
    values = grid_values(grid)
    if not values:
        return 0.0
    return max(max(values) * 0.01, quantile(values, 0.05))


def normalize_grid(grid: list[list[float]], method: str) -> tuple[list[list[float]], float]:
    factor = normalization_factor(grid, method)
    if method == "none":
        return [[float(value) for value in row] for row in grid], factor
    if method == "tic_log1p":
        scale = 1_000_000.0
        return [[math.log1p(float(value) / factor * scale) for value in row] for row in grid], factor
    return [[float(value) / factor for value in row] for row in grid], factor


def difference_type_for_bin(
    reference_raw: float,
    sample_raw: float,
    reference_norm: float,
    sample_norm: float,
    score: float,
) -> str:
    if reference_raw <= 0 and sample_raw > 0:
        return "new_signal"
    if reference_raw > 0 and sample_raw <= 0:
        return "missing_signal"
    if score > 0.75:
        return "uncertain"
    fold = (sample_norm + 1e-12) / (reference_norm + 1e-12)
    if fold >= 1.5:
        return "intensity_increased"
    if fold <= 1 / 1.5:
        return "intensity_decreased"
    return "uncertain"


def difference_region(
    reference_sample: str,
    sample_id: str,
    rt_bins: list[float],
    mz_bins: list[float],
    rt_index: int,
    mz_index: int,
    score: float,
    reference_raw: float,
    sample_raw: float,
    reference_norm: float,
    sample_norm: float,
) -> dict[str, object]:
    raw_values = [reference_raw, sample_raw]
    normalized_values = [reference_norm, sample_norm]
    max_raw = max(raw_values)
    mean_norm = statistics.mean(normalized_values)
    fold_change = (sample_norm + 1e-12) / (reference_norm + 1e-12)
    present = [
        sample
        for sample, value in ((reference_sample, reference_raw), (sample_id, sample_raw))
        if value > 0
    ]
    return {
        "rt": (rt_bins[rt_index] + rt_bins[rt_index + 1]) / 2.0,
        "rt_start": rt_bins[rt_index],
        "rt_end": rt_bins[rt_index + 1],
        "mz": (mz_bins[mz_index] + mz_bins[mz_index + 1]) / 2.0,
        "mz_start": mz_bins[mz_index],
        "mz_end": mz_bins[mz_index + 1],
        "score": score,
        "difference_score": 1.0 - score,
        "reference_intensity": reference_norm,
        "sample_intensity": sample_norm,
        "reference_raw_intensity": reference_raw,
        "sample_raw_intensity": sample_raw,
        "max_intensity": max_raw,
        "sample_presence": present,
        "group_mean_intensity": mean_norm,
        "fold_change": fold_change,
        "difference_type": difference_type_for_bin(
            reference_raw,
            sample_raw,
            reference_norm,
            sample_norm,
            score,
        ),
        "saved_as_feature": False,
    }


def cohort_difference_type(
    raw_values: dict[str, float],
    normalized_values: dict[str, float],
    score: float,
) -> str:
    present = [value for value in raw_values.values() if value > 0]
    if present and len(present) < len(raw_values):
        return "new_signal"
    if score > 0.75:
        return "uncertain"
    positive = [value for value in normalized_values.values() if value > 0]
    if len(positive) >= 2 and (max(positive) + 1e-12) / (min(positive) + 1e-12) >= 1.5:
        return "intensity_increased"
    return "uncertain"


def signal_status(raw_values: dict[str, float]) -> str:
    present_count = sum(1 for value in raw_values.values() if value > 0)
    if present_count == 0:
        return "low_signal_background"
    if present_count < len(raw_values):
        return "presence_absence"
    return "high_confidence"


def cohort_difference_region(
    sample_ids: list[str],
    rt_bins: list[float],
    mz_bins: list[float],
    rt_index: int,
    mz_index: int,
    score: float,
    raw_values: dict[str, float],
    normalized_values: dict[str, float],
) -> dict[str, object]:
    normalized_list = [normalized_values[sample_id] for sample_id in sample_ids]
    raw_list = [raw_values[sample_id] for sample_id in sample_ids]
    positive_norm = [value for value in normalized_list if value > 0]
    mean_norm = statistics.mean(normalized_list) if normalized_list else 0.0
    cv = statistics.pstdev(normalized_list) / mean_norm if mean_norm > 0 else None
    fold_change = (
        (max(positive_norm) + 1e-12) / (min(positive_norm) + 1e-12)
        if len(positive_norm) >= 2
        else None
    )
    return {
        "rt": (rt_bins[rt_index] + rt_bins[rt_index + 1]) / 2.0,
        "rt_start": rt_bins[rt_index],
        "rt_end": rt_bins[rt_index + 1],
        "mz": (mz_bins[mz_index] + mz_bins[mz_index + 1]) / 2.0,
        "mz_start": mz_bins[mz_index],
        "mz_end": mz_bins[mz_index + 1],
        "score": score,
        "difference_score": 1.0 - score,
        "sample_intensities": {
            sample_id: {
                "raw": raw_values[sample_id],
                "normalized": normalized_values[sample_id],
            }
            for sample_id in sample_ids
        },
        "max_intensity": max(raw_list, default=0.0),
        "sample_presence": [
            sample_id for sample_id in sample_ids if raw_values[sample_id] > 0
        ],
        "group_mean_intensity": mean_norm,
        "cv": cv,
        "signal_status": signal_status(raw_values),
        "fold_change": fold_change,
        "difference_type": cohort_difference_type(raw_values, normalized_values, score),
        "saved_as_feature": False,
    }


def apply_abundance_weight_to_regions(regions: list[dict[str, object]]) -> list[dict[str, object]]:
    positive_means = [
        float(region.get("group_mean_intensity") or 0.0)
        for region in regions
        if float(region.get("group_mean_intensity") or 0.0) > 0
    ]
    abundance_reference = quantile(positive_means, 0.90) if positive_means else 0.0
    for region in regions:
        mean_intensity = float(region.get("group_mean_intensity") or 0.0)
        difference_score = float(region.get("difference_score") or 0.0)
        if abundance_reference <= 0 or mean_intensity <= 0:
            abundance_weight = 0.0
        else:
            abundance_weight = min(1.0, math.sqrt(mean_intensity / abundance_reference))
        weighted_difference = difference_score * abundance_weight
        region["abundance_weight"] = abundance_weight
        region["weighted_difference_score"] = weighted_difference
        region["ranking_score"] = weighted_difference
    return regions


def sort_difference_regions(regions: list[dict[str, object]]) -> list[dict[str, object]]:
    apply_abundance_weight_to_regions(regions)
    regions.sort(
        key=lambda item: (
            -float(item.get("weighted_difference_score") or 0.0),
            -float(item.get("difference_score") or 0.0),
            -float(item.get("max_intensity") or 0.0),
        )
    )
    return regions


def top_spatially_distinct_regions(
    regions: list[dict[str, object]],
    rt_radius_bins: int,
    mz_radius_bins: int,
    limit: int = 80,
) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    min_rt_gap = max(1, rt_radius_bins * 2 + 1)
    min_mz_gap = max(1, mz_radius_bins * 2 + 1)
    for region in regions:
        rt_index = int(region.get("rt_index", -10_000))
        mz_index = int(region.get("mz_index", -10_000))
        if any(
            abs(rt_index - int(item.get("rt_index", 10_000))) < min_rt_gap
            and abs(mz_index - int(item.get("mz_index", 10_000))) < min_mz_gap
            for item in selected
        ):
            continue
        selected.append(dict(region))
        if len(selected) >= limit:
            break
    cleaned_regions: list[dict[str, object]] = []
    for region in selected:
        cleaned = dict(region)
        cleaned.pop("rt_index", None)
        cleaned.pop("mz_index", None)
        cleaned_regions.append(cleaned)
    return cleaned_regions


def build_similarity_heatmaps(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    reference_sample: str,
    shifts: dict[str, float],
    rt_min: float,
    rt_max: float,
    mz_min: float,
    mz_max: float,
    rt_bin_count: int = 180,
    mz_bin_count: int = 160,
    normalization_methods: list[str] | None = None,
) -> list[dict[str, object]]:
    methods = normalization_methods or ["max"]
    rt_step = (rt_max - rt_min) / rt_bin_count
    mz_step = (mz_max - mz_min) / mz_bin_count
    rt_neighborhood_bins, mz_neighborhood_bins = neighborhood_radius(rt_step, mz_step)
    rt_bins = [rt_min + i * rt_step for i in range(rt_bin_count + 1)]
    mz_bins = [mz_min + i * mz_step for i in range(mz_bin_count + 1)]
    sample_ids = list(scans_by_sample)
    raw_grids = {
        sample_id: build_intensity_grid(
            scans,
            rt_bins,
            mz_bins,
            shifts.get(sample_id, 0.0),
        )
        for sample_id, scans in scans_by_sample.items()
    }
    signal_floors = {
        sample_id: heatmap_signal_floor(grid)
        for sample_id, grid in raw_grids.items()
    }
    heatmaps: list[dict[str, object]] = []
    if len(sample_ids) >= 2:
        for method_index, method in enumerate(methods):
            normalized_by_sample = {
                sample_id: normalize_grid(raw_grids[sample_id], method)
                for sample_id in sample_ids
            }
            norm_grids = {
                sample_id: normalized_by_sample[sample_id][0]
                for sample_id in sample_ids
            }
            norm_factors = {
                sample_id: normalized_by_sample[sample_id][1]
                for sample_id in sample_ids
            }
            scores = []
            differences = []
            for mz_index in range(mz_bin_count):
                row = []
                for rt_index in range(rt_bin_count):
                    raw_values = {
                        sample_id: (
                            0.0
                            if (
                                local := local_max(
                                    raw_grids[sample_id],
                                    mz_index,
                                    rt_index,
                                    mz_neighborhood_bins,
                                    rt_neighborhood_bins,
                                )
                            ) < signal_floors.get(sample_id, 0.0)
                            else local
                        )
                        for sample_id in sample_ids
                    }
                    if max(raw_values.values(), default=0.0) <= 0:
                        row.append(None)
                        continue
                    values = [
                        (
                            0.0
                            if raw_values[sample_id] <= 0
                            else local_max(
                                norm_grids[sample_id],
                                mz_index,
                                rt_index,
                                mz_neighborhood_bins,
                                rt_neighborhood_bins,
                            )
                        )
                        for sample_id in sample_ids
                    ]
                    mean_value = statistics.mean(values) if values else 0.0
                    cv = statistics.pstdev(values) / mean_value if mean_value > 0 else float("inf")
                    score = 0.0 if not math.isfinite(cv) else 1.0 / (1.0 + cv)
                    row.append(score)
                    if score < 0.75:
                        norm_values = {
                            sample_id: values[index]
                            for index, sample_id in enumerate(sample_ids)
                        }
                        region = cohort_difference_region(
                            sample_ids,
                            rt_bins,
                            mz_bins,
                            rt_index,
                            mz_index,
                            score,
                            raw_values,
                            norm_values,
                        )
                        region["rt_index"] = rt_index
                        region["mz_index"] = mz_index
                        differences.append(region)
                scores.append(row)
            sort_difference_regions(differences)
            heatmap = {
                "comparison_type": "cohort_cv",
                "comparison_key": "cohort:all_samples",
                "comparison_label": "All samples CV similarity",
                "reference_sample": reference_sample,
                "sample_id": "__cohort__",
                "sample_ids": sample_ids,
                "normalization_method": method,
                "normalization_factors_by_sample": norm_factors,
                "signal_floor_by_sample": signal_floors,
                "rt_min": rt_min,
                "rt_max": rt_max,
                "mz_min": mz_min,
                "mz_max": mz_max,
                "rt_bin_count": rt_bin_count,
                "mz_bin_count": mz_bin_count,
                "rt_neighborhood_bins": rt_neighborhood_bins,
                "mz_neighborhood_bins": mz_neighborhood_bins,
                "similarity_method": "local_neighborhood_max",
                "scores": scores,
                "sample_intensity_grids": norm_grids,
                "low_similarity_points": top_spatially_distinct_regions(
                    differences,
                    rt_neighborhood_bins,
                    mz_neighborhood_bins,
                    80,
                ),
            }
            if method_index == 0:
                heatmap["sample_raw_intensity_grids"] = raw_grids
            heatmaps.append(heatmap)
    return heatmaps


def effective_aligned_rt_range(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    shifts: dict[str, float],
    requested_min: float,
    requested_max: float,
) -> tuple[float, float]:
    aligned_rts = [
        scan.rt + float(shifts.get(sample_id, 0.0))
        for sample_id, scans in scans_by_sample.items()
        for scan in scans
    ]
    if not aligned_rts:
        return requested_min, requested_max
    actual_min = min(aligned_rts)
    actual_max = max(aligned_rts)
    rt_min = max(requested_min, actual_min)
    rt_max = min(requested_max, actual_max)
    if rt_max <= rt_min:
        return actual_min, actual_max
    return rt_min, rt_max


def prepare_workbench_payload(
    raw_files: list[LCMSRawFile],
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    reference_sample: str,
    rt_min: float,
    rt_max: float,
    mz_min: float,
    mz_max: float,
    alignment_rt_start: float | None,
    alignment_rt_end: float | None,
    alignment_signal: str,
    alignment_match_window_min: float,
    max_peaks_per_scan: int,
    spectrum_min_intensity: float,
    rt_bin_count: int,
    mz_bin_count: int,
    heatmap_normalization_methods: list[str] | None = None,
) -> dict[str, object]:
    alignment = main_peak_alignment(
        scans_by_sample,
        reference_sample,
        alignment_rt_start,
        alignment_rt_end,
        alignment_signal,
        alignment_match_window_min,
    )
    shifts = alignment["rt_shift_by_sample"]
    assert isinstance(shifts, dict)
    sample_ids = [raw_file.sample_id for raw_file in raw_files]
    effective_rt_min, effective_rt_max = effective_aligned_rt_range(
        scans_by_sample,
        {sample_id: float(value) for sample_id, value in shifts.items()},
        rt_min,
        rt_max,
    )
    chromatograms = {
        sample_id: {
            "raw": chromatogram_points(scans_by_sample[sample_id], 0.0),
            "aligned": chromatogram_points(scans_by_sample[sample_id], float(shifts.get(sample_id, 0.0))),
        }
        for sample_id in sample_ids
    }
    spectra = spectrum_payload(
        scans_by_sample,
        {sample_id: float(value) for sample_id, value in shifts.items()},
        mz_min,
        mz_max,
        max_peaks_per_scan,
        spectrum_min_intensity,
    )
    heatmaps = build_similarity_heatmaps(
        scans_by_sample,
        reference_sample,
        {sample_id: float(value) for sample_id, value in shifts.items()},
        effective_rt_min,
        effective_rt_max,
        mz_min,
        mz_max,
        rt_bin_count,
        mz_bin_count,
        heatmap_normalization_methods,
    )
    return {
        "raw_files": [asdict(raw_file) for raw_file in raw_files],
        "sample_ids": sample_ids,
        "reference_sample": reference_sample,
        "rt_min": effective_rt_min,
        "rt_max": effective_rt_max,
        "requested_rt_min": rt_min,
        "requested_rt_max": rt_max,
        "mz_min": mz_min,
        "mz_max": mz_max,
        "alignment": alignment,
        "chromatograms": chromatograms,
        "spectra": spectra,
        "heatmaps": heatmaps,
        "summary": {
            "scan_count_by_sample": {sample_id: len(scans_by_sample[sample_id]) for sample_id in sample_ids},
            "main_peak_rt_cv": _main_peak_cv(alignment),
        },
    }


def _main_peak_cv(alignment: dict[str, object]) -> float | None:
    apexes = list((alignment.get("apex_by_sample") or {}).values())
    if len(apexes) < 2:
        return None
    values = [float(value) for value in apexes]
    mean = statistics.mean(values)
    return statistics.pstdev(values) / mean if mean else None
