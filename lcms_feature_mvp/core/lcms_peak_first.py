#!/usr/bin/env python3
"""Peak-first LC-MS comparison workflow.

This module intentionally lives beside the existing RT-m/z heatmap workflow.
It implements a coarse-to-fine path: detect TIC peaks first, score each peak by
chromatogram and within-peak spectrum consistency, then expose changed m/z
signals for XIC confirmation.
"""

from __future__ import annotations

import math
import statistics
from bisect import bisect_left
from dataclasses import dataclass
from typing import Iterable

from .lcms_models import LCMSRawFile, LCMSSpectrumScan
from .lcms_parser import raw_file_payload
from .lcms_workbench import chromatogram_points, main_peak_alignment


@dataclass
class PeakFirstParams:
    smoothing_window: int = 3
    baseline_window: int = 41
    min_snr: float = 5.0
    min_prominence_factor: float = 3.0
    min_area: float = 0.0
    min_area_ratio: float = 0.00001
    min_width: float = 0.03
    max_width: float = 3.0
    min_apex_distance_min: float = 0.02
    peak_boundary_fraction: float = 0.08
    peak_body_smoothing_window: int = 11
    min_valley_depth_fraction: float = 0.05
    shoulder_merge_gap_min: float = 0.30
    shoulder_merge_valley_depth_fraction: float = 0.35
    min_sample_presence_count: int = 2
    max_local_shift_min: float = 0.1
    max_local_shift_fraction: float = 0.2
    max_local_shift_cap_min: float = 0.5
    continuous_peak_gap_min: float = 0.08
    isolated_peak_edge_extension_min: float = 0.03
    tic_difference_selection_weight: float = 1.5
    xic_feature_search_margin_min: float = 0.45
    feature_drift_tolerance_min: float = 0.35
    mz_tolerance_ppm: float = 20.0
    mz_tolerance_da: float = 0.5
    mz_tolerance_mode: str = "da"
    top_n_mz: int = 50
    top_n_changed_mz: int = 20
    min_changed_mz_gap_da: float = 1.0
    top_n_peaks: int = 60
    max_spectrum_points_per_scan: int = 200
    min_relative_intensity: float = 0.001
    resample_points: int = 100
    peak_margin_min: float = 0.08
    peak_source: str = "consensus"
    consensus_rt_step_min: float = 0.01
    consensus_max_points: int = 8000
    xic_context_min: float = 3.0
    feature_rt_tolerance_min: float = 0.12
    feature_mz_weight: float = 0.7
    feature_rt_weight: float = 0.3
    feature_min_height_fraction: float = 0.02
    feature_gap_fill_height_fraction: float = 0.005
    feature_abundance_weight: float = 0.8
    global_feature_rt_merge_tolerance_min: float = 0.8
    global_feature_mz_merge_tolerance_factor: float = 2.0


def moving_average(values: list[float], window: int) -> list[float]:
    if window <= 1 or not values:
        return [float(value) for value in values]
    radius = max(1, window // 2)
    smoothed: list[float] = []
    for index in range(len(values)):
        start = max(0, index - radius)
        end = min(len(values), index + radius + 1)
        smoothed.append(sum(values[start:end]) / max(1, end - start))
    return smoothed


def rolling_min(values: list[float], window: int) -> list[float]:
    if window <= 1 or not values:
        return [min(values) if values else 0.0 for _ in values]
    radius = max(1, window // 2)
    result: list[float] = []
    for index in range(len(values)):
        start = max(0, index - radius)
        end = min(len(values), index + radius + 1)
        result.append(min(values[start:end]))
    return result


def quantile(values: list[float], q: float) -> float:
    clean = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not clean:
        return 0.0
    pos = max(0.0, min(1.0, q)) * (len(clean) - 1)
    low = math.floor(pos)
    high = math.ceil(pos)
    if low == high:
        return clean[low]
    return clean[low] * (1 - (pos - low)) + clean[high] * (pos - low)


def score_from_cv(values: list[float]) -> float:
    positive = [float(value) for value in values if float(value) > 0]
    if not positive:
        return 0.0
    mean = statistics.mean(positive)
    if mean <= 0:
        return 0.0
    cv = statistics.pstdev(positive) / mean if len(positive) > 1 else 0.0
    return 1.0 / (1.0 + cv)


def cosine(a: list[float], b: list[float]) -> float:
    if not a or not b:
        return 0.0
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    if na <= 0 or nb <= 0:
        return 0.0
    return max(0.0, min(1.0, dot / (na * nb)))


def align_rt_by_main_or_feature_peak(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    reference_sample: str,
    reference_peak: float | None = None,
    signal: str = "tic",
) -> dict[str, object]:
    alignment = main_peak_alignment(
        scans_by_sample,
        reference_sample=reference_sample,
        signal=signal,
        match_window_min=3.0,
    )
    if reference_peak is None:
        return alignment
    shifts: dict[str, float] = {}
    apex_by_sample: dict[str, float] = {}
    for sample_id, scans in scans_by_sample.items():
        if not scans:
            shifts[sample_id] = 0.0
            apex_by_sample[sample_id] = 0.0
            continue
        nearest = min(scans, key=lambda scan: abs(scan.rt - reference_peak))
        apex_by_sample[sample_id] = nearest.rt
        shifts[sample_id] = reference_peak - nearest.rt
    alignment["rt_shift_by_sample"] = shifts
    alignment["apex_by_sample"] = apex_by_sample
    alignment["reference_apex_rt"] = reference_peak
    alignment["alignment_method"] = "user_feature_peak_shift"
    return alignment


def estimate_baseline_and_noise(
    tic_curve: list[tuple[float, float]],
    params: PeakFirstParams,
) -> dict[str, object]:
    intensities = [float(value) for _, value in tic_curve]
    smoothed = moving_average(intensities, params.smoothing_window)
    baseline = rolling_min(smoothed, params.baseline_window)
    corrected = [max(value - base, 0.0) for value, base in zip(smoothed, baseline)]
    body_corrected = moving_average(corrected, params.peak_body_smoothing_window)
    low = [value for value in corrected if value <= quantile(corrected, 0.4)]
    median = statistics.median(low) if low else 0.0
    mad = statistics.median([abs(value - median) for value in low]) if low else 0.0
    noise = max(median + 1.4826 * mad, quantile(corrected, 0.1), 1e-9)
    return {
        "rt": [rt for rt, _ in tic_curve],
        "raw": intensities,
        "smoothed": smoothed,
        "baseline": baseline,
        "corrected": corrected,
        "body_corrected": body_corrected,
        "noise_level": noise,
    }


def detect_tic_peaks(
    tic_curve: list[tuple[float, float]],
    params: PeakFirstParams,
) -> list[dict[str, object]]:
    estimated = estimate_baseline_and_noise(tic_curve, params)
    rt = list(estimated["rt"])
    corrected = list(estimated["corrected"])
    body_corrected = list(estimated["body_corrected"])
    noise = float(estimated["noise_level"])
    total_area = 0.0
    for i in range(1, len(rt)):
        total_area += (corrected[i - 1] + corrected[i]) * 0.5 * max(rt[i] - rt[i - 1], 0.0)
    peaks: list[dict[str, object]] = []
    min_height = max(noise * params.min_snr, 1e-9)
    if len(rt) < 3:
        return []
    rt_steps = [rt[i] - rt[i - 1] for i in range(1, len(rt)) if rt[i] > rt[i - 1]]
    median_step = statistics.median(rt_steps) if rt_steps else 0.01
    min_apex_distance_points = max(1, int(params.min_apex_distance_min / max(median_step, 1e-9)))
    candidates: list[int] = []
    for index in range(1, len(body_corrected) - 1):
        value = body_corrected[index]
        if value < min_height or value < body_corrected[index - 1] or value < body_corrected[index + 1]:
            continue
        candidates.append(index)
    selected: list[int] = []
    for index in sorted(candidates, key=lambda item: body_corrected[item], reverse=True):
        if any(abs(index - kept) <= min_apex_distance_points for kept in selected):
            continue
        selected.append(index)
    selected.sort()

    # Collapse shoulders inside the same chromatographic peak. If the valley
    # between adjacent apices is shallow, keep the higher apex instead of
    # creating duplicate peak windows with nearly identical RT ranges.
    collapsed: list[int] = []
    for index in selected:
        if not collapsed:
            collapsed.append(index)
            continue
        previous = collapsed[-1]
        start, end = sorted((previous, index))
        valley = min(body_corrected[start : end + 1])
        lower_apex = min(body_corrected[previous], body_corrected[index])
        valley_depth = (lower_apex - valley) / lower_apex if lower_apex > 0 else 0.0
        close_shoulder = (rt[end] - rt[start]) <= params.shoulder_merge_gap_min
        if valley_depth < params.min_valley_depth_fraction or (
            close_shoulder and valley_depth < params.shoulder_merge_valley_depth_fraction
        ):
            if body_corrected[index] > body_corrected[previous]:
                collapsed[-1] = index
        else:
            collapsed.append(index)
    selected = collapsed

    for position, index in enumerate(selected):
        value = body_corrected[index]
        boundary_level = max(noise, value * params.peak_boundary_fraction)
        left = index
        while left > 0 and body_corrected[left] > boundary_level:
            left -= 1
        right = index
        while right < len(body_corrected) - 1 and body_corrected[right] > boundary_level:
            right += 1
        if position > 0:
            previous_apex = selected[position - 1]
            valley_left = min(range(previous_apex, index + 1), key=lambda item: body_corrected[item])
            left = max(left, valley_left)
        if position < len(selected) - 1:
            next_apex = selected[position + 1]
            valley_right = min(range(index, next_apex + 1), key=lambda item: body_corrected[item])
            right = min(right, valley_right)
        if right <= left:
            continue
        apex_index = max(range(left, right + 1), key=lambda item: body_corrected[item])
        apex_height = corrected[apex_index]
        body_height = body_corrected[apex_index]
        width = rt[right] - rt[left]
        area = 0.0
        for i in range(left + 1, right + 1):
            area += (corrected[i - 1] + corrected[i]) * 0.5 * max(rt[i] - rt[i - 1], 0.0)
        edge_level = max(body_corrected[left], body_corrected[right])
        prominence = max(body_height - edge_level, 0.0)
        peaks.append(
            {
                "tic_peak_id": f"TICP_{len(peaks) + 1:04d}",
                "rt_start": rt[left],
                "rt_apex": rt[apex_index],
                "rt_end": rt[right],
                "apex_intensity": apex_height,
                "height": apex_height,
                "body_height": body_height,
                "area": area,
                "width": width,
                "prominence": prominence,
                "snr": body_height / noise if noise > 0 else 0.0,
                "area_ratio": area / total_area if total_area > 0 else 0.0,
                "noise_level": noise,
                "peak_quality": "unclassified",
                "status": "unclassified",
                "source": "consensus_tic",
            }
        )
    peaks.sort(key=lambda item: float(item["area"]), reverse=True)
    for rank, peak in enumerate(peaks, start=1):
        peak["area_rank"] = rank
    peaks.sort(key=lambda item: float(item["rt_apex"]))
    for rank, peak in enumerate(peaks, start=1):
        peak["tic_peak_id"] = f"TICP_{rank:04d}"
        peak["rt_rank"] = rank
    return peaks


def _recalculate_peak_widths(peaks: list[dict[str, object]]) -> None:
    for peak in peaks:
        peak["width"] = max(float(peak["rt_end"]) - float(peak["rt_start"]), 0.0)


def refine_tic_peak_boundaries(
    tic_peaks: list[dict[str, object]],
    params: PeakFirstParams,
) -> list[dict[str, object]]:
    """Make TIC peak windows more integration-friendly without merging peaks."""
    peaks = [dict(peak) for peak in sorted(tic_peaks, key=lambda item: float(item.get("rt_apex") or 0.0))]
    if not peaks:
        return peaks
    # Close tiny gaps between continuous peaks so integrations meet at one boundary.
    for index in range(len(peaks) - 1):
        left = peaks[index]
        right = peaks[index + 1]
        gap = float(right["rt_start"]) - float(left["rt_end"])
        if gap <= params.continuous_peak_gap_min:
            boundary = (float(left["rt_end"]) + float(right["rt_start"])) * 0.5
            if gap < 0:
                boundary = (float(left["rt_apex"]) + float(right["rt_apex"])) * 0.5
            left["rt_end"] = min(boundary, float(right["rt_apex"]))
            right["rt_start"] = max(boundary, float(left["rt_apex"]))
    # For isolated peaks, expand a little toward the local blank space. This
    # catches shoulders/tails that thresholding would otherwise trim too hard.
    for index, peak in enumerate(peaks):
        previous_end = float(peaks[index - 1]["rt_end"]) if index > 0 else None
        next_start = float(peaks[index + 1]["rt_start"]) if index + 1 < len(peaks) else None
        start = float(peak["rt_start"])
        end = float(peak["rt_end"])
        extension = params.isolated_peak_edge_extension_min
        if previous_end is None:
            peak["rt_start"] = max(0.0, start - extension)
        else:
            gap = start - previous_end
            if gap > params.continuous_peak_gap_min:
                peak["rt_start"] = start - min(extension, gap * 0.4)
        if next_start is None:
            peak["rt_end"] = end + extension
        else:
            gap = next_start - end
            if gap > params.continuous_peak_gap_min:
                peak["rt_end"] = end + min(extension, gap * 0.4)
    _recalculate_peak_widths(peaks)
    return peaks


def filter_tic_peaks(
    tic_peaks: list[dict[str, object]],
    params: PeakFirstParams,
) -> list[dict[str, object]]:
    for peak in tic_peaks:
        snr = float(peak.get("snr") or 0.0)
        prominence = float(peak.get("prominence") or 0.0)
        noise = float(peak.get("noise_level") or 0.0)
        area = float(peak.get("area") or 0.0)
        area_ratio = float(peak.get("area_ratio") or 0.0)
        width = float(peak.get("width") or 0.0)
        if snr < max(2.0, params.min_snr * 0.6) or width < params.min_width * 0.5:
            peak["status"] = "noise_peak"
            peak["peak_quality"] = "noise"
        elif (
            snr < params.min_snr
            or prominence < noise * params.min_prominence_factor
            or area < params.min_area
            or area_ratio < params.min_area_ratio
            or width < params.min_width
            or width > params.max_width
        ):
            peak["status"] = "weak_peak"
            peak["peak_quality"] = "weak"
        else:
            peak["status"] = "confirmed_peak"
            peak["peak_quality"] = "confirmed"
    return tic_peaks


def tic_curve_for_sample(
    scans: list[LCMSSpectrumScan],
    rt_shift: float = 0.0,
    signal: str = "tic",
) -> list[tuple[float, float]]:
    if signal == "bpc":
        return [(scan.rt + rt_shift, scan.base_peak_intensity) for scan in scans]
    return [(scan.rt + rt_shift, scan.tic) for scan in scans]


def sample_union_tic_peaks(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    shifts: dict[str, float],
    params: PeakFirstParams,
    signal: str = "tic",
) -> list[dict[str, object]]:
    peaks: list[dict[str, object]] = []
    for sample_id, scans in scans_by_sample.items():
        curve = tic_curve_for_sample(scans, shifts.get(sample_id, 0.0), signal)
        sample_peaks = filter_tic_peaks(refine_tic_peak_boundaries(detect_tic_peaks(curve, params), params), params)
        for peak in sample_peaks:
            if peak.get("status") not in {"confirmed_peak", "weak_peak"}:
                continue
            item = dict(peak)
            item["source"] = "sample_tic"
            item["source_sample_id"] = sample_id
            item["supporting_samples"] = [sample_id]
            peaks.append(item)
    return peaks


def merge_tic_peak_candidates(
    consensus_peaks: list[dict[str, object]],
    sample_peaks: list[dict[str, object]],
    params: PeakFirstParams,
) -> list[dict[str, object]]:
    merged: list[dict[str, object]] = []
    candidates = [dict(peak) for peak in consensus_peaks] + [dict(peak) for peak in sample_peaks]
    candidates.sort(key=lambda item: (float(item.get("rt_apex") or 0.0), -float(item.get("area") or 0.0)))
    for peak in candidates:
        apex = float(peak.get("rt_apex") or 0.0)
        width = max(float(peak.get("width") or params.min_width), params.min_width)
        same_peak_tolerance = max(params.min_apex_distance_min * 2.0, min(params.max_local_shift_min, width * 0.35))
        match: dict[str, object] | None = None
        for existing in merged:
            existing_apex = float(existing.get("rt_apex") or 0.0)
            if abs(apex - existing_apex) <= same_peak_tolerance:
                match = existing
                break
        if match is None:
            peak.setdefault("supporting_samples", [])
            if peak.get("source_sample_id") and peak["source_sample_id"] not in peak["supporting_samples"]:
                peak["supporting_samples"].append(peak["source_sample_id"])
            merged.append(peak)
            continue
        if float(peak.get("area") or 0.0) > float(match.get("area") or 0.0):
            for key in ("rt_start", "rt_apex", "rt_end", "apex_intensity", "height", "body_height", "area", "width", "prominence", "snr", "area_ratio"):
                match[key] = peak.get(key, match.get(key))
        match["rt_start"] = min(float(match["rt_start"]), float(peak["rt_start"]))
        match["rt_end"] = max(float(match["rt_end"]), float(peak["rt_end"]))
        sources = set(match.get("supporting_samples") or [])
        if peak.get("source_sample_id"):
            sources.add(str(peak["source_sample_id"]))
        for sample in peak.get("supporting_samples") or []:
            sources.add(str(sample))
        match["supporting_samples"] = sorted(sources)
        if match.get("source") != "consensus_tic" and peak.get("source") == "consensus_tic":
            match["source"] = "consensus_tic"
    merged = refine_tic_peak_boundaries(merged, params)
    for rank, peak in enumerate(sorted(merged, key=lambda item: float(item.get("rt_apex") or 0.0)), start=1):
        peak["tic_peak_id"] = f"TICP_{rank:04d}"
        peak["rt_rank"] = rank
        peak["width"] = max(float(peak["rt_end"]) - float(peak["rt_start"]), 0.0)
    return sorted(merged, key=lambda item: float(item.get("rt_apex") or 0.0))


def consensus_tic_curve(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    shifts: dict[str, float],
    signal: str = "tic",
    max_points: int = 1600,
    target_step_min: float | None = None,
) -> list[tuple[float, float]]:
    points = [
        (rt, intensity)
        for sample_id, scans in scans_by_sample.items()
        for rt, intensity in tic_curve_for_sample(scans, shifts.get(sample_id, 0.0), signal)
    ]
    if not points:
        return []
    points.sort()
    rt_min, rt_max = points[0][0], points[-1][0]
    if rt_max <= rt_min:
        return points
    if target_step_min and target_step_min > 0:
        bins = min(max_points, max(100, int(math.ceil((rt_max - rt_min) / target_step_min))))
    else:
        bins = min(max_points, max(100, int(math.sqrt(len(points))) * 4))
    sums = [0.0 for _ in range(bins)]
    counts = [0 for _ in range(bins)]
    step = (rt_max - rt_min) / bins
    for rt, intensity in points:
        idx = min(bins - 1, max(0, int((rt - rt_min) / step)))
        sums[idx] += intensity
        counts[idx] += 1
    return [
        (rt_min + (idx + 0.5) * step, sums[idx] / counts[idx])
        for idx in range(bins)
        if counts[idx] > 0
    ]


def extract_peak_window(
    scans: list[LCMSSpectrumScan],
    tic_peak: dict[str, object],
    rt_shift: float = 0.0,
    aligned: bool = True,
    margin: float = 0.05,
    signal: str = "tic",
) -> list[tuple[float, float]]:
    start = float(tic_peak["rt_start"]) - margin
    end = float(tic_peak["rt_end"]) + margin
    values = []
    for scan in scans:
        rt = scan.rt + rt_shift if aligned else scan.rt
        if start <= rt <= end:
            values.append((rt, scan.base_peak_intensity if signal == "bpc" else scan.tic))
    return values


def local_align_peak_by_apex(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    tic_peak: dict[str, object],
    shifts: dict[str, float],
    params: PeakFirstParams,
    signal: str = "tic",
) -> dict[str, object]:
    apex_by_sample: dict[str, float] = {}
    windows: dict[str, list[tuple[float, float]]] = {}
    for sample_id, scans in scans_by_sample.items():
        window = extract_peak_window(scans, tic_peak, shifts.get(sample_id, 0.0), True, params.peak_margin_min, signal)
        windows[sample_id] = window
        apex_by_sample[sample_id] = max(window, key=lambda item: item[1])[0] if window else float(tic_peak["rt_apex"])
    reference_apex = statistics.median(apex_by_sample.values()) if apex_by_sample else float(tic_peak["rt_apex"])
    max_shift = min(
        params.max_local_shift_cap_min,
        max(params.max_local_shift_min, params.max_local_shift_fraction * float(tic_peak.get("width") or 0.1)),
    )
    local_shifts: dict[str, float] = {}
    capped: dict[str, bool] = {}
    for sample_id, apex in apex_by_sample.items():
        shift = reference_apex - apex
        if abs(shift) > max_shift:
            local_shifts[sample_id] = 0.0
            capped[sample_id] = True
        else:
            local_shifts[sample_id] = shift
            capped[sample_id] = False
    return {
        "reference_peak_apex_rt": reference_apex,
        "raw_peak_apex_rt_by_sample": apex_by_sample,
        "local_peak_shift_by_sample": local_shifts,
        "local_shift_capped_by_sample": capped,
        "max_local_shift_allowed": max_shift,
        "local_shift_max": max((abs(v) for v in local_shifts.values()), default=0.0),
        "windows": windows,
    }


def sample_specific_cut_bounds(
    tic_peak: dict[str, object],
    local_shifts: dict[str, float],
    params: PeakFirstParams,
) -> dict[str, dict[str, float]]:
    start = float(tic_peak["rt_start"]) - params.peak_margin_min
    end = float(tic_peak["rt_end"]) + params.peak_margin_min
    return {
        sample_id: {
            "aligned_rt_start": start,
            "aligned_rt_end": end,
            # raw/previously global-aligned scan is compared after adding local
            # shift, so the sample-specific cut frame moves by -local_shift.
            "sample_corrected_rt_start": start - float(local_shift),
            "sample_corrected_rt_end": end - float(local_shift),
            "local_rt_shift": float(local_shift),
        }
        for sample_id, local_shift in local_shifts.items()
    }


def resample_peak_curve(
    rt_array: list[float],
    intensity_array: list[float],
    n_points: int = 100,
) -> list[float]:
    if not rt_array or not intensity_array:
        return [0.0] * n_points
    pairs = sorted(zip(rt_array, intensity_array))
    xs = [item[0] for item in pairs]
    ys = [item[1] for item in pairs]
    if xs[-1] <= xs[0]:
        return [ys[0]] * n_points
    result: list[float] = []
    cursor = 0
    for idx in range(n_points):
        x = xs[0] + (xs[-1] - xs[0]) * idx / max(1, n_points - 1)
        while cursor < len(xs) - 2 and xs[cursor + 1] < x:
            cursor += 1
        x0, x1 = xs[cursor], xs[min(cursor + 1, len(xs) - 1)]
        y0, y1 = ys[cursor], ys[min(cursor + 1, len(ys) - 1)]
        if x1 <= x0:
            result.append(y0)
        else:
            ratio = (x - x0) / (x1 - x0)
            result.append(y0 * (1 - ratio) + y1 * ratio)
    return result


def compute_peak_chromatogram_similarity(
    peak_curves: dict[str, list[float]],
    areas: dict[str, float],
    widths: dict[str, float],
    local_alignment: dict[str, object],
) -> dict[str, float]:
    normalized: dict[str, list[float]] = {}
    for sample_id, curve in peak_curves.items():
        total = sum(curve) or 1.0
        normalized[sample_id] = [value / total for value in curve]
    if not normalized:
        return {"shape_score": 0.0, "area_score": 0.0, "apex_score": 0.0, "width_score": 0.0, "chromatogram_score": 0.0}
    consensus = [
        statistics.median(values)
        for values in zip(*normalized.values())
    ]
    shape_values = [cosine(curve, consensus) for curve in normalized.values()]
    local_shifts = [abs(float(value)) for value in dict(local_alignment.get("local_peak_shift_by_sample") or {}).values()]
    capped = any(dict(local_alignment.get("local_shift_capped_by_sample") or {}).values())
    max_allowed = float(local_alignment.get("max_local_shift_allowed") or 0.1)
    apex_score = max(0.0, 1.0 - (max(local_shifts, default=0.0) / max(max_allowed, 1e-9)))
    if capped:
        apex_score *= 0.5
    shape_score = statistics.mean(shape_values) if shape_values else 0.0
    area_score = score_from_cv(list(areas.values()))
    width_score = score_from_cv(list(widths.values()))
    chromatogram_score = 0.55 * shape_score + 0.25 * area_score + 0.10 * apex_score + 0.10 * width_score
    return {
        "shape_score": shape_score,
        "area_score": area_score,
        "apex_score": apex_score,
        "width_score": width_score,
        "chromatogram_score": chromatogram_score,
    }


def mz_tolerance(mz: float, params: PeakFirstParams) -> float:
    if params.mz_tolerance_mode == "ppm":
        return max(mz * params.mz_tolerance_ppm / 1_000_000.0, 1e-9)
    return params.mz_tolerance_da


def get_summed_spectrum_for_peak(
    scans: list[LCMSSpectrumScan],
    tic_peak: dict[str, object],
    rt_shift: float,
    local_shift: float,
    params: PeakFirstParams,
) -> list[tuple[float, float]]:
    start = float(tic_peak["rt_start"]) - params.peak_margin_min
    end = float(tic_peak["rt_end"]) + params.peak_margin_min
    values: list[tuple[float, float]] = []
    for scan in scans:
        rt = scan.rt + rt_shift + local_shift
        if start <= rt <= end:
            pairs = [
                (float(mz), float(intensity))
                for mz, intensity in zip(scan.mz_array, scan.intensity_array)
                if intensity > 0
            ]
            pairs.sort(key=lambda item: item[1], reverse=True)
            values.extend(pairs[: params.max_spectrum_points_per_scan])
    return values


def get_apex_spectrum_for_peak(
    scans: list[LCMSSpectrumScan],
    tic_peak: dict[str, object],
    rt_shift: float,
    local_shift: float,
    params: PeakFirstParams,
) -> list[tuple[float, float]]:
    apex = float(tic_peak["rt_apex"])
    if not scans:
        return []
    scan = min(scans, key=lambda item: abs((item.rt + rt_shift + local_shift) - apex))
    pairs = [
        (float(mz), float(intensity))
        for mz, intensity in zip(scan.mz_array, scan.intensity_array)
        if intensity > 0
    ]
    pairs.sort(key=lambda item: item[1], reverse=True)
    return pairs[: params.max_spectrum_points_per_scan]


def build_consensus_mz_bins(
    spectra: Iterable[list[tuple[float, float]]],
    params: PeakFirstParams,
) -> list[float]:
    points = sorted(
        (float(mz), float(intensity))
        for spectrum in spectra
        for mz, intensity in spectrum
        if intensity > 0
    )
    clusters: list[list[tuple[float, float]]] = []
    for mz, intensity in points:
        if not clusters:
            clusters.append([(mz, intensity)])
            continue
        center = weighted_mz(clusters[-1])
        if abs(mz - center) <= mz_tolerance(center, params):
            clusters[-1].append((mz, intensity))
        else:
            clusters.append([(mz, intensity)])
    bins = [weighted_mz(cluster) for cluster in clusters]
    return bins


def weighted_mz(points: list[tuple[float, float]]) -> float:
    total = sum(intensity for _, intensity in points)
    if total <= 0:
        return statistics.mean([mz for mz, _ in points]) if points else 0.0
    return sum(mz * intensity for mz, intensity in points) / total


def vectorize_spectrum(
    spectrum: list[tuple[float, float]],
    bins: list[float],
    params: PeakFirstParams,
) -> list[float]:
    vector = [0.0 for _ in bins]
    if not bins:
        return vector
    for mz, intensity in spectrum:
        if intensity <= 0:
            continue
        pos = bisect_left(bins, mz)
        candidates = [pos]
        if pos > 0:
            candidates.append(pos - 1)
        if pos + 1 < len(bins):
            candidates.append(pos + 1)
        best_index = None
        best_delta = float("inf")
        for index in candidates:
            if index < 0 or index >= len(bins):
                continue
            center = bins[index]
            delta = abs(mz - center)
            if delta < best_delta and delta <= mz_tolerance(center, params):
                best_index = index
                best_delta = delta
        if best_index is not None:
            vector[best_index] += intensity
    return vector


def normalize_spectrum_vector(
    spectrum_vector: list[float],
    method: str = "local_tic",
    transform: str = "log1p",
) -> list[float]:
    if method == "local_tic":
        factor = sum(spectrum_vector) or 1.0
        values = [value / factor * 1_000_000.0 for value in spectrum_vector]
    else:
        values = [float(value) for value in spectrum_vector]
    if transform == "log1p":
        return [math.log1p(max(0.0, value)) for value in values]
    return values


def reduce_mz_bins(
    bins: list[float],
    raw_matrix: dict[str, list[float]],
    params: PeakFirstParams,
) -> tuple[list[float], dict[str, list[float]]]:
    if not bins or not raw_matrix:
        return bins, raw_matrix
    totals = {
        index: sum(values[index] for values in raw_matrix.values())
        for index in range(len(bins))
    }
    max_total = max(totals.values(), default=0.0)
    keep = [
        index for index, total in totals.items()
        if max_total <= 0 or total >= max_total * params.min_relative_intensity
    ]
    keep.sort(key=lambda index: totals[index], reverse=True)
    keep = keep[: params.top_n_mz]
    keep.sort()
    return [bins[index] for index in keep], {
        sample_id: [values[index] for index in keep]
        for sample_id, values in raw_matrix.items()
    }


def compute_peak_spectrum_similarity(
    spectrum_matrix: dict[str, list[float]],
    apex_spectrum_matrix: dict[str, list[float]],
    params: PeakFirstParams,
) -> dict[str, float]:
    if not spectrum_matrix:
        return {
            "cosine_score": 0.0,
            "top_mz_overlap_score": 0.0,
            "relative_abundance_score": 0.0,
            "presence_absence_score": 0.0,
            "apex_spectrum_score": 0.0,
            "spectrum_score": 0.0,
        }
    normalized = {
        sample_id: normalize_spectrum_vector(vector)
        for sample_id, vector in spectrum_matrix.items()
    }
    consensus = [statistics.median(values) for values in zip(*normalized.values())] if normalized else []
    cosine_values = [cosine(vector, consensus) for vector in normalized.values()]
    cosine_score = statistics.mean(cosine_values) if cosine_values else 0.0
    top_sets: dict[str, set[int]] = {}
    top_n = min(20, max(1, len(next(iter(spectrum_matrix.values()), []))))
    for sample_id, vector in spectrum_matrix.items():
        top_sets[sample_id] = set(sorted(range(len(vector)), key=lambda i: vector[i], reverse=True)[:top_n])
    consensus_top = set().union(*top_sets.values()) if top_sets else set()
    top_overlap = [
        len(top_set & consensus_top) / max(1, len(consensus_top))
        for top_set in top_sets.values()
    ]
    top_mz_overlap_score = statistics.mean(top_overlap) if top_overlap else 0.0
    raw_values_by_bin = list(zip(*spectrum_matrix.values())) if spectrum_matrix else []
    common_cvs = []
    presence_mismatches = 0
    valid_bins = 0
    for values in raw_values_by_bin:
        positives = [value for value in values if value > 0]
        if not positives:
            continue
        valid_bins += 1
        if 0 < len(positives) < len(values):
            presence_mismatches += 1
        if len(positives) == len(values):
            mean = statistics.mean(values)
            if mean > 0:
                common_cvs.append(statistics.pstdev(values) / mean)
    abundance_cv = statistics.median(common_cvs) if common_cvs else 0.0
    relative_abundance_score = 1.0 / (1.0 + abundance_cv)
    presence_absence_score = 1.0 - (presence_mismatches / valid_bins if valid_bins else 0.0)
    apex_score = 0.0
    if apex_spectrum_matrix:
        apex_norm = {
            sample_id: normalize_spectrum_vector(vector)
            for sample_id, vector in apex_spectrum_matrix.items()
        }
        apex_consensus = [statistics.median(values) for values in zip(*apex_norm.values())] if apex_norm else []
        apex_values = [cosine(vector, apex_consensus) for vector in apex_norm.values()]
        apex_score = statistics.mean(apex_values) if apex_values else 0.0
    spectrum_score = (
        0.45 * cosine_score
        + 0.20 * top_mz_overlap_score
        + 0.15 * relative_abundance_score
        + 0.10 * presence_absence_score
        + 0.10 * apex_score
    )
    return {
        "cosine_score": cosine_score,
        "top_mz_overlap_score": top_mz_overlap_score,
        "relative_abundance_score": relative_abundance_score,
        "presence_absence_score": presence_absence_score,
        "apex_spectrum_score": apex_score,
        "spectrum_score": spectrum_score,
    }


def find_top_changed_mz(
    bins: list[float],
    raw_matrix: dict[str, list[float]],
    params: PeakFirstParams,
) -> list[dict[str, object]]:
    if not bins or not raw_matrix:
        return []
    normalized = {
        sample_id: normalize_spectrum_vector(vector)
        for sample_id, vector in raw_matrix.items()
    }
    changes: list[dict[str, object]] = []
    sample_ids = list(raw_matrix)
    for index, mz in enumerate(bins):
        raw_by_sample = {sample: raw_matrix[sample][index] for sample in sample_ids}
        norm_by_sample = {sample: normalized[sample][index] for sample in sample_ids}
        raw_values = list(raw_by_sample.values())
        norm_values = list(norm_by_sample.values())
        positive = [value for value in raw_values if value > 0]
        if not positive:
            continue
        mean_intensity = statistics.mean(norm_values)
        cv = statistics.pstdev(norm_values) / mean_intensity if mean_intensity > 0 else 0.0
        presence_count = len(positive)
        presence_bonus = 0.5 if 0 < presence_count < len(sample_ids) else 0.0
        fold = (max(norm_values) + 1e-12) / (min([value for value in norm_values if value > 0] or [1e-12]) + 1e-12)
        if presence_count < len(sample_ids):
            diff_type = "new_mz_signal" if presence_count <= len(sample_ids) / 2 else "missing_mz_signal"
        elif cv > 0.25:
            diff_type = "common_but_changed"
        else:
            diff_type = "unstable_low_confidence"
        ranking_score = (cv / (1 + cv)) * math.log1p(mean_intensity) + presence_bonus
        changes.append(
            {
                "mz": mz,
                "raw_intensity_by_sample": raw_by_sample,
                "normalized_intensity_by_sample": norm_by_sample,
                "presence_by_sample": {sample: raw_by_sample[sample] > 0 for sample in sample_ids},
                "presence_count": presence_count,
                "mean_intensity": mean_intensity,
                "cv": cv,
                "max_fold_change": fold,
                "difference_type": diff_type,
                "contribution_to_spectrum_difference": ranking_score,
                "ranking_score": ranking_score,
            }
        )
    changes.sort(key=lambda item: float(item["ranking_score"]), reverse=True)
    selected: list[dict[str, object]] = []
    min_gap = max(params.min_changed_mz_gap_da, params.mz_tolerance_da * 2.0)
    for change in changes:
        mz = float(change["mz"])
        if any(abs(mz - float(item["mz"])) < min_gap for item in selected):
            continue
        selected.append(change)
        if len(selected) >= params.top_n_changed_mz:
            break
    return selected


def extract_xic_points_for_peak(
    scans: list[LCMSSpectrumScan],
    tic_peak: dict[str, object],
    target_mz: float,
    rt_shift: float,
    local_shift: float,
    params: PeakFirstParams,
) -> list[tuple[float, float, float]]:
    margin = max(params.peak_margin_min, params.xic_feature_search_margin_min)
    start = float(tic_peak["rt_start"]) - margin
    end = float(tic_peak["rt_end"]) + margin
    tolerance = mz_tolerance(target_mz, params)
    points: list[tuple[float, float, float]] = []
    for scan in scans:
        aligned_rt = scan.rt + rt_shift + local_shift
        if aligned_rt < start or aligned_rt > end:
            continue
        total = 0.0
        mz_weighted = 0.0
        for mz, intensity in zip(scan.mz_array, scan.intensity_array):
            if abs(float(mz) - target_mz) <= tolerance:
                value = float(intensity)
                total += value
                mz_weighted += float(mz) * value
        observed_mz = mz_weighted / total if total > 0 else target_mz
        points.append((aligned_rt, total, observed_mz))
    points.sort(key=lambda item: item[0])
    return points


def detect_xic_lcms_feature(
    xic_points: list[tuple[float, float, float]],
    target_mz: float,
    sample_id: str,
    raw_file_id: str,
    tic_peak_id: str,
    params: PeakFirstParams,
    force_gap_fill: bool = False,
) -> dict[str, object]:
    if not xic_points:
        return {
            "feature_id": f"{tic_peak_id}:{sample_id}:{target_mz:.5f}:missing",
            "sample_id": sample_id,
            "raw_file_id": raw_file_id,
            "parent_tic_peak_id": tic_peak_id,
            "mz": target_mz,
            "observed_mz": target_mz,
            "rt_start": None,
            "rt_apex": None,
            "rt_end": None,
            "aligned_rt_apex": None,
            "area": 0.0,
            "height": 0.0,
            "signal_to_noise": 0.0,
            "match_status": "missing",
            "gap_filled": True,
        }
    rt = [point[0] for point in xic_points]
    intensity = [point[1] for point in xic_points]
    max_height = max(intensity, default=0.0)
    apex_index = intensity.index(max_height) if max_height > 0 else 0
    low_values = [value for value in intensity if value <= quantile(intensity, 0.4)]
    noise = max(statistics.median(low_values) if low_values else 0.0, quantile(intensity, 0.1), 1e-9)
    min_detect_height = max(noise * 3.0, max_height * params.feature_min_height_fraction)
    is_detected = max_height >= min_detect_height and not force_gap_fill
    boundary_level = max(noise, max_height * (params.feature_gap_fill_height_fraction if not is_detected else params.feature_min_height_fraction))
    left = apex_index
    while left > 0 and intensity[left] > boundary_level:
        left -= 1
    right = apex_index
    while right < len(intensity) - 1 and intensity[right] > boundary_level:
        right += 1
    if right <= left:
        left = max(0, apex_index - 1)
        right = min(len(intensity) - 1, apex_index + 1)
    area = 0.0
    for index in range(left + 1, right + 1):
        area += (intensity[index - 1] + intensity[index]) * 0.5 * max(rt[index] - rt[index - 1], 0.0)
    mz_points = [(point[2], point[1]) for point in xic_points[left : right + 1] if point[1] > 0]
    observed_mz = weighted_mz(mz_points) if mz_points else target_mz
    status = "matched" if is_detected else "gap_filled"
    return {
        "feature_id": f"{tic_peak_id}:{sample_id}:{target_mz:.5f}",
        "sample_id": sample_id,
        "raw_file_id": raw_file_id,
        "parent_tic_peak_id": tic_peak_id,
        "mz": observed_mz,
        "representative_mz": target_mz,
        "mz_tolerance": mz_tolerance(target_mz, params),
        "rt_start": rt[left],
        "rt_apex": rt[apex_index],
        "rt_end": rt[right],
        "aligned_rt_apex": rt[apex_index],
        "area": area,
        "height": max_height,
        "signal_to_noise": max_height / noise if noise > 0 else 0.0,
        "match_status": status,
        "gap_filled": status == "gap_filled",
    }


def mzmine_join_match_score(
    feature: dict[str, object],
    representative_mz: float,
    representative_rt: float,
    params: PeakFirstParams,
) -> float:
    mz_tol = max(mz_tolerance(representative_mz, params), 1e-9)
    rt_tol = max(params.feature_rt_tolerance_min, 1e-9)
    mz_delta = abs(float(feature.get("mz") or representative_mz) - representative_mz)
    rt_value = feature.get("aligned_rt_apex")
    rt_delta = abs(float(rt_value) - representative_rt) if rt_value is not None else rt_tol
    mz_score = max(0.0, 1.0 - mz_delta / mz_tol)
    rt_score = max(0.0, 1.0 - rt_delta / rt_tol)
    weight_sum = max(params.feature_mz_weight + params.feature_rt_weight, 1e-9)
    return (params.feature_mz_weight * mz_score + params.feature_rt_weight * rt_score) / weight_sum


def classify_feature_group_difference(
    areas: dict[str, float],
    matched_count: int,
    sample_count: int,
    similarity_score: float,
) -> str:
    positive = [value for value in areas.values() if value > 0]
    if matched_count == 0 or not positive:
        return "low_confidence"
    if matched_count < sample_count:
        return "presence_absence"
    if similarity_score < 0.60:
        return "area_changed"
    if similarity_score < 0.80:
        return "moderate_difference"
    return "common_feature"


def build_peak_feature_groups(
    candidate_mz_items: list[dict[str, object]],
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    raw_file_by_sample: dict[str, LCMSRawFile],
    tic_peak: dict[str, object],
    global_shifts: dict[str, float],
    local_shifts: dict[str, float],
    params: PeakFirstParams,
) -> list[dict[str, object]]:
    groups: list[dict[str, object]] = []
    sample_ids = list(scans_by_sample)
    tic_peak_id = str(tic_peak.get("tic_peak_id") or "TICP")
    representative_rt = float(tic_peak.get("rt_apex") or 0.0)
    for rank, item in enumerate(candidate_mz_items, start=1):
        target_mz = float(item.get("mz") or 0.0)
        features: dict[str, dict[str, object]] = {}
        match_scores: dict[str, float] = {}
        for sample_id, scans in scans_by_sample.items():
            raw_file = raw_file_by_sample.get(sample_id)
            xic = extract_xic_points_for_peak(
                scans,
                tic_peak,
                target_mz,
                global_shifts.get(sample_id, 0.0),
                float(local_shifts.get(sample_id, 0.0)),
                params,
            )
            feature = detect_xic_lcms_feature(
                xic,
                target_mz,
                sample_id,
                raw_file.raw_file_id if raw_file else "",
                tic_peak_id,
                params,
            )
            features[sample_id] = feature
        detected_features = [
            feature for feature in features.values()
            if feature.get("match_status") == "matched" and feature.get("aligned_rt_apex") is not None
        ]
        if detected_features:
            representative_mz = weighted_mz([
                (float(feature.get("mz") or target_mz), float(feature.get("height") or 0.0))
                for feature in detected_features
            ])
            representative_rt = statistics.median(float(feature["aligned_rt_apex"]) for feature in detected_features)
        else:
            representative_mz = target_mz
        for sample_id, feature in features.items():
            match_scores[sample_id] = mzmine_join_match_score(feature, representative_mz, representative_rt, params)
            mz_delta = abs(float(feature.get("mz") or representative_mz) - representative_mz)
            rt_value = feature.get("aligned_rt_apex")
            rt_delta = abs(float(rt_value) - representative_rt) if rt_value is not None else params.feature_drift_tolerance_min
            rt_shift_only = mz_delta <= mz_tolerance(representative_mz, params) and rt_delta <= params.feature_drift_tolerance_min
            if rt_shift_only:
                feature["rt_shift_corrected_match"] = True
                feature["feature_rt_delta"] = rt_delta
            if match_scores[sample_id] < 0.30 and feature.get("match_status") == "matched" and not rt_shift_only:
                feature["match_status"] = "low_score"
        areas = {sample_id: float(feature.get("area") or 0.0) for sample_id, feature in features.items()}
        heights = {sample_id: float(feature.get("height") or 0.0) for sample_id, feature in features.items()}
        max_area = max(areas.values(), default=0.0)
        normalized = normalize_spectrum_vector([areas[sample_id] for sample_id in sample_ids], method="local_tic", transform="log1p")
        mean_norm = statistics.mean(normalized) if normalized else 0.0
        cv = statistics.pstdev(normalized) / mean_norm if mean_norm > 0 and len(normalized) > 1 else 0.0
        similarity = 1.0 / (1.0 + cv)
        difference = 1.0 - similarity
        matched_count = sum(1 for feature in features.values() if feature.get("match_status") == "matched")
        positive_norm = [value for value in normalized if value > 0]
        fold = (max(normalized) + 1e-12) / (min(positive_norm or [1e-12]) + 1e-12) if normalized else 1.0
        difference_type = classify_feature_group_difference(areas, matched_count, len(sample_ids), similarity)
        abundance_score = math.log1p(max_area)
        abundance_weighted_score = difference * (1.0 + params.feature_abundance_weight * abundance_score)
        groups.append(
            {
                "feature_group_id": f"{tic_peak_id}_FG_{rank:04d}",
                "parent_tic_peak_id": tic_peak_id,
                "representative_mz": representative_mz,
                "representative_rt": representative_rt,
                "sample_count": matched_count,
                "missing_sample_count": len(sample_ids) - matched_count,
                "cv": cv,
                "max_fold_change": fold,
                "max_area": max_area,
                "abundance_score": abundance_score,
                "abundance_weighted_score": abundance_weighted_score,
                "similarity_score": similarity,
                "difference_score": difference,
                "difference_type": difference_type,
                "match_score_by_sample": match_scores,
                "area_by_sample": areas,
                "height_by_sample": heights,
                "features_by_sample": features,
                "rt_correction_by_sample": {
                    sample_id: float(local_shifts.get(sample_id, 0.0))
                    for sample_id in sample_ids
                },
                "alignment_method": "mzmine_join_aligner_with_tic_peak_local_rt_correction",
            }
        )
    groups.sort(
        key=lambda group: (
            float(group.get("abundance_weighted_score") or 0.0),
            float(group.get("difference_score") or 0.0),
            float(group.get("max_fold_change") or 0.0),
        ),
        reverse=True,
    )
    return groups


def feature_groups_to_top_changed_mz(groups: list[dict[str, object]], params: PeakFirstParams) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for group in groups:
        rows.append(
            {
                "mz": group["representative_mz"],
                "feature_group_id": group["feature_group_id"],
                "representative_rt": group["representative_rt"],
                "raw_intensity_by_sample": group["area_by_sample"],
                "normalized_intensity_by_sample": group["area_by_sample"],
                "presence_by_sample": {
                    sample_id: float(area) > 0 and not group["features_by_sample"][sample_id].get("gap_filled")
                    for sample_id, area in group["area_by_sample"].items()
                },
                "presence_count": group["sample_count"],
                "mean_intensity": statistics.mean(list(group["area_by_sample"].values())) if group["area_by_sample"] else 0.0,
                "cv": group["cv"],
                "max_fold_change": group["max_fold_change"],
                "difference_type": group["difference_type"],
                "contribution_to_spectrum_difference": group["difference_score"],
                "ranking_score": group.get("abundance_weighted_score", group["difference_score"]),
                "match_score_by_sample": group["match_score_by_sample"],
            }
        )
    return rows[: params.top_n_changed_mz]


def global_feature_groups_from_peaks(
    peak_results: list[dict[str, object]],
    params: PeakFirstParams,
    limit: int = 500,
) -> list[dict[str, object]]:
    raw_rows: list[dict[str, object]] = []
    max_abundance = 0.0
    for peak in peak_results:
        for group in peak.get("feature_groups") or []:
            max_abundance = max(max_abundance, float(group.get("abundance_score") or 0.0))
    for peak in peak_results:
        for group in peak.get("feature_groups") or []:
            abundance_norm = (float(group.get("abundance_score") or 0.0) / max_abundance) if max_abundance > 0 else 0.0
            ranking = float(group.get("difference_score") or 0.0) * (1.0 + params.feature_abundance_weight * abundance_norm)
            row = {
                "feature_group_id": group.get("feature_group_id"),
                "parent_tic_peak_id": group.get("parent_tic_peak_id"),
                "representative_mz": group.get("representative_mz"),
                "representative_rt": group.get("representative_rt"),
                "similarity_score": group.get("similarity_score"),
                "difference_score": group.get("difference_score"),
                "ranking_score": ranking,
                "abundance_norm": abundance_norm,
                "max_area": group.get("max_area"),
                "max_fold_change": group.get("max_fold_change"),
                "difference_type": group.get("difference_type"),
                "sample_count": group.get("sample_count"),
                "missing_sample_count": group.get("missing_sample_count"),
                "area_by_sample": group.get("area_by_sample"),
                "match_score_by_sample": group.get("match_score_by_sample"),
                "rt_correction_by_sample": group.get("rt_correction_by_sample"),
                "source_feature_group_ids": [group.get("feature_group_id")],
                "source_parent_tic_peak_ids": [group.get("parent_tic_peak_id")],
            }
            raw_rows.append(row)
    raw_rows.sort(
        key=lambda item: (
            float(item.get("ranking_score") or 0.0),
            float(item.get("difference_score") or 0.0),
            float(item.get("max_area") or 0.0),
        ),
        reverse=True,
    )
    clusters: list[dict[str, object]] = []
    for row in raw_rows:
        mz = float(row.get("representative_mz") or 0.0)
        rt = float(row.get("representative_rt") or 0.0)
        mz_tol = max(mz_tolerance(mz, params) * params.global_feature_mz_merge_tolerance_factor, 1e-9)
        rt_tol = max(params.global_feature_rt_merge_tolerance_min, params.feature_rt_tolerance_min)
        matched_cluster: dict[str, object] | None = None
        for cluster in clusters:
            if abs(mz - float(cluster.get("representative_mz") or mz)) <= mz_tol and abs(rt - float(cluster.get("representative_rt") or rt)) <= rt_tol:
                matched_cluster = cluster
                break
        if matched_cluster is None:
            cluster = dict(row)
            cluster["merged_feature_count"] = 1
            cluster["merged_feature_group_ids"] = list(row["source_feature_group_ids"])
            cluster["merged_parent_tic_peak_ids"] = list(row["source_parent_tic_peak_ids"])
            clusters.append(cluster)
            continue
        matched_cluster["merged_feature_count"] = int(matched_cluster.get("merged_feature_count") or 1) + 1
        merged_ids = set(matched_cluster.get("merged_feature_group_ids") or [])
        merged_ids.update(row.get("source_feature_group_ids") or [])
        matched_cluster["merged_feature_group_ids"] = sorted(str(item) for item in merged_ids if item)
        parent_ids = set(matched_cluster.get("merged_parent_tic_peak_ids") or [])
        parent_ids.update(row.get("source_parent_tic_peak_ids") or [])
        matched_cluster["merged_parent_tic_peak_ids"] = sorted(str(item) for item in parent_ids if item)
        if float(row.get("ranking_score") or 0.0) > float(matched_cluster.get("ranking_score") or 0.0):
            keep_ids = matched_cluster["merged_feature_group_ids"]
            keep_parents = matched_cluster["merged_parent_tic_peak_ids"]
            keep_count = matched_cluster["merged_feature_count"]
            matched_cluster.update(row)
            matched_cluster["merged_feature_group_ids"] = keep_ids
            matched_cluster["merged_parent_tic_peak_ids"] = keep_parents
            matched_cluster["merged_feature_count"] = keep_count
        else:
            matched_cluster["ranking_score"] = max(float(matched_cluster.get("ranking_score") or 0.0), float(row.get("ranking_score") or 0.0))
            matched_cluster["difference_score"] = max(float(matched_cluster.get("difference_score") or 0.0), float(row.get("difference_score") or 0.0))
            matched_cluster["max_area"] = max(float(matched_cluster.get("max_area") or 0.0), float(row.get("max_area") or 0.0))
            matched_cluster["max_fold_change"] = max(float(matched_cluster.get("max_fold_change") or 0.0), float(row.get("max_fold_change") or 0.0))
    rows = clusters
    rows.sort(
        key=lambda item: (
            float(item.get("ranking_score") or 0.0),
            float(item.get("difference_score") or 0.0),
            float(item.get("max_area") or 0.0),
        ),
        reverse=True,
    )
    return rows[:limit]


def classify_tic_peak_status(
    chromatogram_scores: dict[str, float],
    spectrum_scores: dict[str, float],
    local_alignment: dict[str, object],
) -> str:
    chromatogram_score = chromatogram_scores.get("chromatogram_score", 0.0)
    spectrum_score = spectrum_scores.get("spectrum_score", 0.0)
    presence_score = spectrum_scores.get("presence_absence_score", 0.0)
    top_overlap = spectrum_scores.get("top_mz_overlap_score", 0.0)
    relative_score = spectrum_scores.get("relative_abundance_score", 0.0)
    if any(dict(local_alignment.get("local_shift_capped_by_sample") or {}).values()):
        return "possible_rt_shift"
    if chromatogram_score >= 0.90 and spectrum_score >= 0.90 and presence_score >= 0.90:
        return "high_consistency"
    if chromatogram_score >= 0.90 and spectrum_score < 0.85:
        return "chromatogram_consistent_spectrum_changed"
    if chromatogram_score < 0.85:
        return "chromatogram_changed"
    if spectrum_score < 0.85 or top_overlap < 0.80 or presence_score < 0.85 or relative_score < 0.80:
        return "need_local_rt_mz_analysis"
    return "moderate_consistency"


def integrate_curve(points: list[tuple[float, float]]) -> float:
    if len(points) < 2:
        return 0.0
    ordered = sorted(points)
    return sum(
        (ordered[i - 1][1] + ordered[i][1]) * 0.5 * max(ordered[i][0] - ordered[i - 1][0], 0.0)
        for i in range(1, len(ordered))
    )


def score_tic_peak_for_selection(
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    tic_peak: dict[str, object],
    shifts: dict[str, float],
    params: PeakFirstParams,
) -> dict[str, object]:
    areas: dict[str, float] = {}
    heights: dict[str, float] = {}
    for sample_id, scans in scans_by_sample.items():
        window = extract_peak_window(scans, tic_peak, shifts.get(sample_id, 0.0), True, params.peak_margin_min, "tic")
        values = [value for _, value in window]
        areas[sample_id] = integrate_curve(window)
        heights[sample_id] = max(values, default=0.0)
    max_area = max(areas.values(), default=0.0)
    positive = [value for value in areas.values() if value > 0]
    if positive:
        min_positive = min(positive)
        fold = (max_area + 1e-12) / (min_positive + 1e-12)
        mean = statistics.mean(positive)
        cv = statistics.pstdev(positive) / mean if len(positive) > 1 and mean > 0 else 0.0
    else:
        fold = 1.0
        cv = 0.0
    presence_count = sum(1 for value in areas.values() if value >= max_area * 0.02 and value > 0)
    similarity = 1.0 / (1.0 + cv)
    difference = 1.0 - similarity
    if presence_count < len(areas) and max_area > 0:
        difference = max(difference, 0.5)
    ranking = math.log1p(max_area) * (1.0 + params.tic_difference_selection_weight * difference)
    return {
        "tic_area_by_sample": areas,
        "tic_height_by_sample": heights,
        "tic_presence_count": presence_count,
        "tic_area_cv": cv,
        "tic_area_fold": fold,
        "tic_difference_score": difference,
        "tic_selection_score": ranking,
    }


def prepare_peak_first_payload(
    raw_files: list[LCMSRawFile],
    scans_by_sample: dict[str, list[LCMSSpectrumScan]],
    project_id: str,
    reference_sample: str,
    params: PeakFirstParams | None = None,
) -> dict[str, object]:
    params = params or PeakFirstParams()
    alignment = align_rt_by_main_or_feature_peak(scans_by_sample, reference_sample, signal="bpc")
    shifts = {str(k): float(v) for k, v in dict(alignment["rt_shift_by_sample"]).items()}
    consensus_curve = consensus_tic_curve(
        scans_by_sample,
        shifts,
        "tic",
        max_points=params.consensus_max_points,
        target_step_min=params.consensus_rt_step_min,
    )
    consensus_peaks = filter_tic_peaks(refine_tic_peak_boundaries(detect_tic_peaks(consensus_curve, params), params), params)
    sample_peaks = sample_union_tic_peaks(scans_by_sample, shifts, params, "tic")
    tic_peaks = filter_tic_peaks(merge_tic_peak_candidates(consensus_peaks, sample_peaks, params), params)
    for peak in tic_peaks:
        peak.update(score_tic_peak_for_selection(scans_by_sample, peak, shifts, params))
    confirmed = sorted(
        [peak for peak in tic_peaks if peak["status"] == "confirmed_peak"],
        key=lambda item: (
            float(item.get("tic_selection_score") or 0.0),
            float(item.get("area") or 0.0),
        ),
        reverse=True,
    )[: params.top_n_peaks]
    peak_results: list[dict[str, object]] = []
    raw_file_by_sample = {raw_file.sample_id: raw_file for raw_file in raw_files}
    for peak in confirmed:
        local = local_align_peak_by_apex(scans_by_sample, peak, shifts, params)
        local_shifts = dict(local["local_peak_shift_by_sample"])
        cut_bounds = sample_specific_cut_bounds(peak, local_shifts, params)
        curves: dict[str, list[float]] = {}
        areas: dict[str, float] = {}
        widths: dict[str, float] = {}
        for sample_id, scans in scans_by_sample.items():
            window = extract_peak_window(scans, peak, shifts.get(sample_id, 0.0) + float(local_shifts.get(sample_id, 0.0)), True, params.peak_margin_min)
            curves[sample_id] = resample_peak_curve([rt for rt, _ in window], [value for _, value in window], params.resample_points)
            areas[sample_id] = integrate_curve(window)
            widths[sample_id] = float(peak.get("width") or 0.0)
        chrom_scores = compute_peak_chromatogram_similarity(curves, areas, widths, local)
        summed_spectra = {
            sample_id: get_summed_spectrum_for_peak(
                scans,
                peak,
                shifts.get(sample_id, 0.0),
                float(local_shifts.get(sample_id, 0.0)),
                params,
            )
            for sample_id, scans in scans_by_sample.items()
        }
        bins = build_consensus_mz_bins(summed_spectra.values(), params)
        raw_matrix = {
            sample_id: vectorize_spectrum(spectrum, bins, params)
            for sample_id, spectrum in summed_spectra.items()
        }
        bins, raw_matrix = reduce_mz_bins(bins, raw_matrix, params)
        apex_spectra = {
            sample_id: get_apex_spectrum_for_peak(
                scans,
                peak,
                shifts.get(sample_id, 0.0),
                float(local_shifts.get(sample_id, 0.0)),
                params,
            )
            for sample_id, scans in scans_by_sample.items()
        }
        apex_matrix = {
            sample_id: vectorize_spectrum(spectrum, bins, params)
            for sample_id, spectrum in apex_spectra.items()
        }
        spectrum_scores = compute_peak_spectrum_similarity(raw_matrix, apex_matrix, params)
        spectrum_changed = find_top_changed_mz(bins, raw_matrix, params)
        feature_groups = build_peak_feature_groups(
            spectrum_changed,
            scans_by_sample,
            raw_file_by_sample,
            peak,
            shifts,
            local_shifts,
            params,
        )
        top_changed = feature_groups_to_top_changed_mz(feature_groups, params) or spectrum_changed
        status = classify_tic_peak_status(chrom_scores, spectrum_scores, local)
        consistency = (
            0.15 * chrom_scores["shape_score"]
            + 0.10 * chrom_scores["area_score"]
            + 0.05 * chrom_scores["apex_score"]
            + 0.05 * chrom_scores["width_score"]
            + 0.65 * spectrum_scores["spectrum_score"]
        )
        peak_results.append(
            {
                **peak,
                **chrom_scores,
                **spectrum_scores,
                "peak_consistency_score": consistency,
                "status": status,
                "local_alignment": {k: v for k, v in local.items() if k != "windows"},
                "sample_cut_bounds": cut_bounds,
                "areas_by_sample": areas,
                "top_changed_mz": top_changed,
                "spectrum_top_changed_mz": spectrum_changed,
                "feature_groups": feature_groups,
                "feature_alignment_method": "mzmine_join_aligner_with_tic_peak_local_rt_correction",
                "spectrum_mz_bins": bins,
                "spectrum_raw_matrix": raw_matrix,
            }
        )
    peak_results.sort(key=lambda item: float(item.get("rt_apex") or 0.0))
    global_feature_groups = global_feature_groups_from_peaks(peak_results, params)
    return {
        "module": "LCMSPeakFirstCompare",
        "project_id": project_id,
        "raw_files": [raw_file_payload(raw_file) for raw_file in raw_files],
        "sample_ids": list(scans_by_sample),
        "reference_sample": reference_sample,
        "alignment": alignment,
        "params": params.__dict__,
        "chromatograms": {
            sample_id: {
                "raw": chromatogram_points(scans, 0.0),
                "aligned": chromatogram_points(scans, shifts.get(sample_id, 0.0)),
            }
            for sample_id, scans in scans_by_sample.items()
        },
        "tic_peaks": tic_peaks,
        "peak_results": peak_results,
        "global_feature_groups": global_feature_groups,
    }
