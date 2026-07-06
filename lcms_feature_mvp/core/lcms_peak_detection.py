#!/usr/bin/env python3
"""Lightweight XIC peak detection and trapezoid integration."""

from __future__ import annotations

import statistics

from .lcms_models import LCMSFeature, XICTrace


def moving_average(values: list[float], window: int) -> list[float]:
    if window <= 1 or len(values) <= 2:
        return list(values)
    if window % 2 == 0:
        window += 1
    radius = window // 2
    padded = [values[0]] * radius + values + [values[-1]] * radius
    result: list[float] = []
    for index in range(radius, radius + len(values)):
        segment = padded[index - radius:index + radius + 1]
        result.append(sum(segment) / len(segment))
    return result


def trapezoid_area(rt: list[float], signal: list[float]) -> float:
    area = 0.0
    for t0, t1, y0, y1 in zip(rt, rt[1:], signal, signal[1:]):
        area += max(0.0, t1 - t0) * (max(0.0, y0) + max(0.0, y1)) / 2.0
    return area


def robust_noise(values: list[float]) -> float:
    if not values:
        return 1.0
    ordered = sorted(values)
    low_count = max(3, len(values) // 5)
    low_values = ordered[:low_count]
    if len(low_values) >= 2:
        spread = statistics.pstdev(low_values)
    else:
        spread = low_values[0] if low_values else 1.0
    return max(spread, statistics.median(low_values) * 0.1 if low_values else 0.0, 1.0)


def detect_xic_peaks(
    xic: XICTrace,
    smoothing_window: int = 5,
    min_peak_height: float = 2000.0,
    min_peak_area: float = 100.0,
    min_snr: float = 3.0,
    min_peak_width: float = 0.03,
    max_peak_width: float = 1.0,
    relative_foot_fraction: float = 0.05,
) -> list[LCMSFeature]:
    rt = xic.rt_array
    raw = xic.intensity_array
    if len(rt) < 3:
        return []
    signal = moving_average(raw, smoothing_window)
    noise = robust_noise(signal)
    candidates: list[int] = []
    for index in range(1, len(signal) - 1):
        if signal[index] >= signal[index - 1] and signal[index] > signal[index + 1]:
            if signal[index] >= min_peak_height and signal[index] / noise >= min_snr:
                candidates.append(index)

    features: list[LCMSFeature] = []
    occupied: list[tuple[int, int]] = []
    for apex_idx in sorted(candidates, key=lambda idx: signal[idx], reverse=True):
        foot_level = max(noise * 2.0, signal[apex_idx] * relative_foot_fraction)
        left = apex_idx
        while left > 0 and signal[left] > foot_level:
            left -= 1
        right = apex_idx
        while right < len(signal) - 1 and signal[right] > foot_level:
            right += 1
        if any(not (right < used_left or left > used_right) for used_left, used_right in occupied):
            continue
        width = rt[right] - rt[left]
        if width < min_peak_width or width > max_peak_width:
            continue
        segment_rt = rt[left:right + 1]
        segment_raw = raw[left:right + 1]
        baseline = min(segment_raw[0], segment_raw[-1])
        corrected = [max(0.0, value - baseline) for value in segment_raw]
        area = trapezoid_area(segment_rt, corrected)
        height = max(corrected, default=0.0)
        if area < min_peak_area:
            continue
        feature = LCMSFeature(
            feature_id=f"{xic.sample_id}:{xic.target_mz:.5f}:{rt[apex_idx]:.4f}",
            sample_id=xic.sample_id,
            raw_file_id=xic.raw_file_id,
            mz=xic.target_mz,
            mz_tolerance=xic.mz_tolerance,
            rt_start=rt[left],
            rt_apex=rt[apex_idx],
            rt_end=rt[right],
            aligned_rt_apex=None,
            area=area,
            height=height,
            signal_to_noise=height / noise,
        )
        features.append(feature)
        occupied.append((left, right))
    return sorted(features, key=lambda item: item.area, reverse=True)

