#!/usr/bin/env python3
"""Cross-sample LC-MS feature matching."""

from __future__ import annotations

from .lcms_models import LCMSFeature
from .lcms_xic import ppm_window


def feature_rt(feature: LCMSFeature) -> float:
    return feature.aligned_rt_apex if feature.aligned_rt_apex is not None else feature.rt_apex


def match_features(
    features: list[LCMSFeature],
    sample_ids: list[str],
    mz_tolerance_ppm: float = 20.0,
    rt_tolerance_min: float = 0.2,
    min_samples_present: int = 1,
) -> list[list[LCMSFeature]]:
    """Greedy grouping by representative m/z and apex RT."""
    groups: list[list[LCMSFeature]] = []
    for feature in sorted(features, key=lambda item: (item.mz, feature_rt(item), item.sample_id)):
        best_index: int | None = None
        best_score = float("inf")
        for index, group in enumerate(groups):
            if any(member.sample_id == feature.sample_id for member in group):
                continue
            rep_mz = sum(member.mz for member in group) / len(group)
            rep_rt = sum(feature_rt(member) for member in group) / len(group)
            mz_delta = abs(feature.mz - rep_mz)
            rt_delta = abs(feature_rt(feature) - rep_rt)
            if mz_delta <= ppm_window(rep_mz, mz_tolerance_ppm) and rt_delta <= rt_tolerance_min:
                score = mz_delta / max(ppm_window(rep_mz, mz_tolerance_ppm), 1e-12) + rt_delta / rt_tolerance_min
                if score < best_score:
                    best_index = index
                    best_score = score
        if best_index is None:
            groups.append([feature])
        else:
            groups[best_index].append(feature)

    groups = [group for group in groups if len({feature.sample_id for feature in group}) >= min_samples_present]
    groups.sort(key=lambda group: (sum(feature_rt(item) for item in group) / len(group), sum(item.mz for item in group) / len(group)))
    for index, group in enumerate(groups, start=1):
        group_id = f"FG{index:04d}"
        for feature in group:
            feature.feature_group_id = group_id
            feature.match_status = "matched" if len(group) > 1 else "single_sample"
    return groups

