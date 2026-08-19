#!/usr/bin/env python3
"""Feature matrix construction and difference classification."""

from __future__ import annotations

import math
import statistics

from .lcms_feature_matching import feature_rt
from .lcms_models import LCMSFeature, LCMSFeatureGroup


def coefficient_of_variation(values: list[float]) -> float | None:
    positives = [value for value in values if value > 0]
    if len(positives) < 2:
        return None
    mean = statistics.mean(positives)
    return statistics.pstdev(positives) / mean if mean > 0 else None


def fold_change(target: float, reference: float, missing_value: float = 1.0) -> float | None:
    if target <= 0 and reference <= 0:
        return None
    return max(target, missing_value) / max(reference, missing_value)


def classify_group(
    area_by_sample: dict[str, float],
    reference_samples: set[str],
    target_samples: set[str],
    fold_change_threshold: float = 2.0,
    missing_threshold: float = 1.0,
    rt_shift_flag: bool = False,
) -> tuple[str, float, float | None]:
    ref_values = [area_by_sample.get(sample, 0.0) for sample in reference_samples]
    target_values = [area_by_sample.get(sample, 0.0) for sample in target_samples]
    ref_present = any(value > missing_threshold for value in ref_values)
    target_present = any(value > missing_threshold for value in target_values)
    ref_mean = statistics.mean(ref_values) if ref_values else 0.0
    target_mean = statistics.mean(target_values) if target_values else 0.0
    fc = fold_change(target_mean, ref_mean)
    if rt_shift_flag and ref_present and target_present:
        return "rt_shift", 0.75, fc
    if target_present and not ref_present:
        return "new_feature", 0.9, fc
    if ref_present and not target_present:
        return "missing_feature", 0.9, fc
    if fc is not None and fc >= fold_change_threshold:
        return "increased", min(0.95, 0.5 + math.log2(fc) / 6.0), fc
    if fc is not None and fc <= 1.0 / fold_change_threshold:
        return "decreased", min(0.95, 0.5 + abs(math.log2(fc)) / 6.0), fc
    if ref_present and target_present:
        return "common_feature", 0.7, fc
    return "uncertain", 0.3, fc


def build_feature_groups(
    matched_groups: list[list[LCMSFeature]],
    sample_ids: list[str],
    reference_samples: set[str],
    target_samples: set[str],
    rt_shift_threshold_min: float = 0.2,
    fold_change_threshold: float = 2.0,
) -> list[LCMSFeatureGroup]:
    groups: list[LCMSFeatureGroup] = []
    for index, members in enumerate(matched_groups, start=1):
        group_id = members[0].feature_group_id or f"FG{index:04d}"
        representative_mz = sum(member.mz for member in members) / len(members)
        representative_rt = sum(feature_rt(member) for member in members) / len(members)
        area_by_sample = {sample_id: 0.0 for sample_id in sample_ids}
        sample_presence = {sample_id: False for sample_id in sample_ids}
        for member in members:
            area_by_sample[member.sample_id] = max(area_by_sample.get(member.sample_id, 0.0), member.area)
            sample_presence[member.sample_id] = True
        ref_rts = [feature_rt(member) for member in members if member.sample_id in reference_samples]
        target_rts = [feature_rt(member) for member in members if member.sample_id in target_samples]
        rt_shift_flag = bool(ref_rts and target_rts and abs(statistics.mean(target_rts) - statistics.mean(ref_rts)) > rt_shift_threshold_min)
        difference_type, confidence, fc = classify_group(
            area_by_sample,
            reference_samples,
            target_samples,
            fold_change_threshold=fold_change_threshold,
            rt_shift_flag=rt_shift_flag,
        )
        positives = [area for area in area_by_sample.values() if area > 0]
        max_fc = None
        if positives:
            max_fc = max(positives) / max(min(positives), 1.0)
        groups.append(
            LCMSFeatureGroup(
                feature_group_id=group_id,
                representative_mz=representative_mz,
                representative_rt=representative_rt,
                sample_count=sum(1 for present in sample_presence.values() if present),
                missing_sample_count=sum(1 for present in sample_presence.values() if not present),
                cv=coefficient_of_variation(list(area_by_sample.values())),
                max_fold_change=max_fc if max_fc is not None else fc,
                difference_type=difference_type,
                confidence_score=confidence,
                sample_presence=sample_presence,
                area_by_sample=area_by_sample,
            )
        )
    return groups

