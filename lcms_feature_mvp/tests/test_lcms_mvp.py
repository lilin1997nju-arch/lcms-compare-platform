from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from core.lcms_difference import build_feature_groups  # noqa: E402
from core.lcms_feature_matching import match_features  # noqa: E402
from core.lcms_models import LCMSSpectrumScan  # noqa: E402
from core.lcms_parser import mock_scans_from_raw  # noqa: E402
from core.lcms_peak_detection import detect_xic_peaks  # noqa: E402
from core.lcms_peak_first import (  # noqa: E402
    PeakFirstParams,
    abundance_profile_metrics,
    build_consensus_mz_bins,
    detect_xic_lcms_feature,
    detect_tic_peaks,
    extract_xic_points_for_peak,
    filter_tic_peaks,
    find_top_changed_mz,
    global_feature_groups_from_peaks,
    local_align_peak_by_apex,
    prepare_peak_first_payload,
)
from core.lcms_workbench import build_similarity_heatmaps, main_peak_alignment, sort_difference_regions, spectrum_payload  # noqa: E402
from core.lcms_xic import extract_xic, screen_candidate_mz  # noqa: E402
from run_xcalibur_workbench import (  # noqa: E402
    auto_feature_regions_from_payload,
    feature_matrix_header,
    feature_matrix_rows_from_regions,
    write_workbench_sqlite,
)
from run_peak_first_compare import PEAK_FIRST_TEMPLATE  # noqa: E402
from serve_xcalibur_workbench import (  # noqa: E402
    build_feature_matrix,
    feature_matrix_csv,
    heatmap_window,
    read_saved_features,
    replace_saved_features,
)
from serve_peak_first_compare import single_centroid_intensity  # noqa: E402


def gaussian(x: float, center: float, width: float, height: float) -> float:
    return height * pow(2.718281828, -0.5 * pow((x - center) / width, 2))


def synthetic_chrom_scans(sample_id: str, main_rt: float, impurity_rt: float, impurity_height: float) -> list[LCMSSpectrumScan]:
    scans: list[LCMSSpectrumScan] = []
    for index in range(121):
        rt = index * 0.1
        main = gaussian(rt, main_rt, 0.18, 1000.0)
        impurity = gaussian(rt, impurity_rt, 0.18, impurity_height)
        intensity = main + impurity + 5.0
        scans.append(
            LCMSSpectrumScan(
                scan_id=f"{sample_id}_{index}",
                raw_file_id=f"{sample_id}.raw",
                sample_id=sample_id,
                rt=rt,
                ms_level=1,
                mz_array=[500.0],
                intensity_array=[intensity],
                tic=intensity,
                base_peak_mz=500.0,
                base_peak_intensity=intensity,
            )
        )
    return scans


def scaled_single_mz_scans(sample_id: str, scale: float) -> list[LCMSSpectrumScan]:
    scans: list[LCMSSpectrumScan] = []
    for index in range(21):
        rt = index * 0.1
        intensity = gaussian(rt, 1.0, 0.2, 1000.0) * scale
        scans.append(
            LCMSSpectrumScan(
                scan_id=f"{sample_id}_{index}",
                raw_file_id=f"{sample_id}.raw",
                sample_id=sample_id,
                rt=rt,
                ms_level=1,
                mz_array=[500.0],
                intensity_array=[intensity],
                tic=intensity,
                base_peak_mz=500.0,
                base_peak_intensity=intensity,
            )
        )
    return scans


def shifted_single_feature_scans(sample_id: str, rt_center: float, mz: float) -> list[LCMSSpectrumScan]:
    scans: list[LCMSSpectrumScan] = []
    for index in range(21):
        rt = index * 0.1
        intensity = gaussian(rt, rt_center, 0.12, 1000.0)
        scans.append(
            LCMSSpectrumScan(
                scan_id=f"{sample_id}_{index}",
                raw_file_id=f"{sample_id}.raw",
                sample_id=sample_id,
                rt=rt,
                ms_level=1,
                mz_array=[mz],
                intensity_array=[intensity],
                tic=intensity,
                base_peak_mz=mz,
                base_peak_intensity=intensity,
            )
        )
    return scans


def composite_peak_scans(
    sample_id: str,
    shift: float,
    first_height: float,
    second_height: float,
) -> list[LCMSSpectrumScan]:
    scans: list[LCMSSpectrumScan] = []
    for index in range(401):
        rt = 4.0 + index * 0.005
        intensity = (
            gaussian(rt, 4.80 + shift, 0.07, first_height)
            + gaussian(rt, 5.15 + shift, 0.09, second_height)
            + 5.0
        )
        scans.append(
            LCMSSpectrumScan(
                scan_id=f"{sample_id}_{index}",
                raw_file_id=f"{sample_id}.raw",
                sample_id=sample_id,
                rt=rt,
                ms_level=1,
                mz_array=[500.0],
                intensity_array=[intensity],
                tic=intensity,
                base_peak_mz=500.0,
                base_peak_intensity=intensity,
            )
        )
    return scans


class LCMSMvpTests(unittest.TestCase):
    def test_peak_first_template_marks_selected_feature_rt_and_uses_one_tic_band(self) -> None:
        self.assertIn("function drawSelectedFeatureRt", PEAK_FIRST_TEMPLATE)
        self.assertIn("Purple dashed line marks selected Feature RT", PEAK_FIRST_TEMPLATE)
        self.assertIn('ctx.fillStyle="rgba(59,130,246,.11)"', PEAK_FIRST_TEMPLATE)
        self.assertIn("data-saved-feature", PEAK_FIRST_TEMPLATE)
        self.assertIn("const MAX_FEATURE_FOLD = 1000", PEAK_FIRST_TEMPLATE)
        self.assertIn("Math.log(MAX_FEATURE_FOLD)", PEAK_FIRST_TEMPLATE)
        self.assertIn("capped |ln(test/reference)|", PEAK_FIRST_TEMPLATE)
        self.assertIn("function featureMapPairAbundance", PEAK_FIRST_TEMPLATE)
        self.assertIn("function featureMapAbundanceScale", PEAK_FIRST_TEMPLATE)
        self.assertIn("const radius=featureMapHalfSide(item,abundanceScale)", PEAK_FIRST_TEMPLATE)
        self.assertIn("square area uses the higher normalized abundance in the selected pair on a robust ln scale", PEAK_FIRST_TEMPLATE)
        self.assertNotIn("abundanceNorm*4+Math.sqrt(rankNorm)*3+foldNorm*5", PEAK_FIRST_TEMPLATE)
        self.assertIn("function logFoldTicks", PEAK_FIRST_TEMPLATE)
        self.assertIn("[1,10,100,1000]", PEAK_FIRST_TEMPLATE)
        self.assertIn("requestedComparison", PEAK_FIRST_TEMPLATE)
        self.assertNotIn("msmsReportLink", PEAK_FIRST_TEMPLATE)
        self.assertNotIn('href="/lcms_msms_report.html"', PEAK_FIRST_TEMPLATE)
        self.assertIn("selectedFeatureGroupId", PEAK_FIRST_TEMPLATE)
        self.assertIn("componentContainsFeature", PEAK_FIRST_TEMPLATE)
        self.assertIn('data-selection-target="${parentSelected?', PEAK_FIRST_TEMPLATE)
        self.assertIn("scrollGlobalSelectionIntoView", PEAK_FIRST_TEMPLATE)
        self.assertIn("container.scrollTop", PEAK_FIRST_TEMPLATE)
        self.assertIn("feature-analysis-grid", PEAK_FIRST_TEMPLATE)
        self.assertIn('id="featureMs2Canvas"', PEAK_FIRST_TEMPLATE)
        self.assertIn("function ensureMsmsData", PEAK_FIRST_TEMPLATE)
        self.assertNotIn('id="globalOpenModificationModule"', PEAK_FIRST_TEMPLATE)
        self.assertNotIn("function renderGlobalOpenModifications", PEAK_FIRST_TEMPLATE)
        self.assertNotIn("已知修饰优先汇总与开放质量筛查", PEAK_FIRST_TEMPLATE)
        self.assertIn('apiUrl("/api/msms")', PEAK_FIRST_TEMPLATE)
        self.assertIn("function selectedMs2Evidence", PEAK_FIRST_TEMPLATE)
        self.assertIn("function componentCandidateEvidence", PEAK_FIRST_TEMPLATE)
        self.assertIn("function componentDisplayLabel", PEAK_FIRST_TEMPLATE)
        self.assertIn("low_evidence_sequence_candidate", PEAK_FIRST_TEMPLATE)
        self.assertIn("displaying the best stored covering scan", PEAK_FIRST_TEMPLATE)
        self.assertIn("component_isotope_fit_score", PEAK_FIRST_TEMPLATE)
        self.assertIn("isotope fit ${nice(isotopeFit,2)}", PEAK_FIRST_TEMPLATE)
        self.assertIn('"monoisotope_corrected"', PEAK_FIRST_TEMPLATE)
        self.assertIn("function drawFeatureMs2", PEAK_FIRST_TEMPLATE)
        self.assertIn('queryParams.get("feature_group_id")', PEAK_FIRST_TEMPLATE)
        self.assertIn("b ions", PEAK_FIRST_TEMPLATE)
        self.assertIn("y ions", PEAK_FIRST_TEMPLATE)
        self.assertIn("feature_group_id", PEAK_FIRST_TEMPLATE)
        self.assertIn('id="globalSequenceTracks"', PEAK_FIRST_TEMPLATE)
        self.assertIn('id="sequenceOverlapDetails"', PEAK_FIRST_TEMPLATE)
        self.assertIn("function residueHasDirectionOverlap", PEAK_FIRST_TEMPLATE)
        self.assertIn("function sequenceEvidenceSegment", PEAK_FIRST_TEMPLATE)
        self.assertIn("选中区域（无方向冲突）", PEAK_FIRST_TEMPLATE)
        self.assertIn("function renderSequenceOverlapDetails", PEAK_FIRST_TEMPLATE)
        self.assertIn("data-overlap-position", PEAK_FIRST_TEMPLATE)
        self.assertIn("diff-overlap-selected", PEAK_FIRST_TEMPLATE)
        self.assertIn("sequence-overlap-modified-site", PEAK_FIRST_TEMPLATE)
        self.assertIn("sequence-overlap-fold", PEAK_FIRST_TEMPLATE)
        self.assertIn("state.sequenceOverlapSelection=null", PEAK_FIRST_TEMPLATE)
        self.assertIn("红蓝重叠：不同 Feature 方向相反", PEAK_FIRST_TEMPLATE)
        self.assertIn("9. 修饰水平与蛋白形式差异定量", PEAK_FIRST_TEMPLATE)
        self.assertIn("10. 差异组分的序列与三级结构定位", PEAK_FIRST_TEMPLATE)
        self.assertIn("function renderModificationQuantitation", PEAK_FIRST_TEMPLATE)
        self.assertIn('id="sequenceHigherLegend"', PEAK_FIRST_TEMPLATE)
        self.assertIn('id="sequenceLowerLegend"', PEAK_FIRST_TEMPLATE)
        self.assertIn("function updateSequenceDirectionLegend", PEAK_FIRST_TEMPLATE)
        self.assertIn("红色：${sampleShort(pair.test)} 高于 ${sampleShort(pair.reference)}", PEAK_FIRST_TEMPLATE)
        self.assertIn("蓝色：${sampleShort(pair.reference)} 高于 ${sampleShort(pair.test)}", PEAK_FIRST_TEMPLATE)
        self.assertIn("function globalSequenceLocations", PEAK_FIRST_TEMPLATE)
        self.assertIn("function preferredStructureChains", PEAK_FIRST_TEMPLATE)
        self.assertIn("PROJECT_STRUCTURE_DEFAULTS", PEAK_FIRST_TEMPLATE)
        self.assertIn('id="saveStructureFile"', PEAK_FIRST_TEMPLATE)
        self.assertIn("/api/task-reference/upload", PEAK_FIRST_TEMPLATE)
        self.assertIn("function loadUploadedStructure(file,save=true)", PEAK_FIRST_TEMPLATE)
        self.assertIn('pdbId:"3V4P"', PEAK_FIRST_TEMPLATE)
        self.assertIn('pdbId:"5I1C"', PEAK_FIRST_TEMPLATE)
        self.assertIn("同源模板，非药物精确结构", PEAK_FIRST_TEMPLATE)
        self.assertNotIn("colored bands are sample-specific integration windows", PEAK_FIRST_TEMPLATE)

    def test_mock_raw_parser_creates_scan_level_arrays(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "SCX-HPLC-MS_Intact_MabThera_Untreated_1.raw"
            path.write_bytes(b"placeholder")
            raw_file, scans = mock_scans_from_raw(path)

        self.assertEqual(raw_file.sample_id, "MabThera_1")
        self.assertGreater(raw_file.scan_count, 100)
        self.assertGreater(scans[0].tic, 0)
        self.assertEqual(len(scans[0].mz_array), len(scans[0].intensity_array))
        self.assertEqual(raw_file.parser_status, "mock_from_vendor_raw")

    def test_candidate_screening_extracts_known_mz(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "SCX-HPLC-MS_Intact_MabThera_Untreated_1.raw"
            path.write_bytes(b"placeholder")
            _, scans = mock_scans_from_raw(path)
        candidates = screen_candidate_mz(scans, 4.0, 8.8, intensity_threshold=1500, min_scan_count=4)
        mzs = [candidate.candidate_mz for candidate in candidates]

        self.assertTrue(any(abs(mz - 548.3124) < 0.02 for mz in mzs))
        self.assertTrue(all(candidate.scan_count >= 4 for candidate in candidates))

    def test_xic_peak_detection_integrates_area(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "SCX-HPLC-MS_Intact_MabThera_Untreated_1.raw"
            path.write_bytes(b"placeholder")
            _, scans = mock_scans_from_raw(path)
        xic = extract_xic(scans, 548.3124, 20)
        features = detect_xic_peaks(xic, min_peak_height=2500, min_peak_area=200, min_snr=3)

        self.assertGreaterEqual(len(features), 1)
        self.assertAlmostEqual(features[0].rt_apex, 4.8, delta=0.15)
        self.assertGreater(features[0].area, 0)

    def test_feature_matching_and_difference_classification(self) -> None:
        all_features = []
        sample_ids = []
        with tempfile.TemporaryDirectory() as tmp:
            for name in [
                "SCX-HPLC-MS_Intact_MabThera_Untreated_1.raw",
                "SCX-HPLC-MS_Intact_Reditux_Untreated_1.raw",
            ]:
                path = Path(tmp) / name
                path.write_bytes(b"placeholder")
                raw_file, scans = mock_scans_from_raw(path)
                sample_ids.append(raw_file.sample_id)
                xic = extract_xic(scans, 732.4411, 20)
                all_features.extend(detect_xic_peaks(xic, min_peak_height=2500, min_peak_area=200, min_snr=3)[:1])

        matched = match_features(all_features, sample_ids, mz_tolerance_ppm=20, rt_tolerance_min=0.2)
        groups = build_feature_groups(
            matched,
            sample_ids,
            reference_samples={"MabThera_1"},
            target_samples={"Reditux_1"},
            fold_change_threshold=1.5,
        )

        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0].difference_type, "increased")
        self.assertGreater(groups[0].area_by_sample["Reditux_1"], groups[0].area_by_sample["MabThera_1"])

    def test_main_peak_alignment_matches_reference_peak_window(self) -> None:
        scans_by_sample = {
            "reference": synthetic_chrom_scans("reference", main_rt=5.0, impurity_rt=8.0, impurity_height=350.0),
            "sample": synthetic_chrom_scans("sample", main_rt=5.4, impurity_rt=8.4, impurity_height=1500.0),
        }

        alignment = main_peak_alignment(
            scans_by_sample,
            reference_sample="reference",
            rt_start=2.0,
            rt_end=10.0,
            signal="bpc",
            match_window_min=1.0,
        )

        self.assertAlmostEqual(alignment["apex_by_sample"]["reference"], 5.0, delta=0.11)
        self.assertAlmostEqual(alignment["apex_by_sample"]["sample"], 5.4, delta=0.11)
        self.assertAlmostEqual(alignment["rt_shift_by_sample"]["sample"], -0.4, delta=0.11)
        self.assertGreaterEqual(alignment["alignment_landmark_count_by_sample"]["sample"], 2)
        self.assertEqual(alignment["main_peak_by_sample"]["sample"]["match_method"], "multi_peak_landmark_median")

    def test_local_tic_alignment_uses_composite_profile_and_keeps_reference_fixed(self) -> None:
        scans_by_sample = {
            "reference": composite_peak_scans("reference", 0.0, 1100.0, 900.0),
            # The second shoulder is now the highest point. Highest-point
            # alignment would pair it with the reference's first subpeak.
            "sample": composite_peak_scans("sample", 0.12, 900.0, 1050.0),
        }
        tic_peak = {
            "tic_peak_id": "TICP_TEST",
            "rt_start": 4.55,
            "rt_apex": 4.80,
            "rt_end": 5.45,
            "width": 0.90,
        }
        local = local_align_peak_by_apex(
            scans_by_sample,
            tic_peak,
            shifts={"reference": 0.0, "sample": 0.0},
            params=PeakFirstParams(max_local_shift_cap_min=0.3, max_local_shift_fraction=0.3),
            reference_sample="reference",
        )

        self.assertEqual(local["alignment_method"], "reference_profile_cross_correlation")
        self.assertEqual(local["local_peak_shift_by_sample"]["reference"], 0.0)
        self.assertAlmostEqual(local["local_peak_shift_by_sample"]["sample"], -0.12, delta=0.015)
        self.assertGreater(
            local["profile_score_by_sample"]["sample"],
            local["profile_score_before_by_sample"]["sample"] + 0.05,
        )
        self.assertFalse(local["local_shift_rejected_by_sample"]["sample"])

    def test_local_tic_alignment_rejects_flat_window_with_boundary_peak(self) -> None:
        scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
        for sample_id, center in (("reference", 5.42), ("sample", 5.34)):
            scans: list[LCMSSpectrumScan] = []
            for index in range(401):
                rt = 4.0 + index * 0.005
                intensity = gaussian(rt, center, 0.07, 1000.0) + 5.0
                scans.append(
                    LCMSSpectrumScan(
                        scan_id=f"{sample_id}_{index}",
                        raw_file_id=f"{sample_id}.raw",
                        sample_id=sample_id,
                        rt=rt,
                        ms_level=1,
                        mz_array=[500.0],
                        intensity_array=[intensity],
                        tic=intensity,
                        base_peak_mz=500.0,
                        base_peak_intensity=intensity,
                    )
                )
            scans_by_sample[sample_id] = scans
        local = local_align_peak_by_apex(
            scans_by_sample,
            {
                "tic_peak_id": "TICP_FLAT",
                "rt_start": 4.60,
                "rt_apex": 5.10,
                "rt_end": 5.20,
                "width": 0.60,
            },
            shifts={"reference": 0.0, "sample": 0.0},
            params=PeakFirstParams(max_local_shift_cap_min=0.25, max_local_shift_fraction=0.4),
            reference_sample="reference",
        )

        self.assertFalse(local["local_alignment_quality_gate_passed"])
        self.assertEqual(local["local_peak_shift_by_sample"]["sample"], 0.0)
        self.assertTrue(local["local_shift_rejected_by_sample"]["sample"])
        self.assertIn(
            "apex_at_window_edge",
            local["peak_profile_quality_by_sample"]["reference"]["reason"],
        )

    def test_heatmap_normalization_reduces_scale_bias(self) -> None:
        scans_by_sample = {
            "reference": scaled_single_mz_scans("reference", 1.0),
            "sample": scaled_single_mz_scans("sample", 10.0),
        }

        heatmaps = build_similarity_heatmaps(
            scans_by_sample,
            reference_sample="reference",
            shifts={"reference": 0.0, "sample": 0.0},
            rt_min=0.0,
            rt_max=2.0,
            mz_min=499.0,
            mz_max=501.0,
            rt_bin_count=20,
            mz_bin_count=1,
            normalization_methods=["none", "max"],
        )

        by_method = {heatmap["normalization_method"]: heatmap for heatmap in heatmaps}
        raw_score = by_method["none"]["scores"][0][10]
        normalized_score = by_method["max"]["scores"][0][10]

        self.assertLess(raw_score, 0.6)
        self.assertGreater(normalized_score, 0.95)
        first_region = by_method["none"]["low_similarity_points"][0]
        self.assertIn("fold_change", first_region)
        self.assertIn("difference_type", first_region)

    def test_heatmap_similarity_tolerates_local_rt_mz_drift(self) -> None:
        scans_by_sample = {
            "reference": shifted_single_feature_scans("reference", 1.0, 500.0),
            "sample": shifted_single_feature_scans("sample", 1.08, 500.6),
        }

        heatmaps = build_similarity_heatmaps(
            scans_by_sample,
            reference_sample="reference",
            shifts={"reference": 0.0, "sample": 0.0},
            rt_min=0.0,
            rt_max=2.0,
            mz_min=499.0,
            mz_max=501.0,
            rt_bin_count=20,
            mz_bin_count=4,
            normalization_methods=["max"],
        )

        heatmap = heatmaps[0]
        self.assertEqual(heatmap["similarity_method"], "local_neighborhood_max")
        self.assertGreater(heatmap["scores"][2][10], 0.9)
        self.assertLess(len(heatmap["low_similarity_points"]), 5)

    def test_heatmap_empty_bins_are_null_not_colored(self) -> None:
        scans_by_sample = {
            "reference": shifted_single_feature_scans("reference", 1.0, 500.0),
            "sample": shifted_single_feature_scans("sample", 1.0, 500.0),
        }

        heatmaps = build_similarity_heatmaps(
            scans_by_sample,
            reference_sample="reference",
            shifts={"reference": 0.0, "sample": 0.0},
            rt_min=0.0,
            rt_max=4.0,
            mz_min=499.0,
            mz_max=503.0,
            rt_bin_count=40,
            mz_bin_count=4,
            normalization_methods=["max"],
        )

        self.assertIsNone(heatmaps[0]["scores"][3][35])

    def test_cohort_cv_heatmap_is_generated_for_multiple_samples(self) -> None:
        scans_by_sample = {
            "reference": scaled_single_mz_scans("reference", 1.0),
            "biosimilar_a": scaled_single_mz_scans("biosimilar_a", 1.0),
            "biosimilar_b": scaled_single_mz_scans("biosimilar_b", 5.0),
        }

        heatmaps = build_similarity_heatmaps(
            scans_by_sample,
            reference_sample="reference",
            shifts={"reference": 0.0, "biosimilar_a": 0.0, "biosimilar_b": 0.0},
            rt_min=0.0,
            rt_max=2.0,
            mz_min=499.0,
            mz_max=501.0,
            rt_bin_count=20,
            mz_bin_count=1,
            normalization_methods=["none", "max"],
        )

        cohorts = {
            heatmap["normalization_method"]: heatmap
            for heatmap in heatmaps
            if heatmap["comparison_type"] == "cohort_cv"
        }

        self.assertEqual(set(cohorts), {"none", "max"})
        self.assertEqual(cohorts["none"]["comparison_key"], "cohort:all_samples")
        self.assertLess(cohorts["none"]["scores"][0][10], 0.75)
        self.assertGreater(cohorts["max"]["scores"][0][10], 0.95)
        first_region = cohorts["none"]["low_similarity_points"][0]
        self.assertIn("sample_intensities", first_region)
        self.assertEqual(set(first_region["sample_intensities"]), set(scans_by_sample))

    def test_difference_ranking_downweights_low_abundance_regions(self) -> None:
        regions = [
            {
                "rt_index": 1,
                "mz_index": 1,
                "difference_score": 0.5,
                "max_intensity": 100.0,
                "group_mean_intensity": 0.05,
            },
            {
                "rt_index": 10,
                "mz_index": 10,
                "difference_score": 0.45,
                "max_intensity": 10000.0,
                "group_mean_intensity": 5.0,
            },
        ]

        sort_difference_regions(regions)

        self.assertEqual(regions[0]["rt_index"], 10)
        self.assertLess(regions[1]["weighted_difference_score"], regions[1]["difference_score"])
        self.assertGreater(regions[0]["weighted_difference_score"], regions[1]["weighted_difference_score"])

    def test_peak_first_payload_scores_tic_peaks_and_changed_mz(self) -> None:
        raw_files = []
        scans_by_sample = {
            "reference": shifted_single_feature_scans("reference", 1.0, 500.0),
            "sample": shifted_single_feature_scans("sample", 1.0, 500.0),
            "changed": shifted_single_feature_scans("changed", 1.0, 501.0),
        }
        params = PeakFirstParams(
            min_snr=2,
            min_area_ratio=0.0,
            min_prominence_factor=1.0,
            min_width=0.01,
            max_width=2.0,
            top_n_peaks=5,
            top_n_mz=20,
            top_n_changed_mz=5,
            mz_tolerance_da=0.2,
        )

        payload = prepare_peak_first_payload(
            raw_files,
            scans_by_sample,
            project_id="test_peak_first",
            reference_sample="reference",
            params=params,
        )

        self.assertEqual(payload["module"], "LCMSPeakFirstCompare")
        self.assertGreaterEqual(len(payload["peak_results"]), 1)
        first_peak = payload["peak_results"][0]
        self.assertIn("spectrum_score", first_peak)
        self.assertIn("peak_consistency_score", first_peak)
        self.assertTrue(first_peak["top_changed_mz"])
        self.assertIn("difference_type", first_peak["top_changed_mz"][0])
        self.assertEqual(payload["feature_area_normalization"]["method"], "total_tic_area_to_median")

    def test_peak_first_abundance_types_use_true_fold_change(self) -> None:
        params = PeakFirstParams(
            feature_presence_relative_area_fraction=0.02,
            feature_common_fold_change_threshold=2.0,
            feature_strong_fold_change_threshold=4.0,
        )
        detected = {"reference": True, "sample": True}

        common = abundance_profile_metrics({"reference": 100.0, "sample": 199.0}, detected, params)
        moderate = abundance_profile_metrics({"reference": 100.0, "sample": 200.0}, detected, params)
        changed = abundance_profile_metrics({"reference": 100.0, "sample": 400.0}, detected, params)
        absent = abundance_profile_metrics({"reference": 100.0, "sample": 2.0}, detected, params)

        self.assertEqual(common["difference_type"], "common_feature")
        self.assertAlmostEqual(common["max_fold_change"], 1.99)
        self.assertEqual(moderate["difference_type"], "moderate_difference")
        self.assertEqual(changed["difference_type"], "area_changed")
        self.assertEqual(absent["difference_type"], "presence_absence")
        self.assertFalse(absent["presence_by_sample"]["sample"])

        extreme = abundance_profile_metrics(
            {"reference": 1.0, "sample": 1_000_000.0},
            detected,
            params,
        )
        self.assertEqual(extreme["max_fold_change"], 1000.0)
        self.assertEqual(extreme["raw_max_fold_change"], 1_000_000.0)

        partial = abundance_profile_metrics(
            {"reference": 100.0, "sample": 156.0},
            {"reference": True, "sample": False},
            params,
        )
        self.assertEqual(partial["difference_type"], "common_feature")
        self.assertEqual(partial["presence_count"], 2)
        self.assertEqual(partial["detection_count"], 1)
        self.assertEqual(partial["quantitation_confidence"], "partial_detection")
        self.assertEqual(partial["ranking_confidence_weight"], 0.35)
        self.assertAlmostEqual(partial["similarity_score"], 100.0 / 156.0)

        undetected = abundance_profile_metrics(
            {"reference": 100.0, "sample": 156.0},
            {"reference": False, "sample": False},
            params,
        )
        self.assertEqual(undetected["difference_type"], "low_confidence")
        self.assertEqual(undetected["ranking_confidence_weight"], 0.0)

    def test_peak_first_xic_detection_uses_configured_snr(self) -> None:
        xic = [(index * 0.1, value, 500.0) for index, value in enumerate([1.0, 1.0, 4.0, 1.0, 1.0])]
        strict = detect_xic_lcms_feature(
            xic,
            target_mz=500.0,
            sample_id="sample",
            raw_file_id="sample.mzML",
            tic_peak_id="TICP_0001",
            params=PeakFirstParams(min_snr=5.0),
        )
        permissive = detect_xic_lcms_feature(
            xic,
            target_mz=500.0,
            sample_id="sample",
            raw_file_id="sample.mzML",
            tic_peak_id="TICP_0001",
            params=PeakFirstParams(min_snr=2.0),
        )

        self.assertEqual(strict["match_status"], "gap_filled")
        self.assertEqual(permissive["match_status"], "matched")
        self.assertEqual(permissive["observed_scan_count"], 5)
        self.assertEqual(permissive["max_consecutive_observed_scans"], 5)

    def test_dynamic_background_candidates_are_thresholded_not_top15_limited(self) -> None:
        bins = [300.0 + index for index in range(40)]
        strong = [10_000.0 for _ in range(20)]
        weak = [100.0 for _ in range(20)]
        raw_matrix = {
            "reference": strong + weak,
            "sample": [value * (2.0 if index % 2 else 0.5) for index, value in enumerate(strong + weak)],
        }
        fixed = find_top_changed_mz(
            bins,
            raw_matrix,
            PeakFirstParams(top_n_changed_mz=15),
        )
        dynamic = find_top_changed_mz(
            bins,
            raw_matrix,
            PeakFirstParams(
                top_n_changed_mz=15,
                dynamic_background_candidates=True,
                candidate_spectral_noise_multiplier=3.0,
                candidate_min_local_tic_ppm=0.0,
                candidate_max_per_tic=200,
            ),
        )

        self.assertEqual(len(fixed), 15)
        self.assertEqual(len(dynamic), 20)
        self.assertTrue(all(float(item["candidate_spectral_noise_ratio"]) >= 3.0 for item in dynamic))

    def test_centroid_consensus_keeps_z3_isotopes_as_real_peaks(self) -> None:
        params = PeakFirstParams(mz_tolerance_ppm=10.0, mz_tolerance_da=0.16)
        spectra = [
            [(627.6874, 100.0), (628.0218, 80.0), (628.3558, 60.0)],
            [(627.6878, 90.0), (628.0221, 70.0), (628.3561, 50.0)],
        ]

        bins = build_consensus_mz_bins(spectra, params)

        self.assertEqual(len(bins), 3)
        self.assertAlmostEqual(bins[0], 627.68759, places=4)
        self.assertAlmostEqual(bins[1] - bins[0], 1.00335483507 / 3, places=3)

    def test_single_centroid_xic_does_not_average_neighbor_inside_da_window(self) -> None:
        scan = LCMSSpectrumScan(
            scan_id="scan=1",
            raw_file_id="sample.mzML",
            sample_id="sample",
            rt=10.0,
            ms_level=1,
            mz_array=[500.0000, 500.1200],
            intensity_array=[100.0, 900.0],
            tic=1000.0,
            base_peak_mz=500.12,
            base_peak_intensity=900.0,
        )
        points = extract_xic_points_for_peak(
            [scan],
            {"rt_start": 9.9, "rt_end": 10.1},
            target_mz=500.0,
            rt_shift=0.0,
            local_shift=0.0,
            params=PeakFirstParams(mz_tolerance_ppm=10.0, mz_tolerance_da=0.16),
        )

        self.assertEqual(len(points), 1)
        self.assertAlmostEqual(points[0][1], 100.0)
        self.assertAlmostEqual(points[0][2], 500.0)

    def test_single_centroid_xic_does_not_substitute_neighbor_when_target_is_missing(self) -> None:
        scan = LCMSSpectrumScan(
            scan_id="scan=1",
            raw_file_id="sample.mzML",
            sample_id="sample",
            rt=10.0,
            ms_level=1,
            mz_array=[500.1200],
            intensity_array=[900.0],
            tic=900.0,
            base_peak_mz=500.12,
            base_peak_intensity=900.0,
        )
        points = extract_xic_points_for_peak(
            [scan],
            {"rt_start": 9.9, "rt_end": 10.1},
            target_mz=500.0,
            rt_shift=0.0,
            local_shift=0.0,
            params=PeakFirstParams(mz_tolerance_ppm=10.0, mz_tolerance_da=0.16),
        )

        self.assertEqual(len(points), 1)
        self.assertEqual(points[0][1], 0.0)
        self.assertEqual(points[0][2], 500.0)

    def test_browser_xic_does_not_plot_neighbor_as_target_signal(self) -> None:
        params = {"mz_tolerance_ppm": 10.0, "mz_tolerance_da": 0.16, "mz_tolerance_mode": "da"}

        exact = single_centroid_intensity([500.0, 500.12], [100.0, 900.0], 500.0, 0.16, params)
        missing = single_centroid_intensity([500.12], [900.0], 500.0, 0.16, params)

        self.assertEqual(exact, 100.0)
        self.assertEqual(missing, 0.0)

    def test_global_feature_merge_keeps_type_and_fold_from_same_row(self) -> None:
        common = {
            "feature_group_id": "FG_common",
            "parent_tic_peak_id": "TICP_1",
            "representative_mz": 500.0,
            "representative_rt": 10.0,
            "similarity_score": 0.8,
            "difference_score": 0.2,
            "abundance_weighted_score": 0.3,
            "abundance_score": 10.0,
            "max_area": 1000.0,
            "max_fold_change": 1.25,
            "difference_type": "common_feature",
        }
        moderate = {
            **common,
            "feature_group_id": "FG_moderate",
            "parent_tic_peak_id": "TICP_2",
            "representative_mz": 500.004,
            "representative_rt": 10.2,
            "abundance_weighted_score": 0.2,
            "max_fold_change": 3.0,
            "difference_type": "moderate_difference",
        }

        rows = global_feature_groups_from_peaks(
            [{"feature_groups": [common]}, {"feature_groups": [moderate]}],
            PeakFirstParams(mz_tolerance_ppm=10.0, mz_tolerance_da=0.16),
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["difference_type"], "common_feature")
        self.assertEqual(rows[0]["max_fold_change"], 1.25)

    def test_peak_first_splits_adjacent_tic_peaks_by_valley(self) -> None:
        curve = []
        for index in range(401):
            rt = index * 0.01
            intensity = (
                gaussian(rt, 1.00, 0.04, 1000)
                + gaussian(rt, 1.22, 0.04, 850)
                + gaussian(rt, 1.55, 0.04, 900)
                + 10
            )
            curve.append((rt, intensity))
        params = PeakFirstParams(
            min_snr=2,
            min_prominence_factor=1,
            min_area_ratio=0,
            min_width=0.02,
            max_width=0.4,
            min_apex_distance_min=0.08,
            peak_boundary_fraction=0.02,
            min_valley_depth_fraction=0.10,
        )

        peaks = filter_tic_peaks(detect_tic_peaks(curve, params), params)
        confirmed = [peak for peak in peaks if peak["status"] == "confirmed_peak"]

        self.assertGreaterEqual(len(confirmed), 3)
        self.assertTrue(all(float(peak["width"]) < 0.4 for peak in confirmed[:3]))
        windows = {(round(float(peak["rt_start"]), 3), round(float(peak["rt_end"]), 3)) for peak in confirmed}
        self.assertEqual(len(windows), len(confirmed))

    def test_spectrum_payload_keeps_scan_trace_fields(self) -> None:
        scans_by_sample = {"reference": scaled_single_mz_scans("reference", 1.0)}
        payload = spectrum_payload(
            scans_by_sample,
            shifts={"reference": 0.05},
            mz_min=499.0,
            mz_max=501.0,
            max_peaks_per_scan=10,
            min_intensity=0.0,
        )
        first_scan = payload["reference"][0]

        self.assertEqual(first_scan["scan_id"], "reference_0")
        self.assertEqual(first_scan["raw_file_id"], "reference.raw")
        self.assertAlmostEqual(first_scan["aligned_rt"], 0.05, delta=1e-9)

    def test_auto_feature_regions_build_matrix_rows(self) -> None:
        scans_by_sample = {
            "reference": scaled_single_mz_scans("reference", 1.0),
            "biosimilar_a": scaled_single_mz_scans("biosimilar_a", 1.0),
            "biosimilar_b": scaled_single_mz_scans("biosimilar_b", 5.0),
        }
        heatmaps = build_similarity_heatmaps(
            scans_by_sample,
            reference_sample="reference",
            shifts={"reference": 0.0, "biosimilar_a": 0.0, "biosimilar_b": 0.0},
            rt_min=0.0,
            rt_max=2.0,
            mz_min=499.0,
            mz_max=501.0,
            rt_bin_count=20,
            mz_bin_count=1,
            normalization_methods=["none", "max"],
        )
        payload = {
            "sample_ids": list(scans_by_sample),
            "default_heatmap_normalization": "none",
            "heatmaps": heatmaps,
        }

        regions = auto_feature_regions_from_payload(payload, top_n_per_heatmap=2)
        header = feature_matrix_header(list(scans_by_sample))
        rows = feature_matrix_rows_from_regions(regions, list(scans_by_sample))

        self.assertGreaterEqual(len(regions), 2)
        self.assertTrue(all(region["normalization_method"] == "none" for region in regions))
        self.assertIn("biosimilar_b raw", header)
        self.assertEqual(len(rows[0]), len(header))
        self.assertIn("sample_intensities", regions[0])
        self.assertEqual(set(regions[0]["sample_intensities"]), set(scans_by_sample))

    def test_saved_lcms_features_survive_workbench_regeneration(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db_path = Path(tmp) / "lcms_workbench.sqlite"
            feature = {
                "region_id": "LCMSF_persist_1",
                "aligned_rt": 1.25,
                "mz": 500.0,
                "sample_intensities": {"reference": {"raw": 10.0, "normalized": 1.0}},
            }

            self.assertEqual(replace_saved_features(db_path, [feature]), 1)
            payload = {
                "sample_ids": ["reference"],
                "spectra": {"reference": []},
                "heatmaps": [],
                "reference_sample": "reference",
                "default_heatmap_normalization": "max",
            }
            write_workbench_sqlite(db_path, payload)
            saved = read_saved_features(db_path)

        self.assertEqual(len(saved), 1)
        self.assertEqual(saved[0]["region_id"], "LCMSF_persist_1")
        self.assertEqual(saved[0]["sample_intensities"]["reference"]["raw"], 10.0)

    def test_saved_feature_matrix_includes_trace_fields(self) -> None:
        features = [
            {
                "region_id": "LCMSF_trace_1",
                "comparison_label": "reference vs biosimilar",
                "comparison_type": "pair",
                "normalization_method": "max",
                "aligned_rt": 1.25,
                "mz": 500.0,
                "sample_presence": ["reference"],
                "sample_intensities": {
                    "reference": {"raw": 10.0, "normalized": 1.0},
                    "biosimilar": {"raw": 0.0, "normalized": 0.0},
                },
                "source_scans": {
                    "reference": {
                        "scan_id": "reference_12",
                        "raw_rt": 1.2,
                        "aligned_rt": 1.25,
                    }
                },
            }
        ]

        matrix = build_feature_matrix(features, ["reference", "biosimilar"])
        csv_text = feature_matrix_csv(matrix)

        self.assertEqual(matrix["feature_count"], 1)
        self.assertIn("reference scan_id", matrix["columns"])
        self.assertEqual(matrix["rows"][0]["reference scan_id"], "reference_12")
        self.assertEqual(matrix["rows"][0]["biosimilar present"], "missing")
        self.assertIn("LCMSF_trace_1", csv_text)
        self.assertIn("reference_12", csv_text)

    def test_heatmap_window_slices_scores_and_intensity_grids(self) -> None:
        heatmap = {
            "comparison_type": "cohort_cv",
            "comparison_key": "cohort:all_samples",
            "normalization_method": "max",
            "rt_min": 0.0,
            "rt_max": 4.0,
            "mz_min": 100.0,
            "mz_max": 104.0,
            "rt_bin_count": 4,
            "mz_bin_count": 4,
            "scores": [
                [0.1, 0.2, 0.3, 0.4],
                [0.5, 0.6, 0.7, 0.8],
                [0.9, 1.0, 0.9, 0.8],
                [0.7, 0.6, 0.5, 0.4],
            ],
            "sample_intensity_grids": {
                "reference": [
                    [1, 2, 3, 4],
                    [5, 6, 7, 8],
                    [9, 10, 11, 12],
                    [13, 14, 15, 16],
                ]
            },
            "low_similarity_points": [
                {"rt": 1.5, "mz": 101.5, "score": 0.6},
                {"rt": 3.5, "mz": 103.5, "score": 0.4},
            ],
        }

        window = heatmap_window(heatmap, rt_min=1.0, rt_max=2.0, mz_min=101.0, mz_max=102.0)

        self.assertEqual(window["rt_bin_count"], 2)
        self.assertEqual(window["mz_bin_count"], 2)
        self.assertEqual(window["source_rt_bin_offset"], 1)
        self.assertEqual(window["source_mz_bin_offset"], 1)
        self.assertEqual(window["scores"], [[0.6, 0.7], [1.0, 0.9]])
        self.assertEqual(window["sample_intensity_grids"]["reference"], [[6, 7], [10, 11]])
        self.assertEqual(len(window["low_similarity_points"]), 1)
        self.assertEqual(window["low_similarity_points"][0]["rt"], 1.5)


if __name__ == "__main__":
    unittest.main()
