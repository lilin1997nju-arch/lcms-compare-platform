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
from core.lcms_peak_first import PeakFirstParams, detect_tic_peaks, filter_tic_peaks, prepare_peak_first_payload  # noqa: E402
from core.lcms_workbench import build_similarity_heatmaps, main_peak_alignment, sort_difference_regions, spectrum_payload  # noqa: E402
from core.lcms_xic import extract_xic, screen_candidate_mz  # noqa: E402
from run_xcalibur_workbench import (  # noqa: E402
    auto_feature_regions_from_payload,
    feature_matrix_header,
    feature_matrix_rows_from_regions,
    write_workbench_sqlite,
)
from serve_xcalibur_workbench import (  # noqa: E402
    build_feature_matrix,
    feature_matrix_csv,
    heatmap_window,
    read_saved_features,
    replace_saved_features,
)


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


class LCMSMvpTests(unittest.TestCase):
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
