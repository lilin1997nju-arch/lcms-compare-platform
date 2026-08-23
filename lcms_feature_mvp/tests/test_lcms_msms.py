from __future__ import annotations

import base64
import math
import struct
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from core.lcms_models import LCMSSpectrumScan  # noqa: E402
from core.lcms_msms import (  # noqa: E402
    ISOTOPE_MASS_DIFF,
    PROTON,
    TARGETED_N_GLYCANS,
    annotate_feature_groups,
    build_modification_level_quantitation,
    build_proteolytic_modification_families,
    build_peptide_form_comparisons,
    build_feature_ms2_evidence,
    build_global_known_modification_summary,
    build_ms1_component_groups,
    build_modified_peptide_findings,
    differential_annotations,
    fit_theoretical_isotope_envelope,
    format_modification_alternatives,
    generate_candidates,
    generate_sequence_inference_candidates,
    read_fasta,
    search_component_consensus_scans,
    search_component_sequence_tag_scans,
    search_feature_glycopeptide_scans,
    search_feature_open_mass_scans,
    search_global_open_modification_scans,
    search_scans,
    search_selected_feature_consensus_scans,
    theoretical_averagine_isotope_distribution,
    theoretical_fragments,
)
from core.lcms_parser import read_mzml  # noqa: E402
from run_msms_compare import component_search_excluded_feature_ids  # noqa: E402


def encoded(values: list[float]) -> str:
    return base64.b64encode(struct.pack(f"<{len(values)}d", *values)).decode("ascii")


class LCMSMSMVPTests(unittest.TestCase):
    def test_ambiguous_modification_sites_use_compact_residue_aware_labels(self) -> None:
        label = format_modification_alternatives(
            [
                "PyroGlu-Q@1; Glycation@12",
                "PyroGlu-Q@1; Glycation@13",
            ],
            "QVQLVQSGAEVKKPGASVK",
        )
        deamidation = format_modification_alternatives(
            ["Deamidation@12", "Deamidation@16", "Deamidation@4"],
            "SYGNTYLSWYLQKPGQSPQLLIYGISNR",
        )
        competing = format_modification_alternatives(
            ["Glycation@5", "N-terminal glycation@1"],
            "ADYEKHKVYACEVTHQGLSSPVTK",
        )

        self.assertEqual(label, "PyroGlu-Q@1；Glycation@K12/K13（位点未区分）")
        self.assertEqual(deamidation, "Deamidation@N4/Q12/Q16（位点未区分）")
        self.assertIn("方案A：Glycation@K5", competing)
        self.assertIn("方案B：N-terminal glycation@1", competing)
        self.assertTrue(competing.endswith("（方案未区分）"))

    def test_averagine_envelope_fit_recovers_missing_monoisotopic_peak(self) -> None:
        neutral_mass = 1880.0
        theoretical = theoretical_averagine_isotope_distribution(neutral_mass, 10)
        fit = fit_theoretical_isotope_envelope(
            [value * 1000.0 for value in theoretical[1:5]],
            neutral_mass + ISOTOPE_MASS_DIFF,
        )

        self.assertAlmostEqual(sum(theoretical), 1.0, places=8)
        self.assertEqual(fit["monoisotopic_offset"], 1)
        self.assertGreater(fit["cosine_score"], 0.99)
        self.assertGreater(fit["fit_score"], 0.95)

    def test_component_search_can_rescue_d_level_feature_links(self) -> None:
        excluded = component_search_excluded_feature_ids([
            {
                "feature_or_candidate_id": "FG_B",
                "confidence": "B_high_confidence_inferred",
            },
            {
                "feature_or_candidate_id": "FG_C",
                "confidence": "C_tentative",
            },
            {
                "feature_or_candidate_id": "FG_D",
                "confidence": "D_low_evidence",
            },
        ])

        self.assertEqual(excluded, {"FG_B", "FG_C"})

    def test_component_mass_offset_search_uses_joint_sequence_tags(self) -> None:
        candidates = generate_sequence_inference_candidates(
            {"HC": "PEPMIDEK"},
            max_missed_cleavages=0,
            max_terminal_trim=0,
        )
        target = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and candidate.sequence == "PEPMIDEK"
        )
        mass_delta = 15.99491462
        component_mass = target.neutral_mass + mass_delta
        site_index = 3
        theoretical = {
            label: mz
            for label, mz, _, __ in theoretical_fragments(target)
            if "^" not in label
        }

        def fragment_mz(label: str) -> float:
            series = label[0]
            ordinal = int(label[1:])
            contains_site = (
                (series == "b" and ordinal > site_index)
                or (
                    series == "y"
                    and ordinal >= len(target.sequence) - site_index
                )
            )
            return theoretical[label] + (
                mass_delta if contains_site else 0.0
            )

        def scan(
            scan_id: str,
            precursor_mz: float,
            charge: int,
            labels: list[str],
        ) -> LCMSSpectrumScan:
            mz_values = [fragment_mz(label) for label in labels]
            return LCMSSpectrumScan(
                scan_id=scan_id,
                raw_file_id="raw",
                sample_id="sample",
                rt=10.0,
                ms_level=2,
                mz_array=mz_values,
                intensity_array=[
                    1000.0 + index * 100.0
                    for index in range(len(mz_values))
                ],
                tic=10000.0,
                base_peak_mz=mz_values[-1],
                base_peak_intensity=2000.0,
                precursor_mz=precursor_mz,
                precursor_charge=charge,
                precursor_intensity=10000.0,
                activation_method="HCD",
                collision_energy=25.0,
            )

        z2_mz = component_mass / 2 + PROTON
        z3_mz = component_mass / 3 + PROTON
        scans = [
            scan(
                "scan=z2",
                z2_mz,
                2,
                ["b1", "b2", "b3", "y1", "y2", "y3"],
            ),
            scan(
                "scan=z3",
                z3_mz,
                3,
                ["b4", "b5", "b6", "y4", "y5", "y6"],
            ),
        ]
        component = {
            "component_group_id": "MS1INFER_TAG",
            "inferred_component": True,
            "component_confidence": "MS1_high",
            "component_neutral_mass": component_mass,
            "component_representative_charge": 2,
            "members": [
                {
                    "feature_group_id": "FG_Z2",
                    "representative_mz": z2_mz,
                    "true_peak_mz": z2_mz,
                    "representative_rt": 10.0,
                    "component_charge": 2,
                    "component_isotope_offset": 0,
                },
                {
                    "feature_group_id": "FG_Z3",
                    "representative_mz": z3_mz,
                    "true_peak_mz": z3_mz,
                    "representative_rt": 10.0,
                    "component_charge": 3,
                    "component_isotope_offset": 0,
                },
            ],
        }

        accepted, region = search_component_sequence_tag_scans(
            {"alignment": {"rt_shift_by_sample": {"sample": 0.0}}},
            scans,
            candidates,
            [component],
        )

        self.assertEqual(region, [])
        self.assertEqual(len(accepted), 1)
        self.assertEqual(accepted[0]["sequence"], "PEPMIDEK")
        self.assertEqual(
            accepted[0]["sequence_inference_level"],
            "backbone_sequence_supported",
        )
        self.assertAlmostEqual(accepted[0]["mass_delta"], mass_delta, places=4)
        self.assertIn(
            "compatible with Oxidation",
            accepted[0]["modification_text"],
        )
        self.assertGreaterEqual(accepted[0]["sequence_tag_length"], 3)
        self.assertGreaterEqual(
            accepted[0]["complementary_ion_pair_count"],
            1,
        )

    def test_short_backbone_with_large_unknown_delta_is_not_accepted(self) -> None:
        candidates = generate_sequence_inference_candidates(
            {"HC": "PEPTIK"},
            max_missed_cleavages=0,
            max_terminal_trim=0,
        )
        target = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and candidate.sequence == "PEPTIK"
        )
        mass_delta = 123.456
        component_mass = target.neutral_mass + mass_delta
        fragment_map = {
            label: mz
            for label, mz, _, __ in theoretical_fragments(target)
            if "^" not in label
        }

        def shifted_mz(label: str) -> float:
            ordinal = int(label[1:])
            contains_site = (
                (label.startswith("b") and ordinal > 2)
                or (label.startswith("y") and ordinal >= 4)
            )
            return fragment_map[label] + (
                mass_delta if contains_site else 0.0
            )

        scans = []
        for charge, labels in (
            (2, ["b1", "b2", "b3", "y1", "y2"]),
            (3, ["b4", "b5", "y3", "y4", "y5"]),
        ):
            mz_values = [shifted_mz(label) for label in labels]
            scans.append(LCMSSpectrumScan(
                scan_id=f"short=z{charge}",
                raw_file_id="raw",
                sample_id="sample",
                rt=8.0,
                ms_level=2,
                mz_array=mz_values,
                intensity_array=[1800.0] * len(mz_values),
                tic=1800.0 * len(mz_values),
                base_peak_mz=mz_values[0],
                base_peak_intensity=1800.0,
                precursor_mz=component_mass / charge + PROTON,
                precursor_charge=charge,
                precursor_intensity=9000.0,
                activation_method="HCD",
                collision_energy=25.0,
            ))
        component = {
            "component_group_id": "MS1INFER_SHORT_DELTA",
            "inferred_component": True,
            "component_neutral_mass": component_mass,
            "component_representative_charge": 2,
            "members": [
                {
                    "feature_group_id": f"SHORT_Z{charge}",
                    "true_peak_mz": component_mass / charge + PROTON,
                    "representative_mz": component_mass / charge + PROTON,
                    "representative_rt": 8.0,
                    "component_charge": charge,
                    "component_isotope_offset": 0,
                }
                for charge in (2, 3)
            ],
        }

        accepted, region = search_component_sequence_tag_scans(
            {"alignment": {"rt_shift_by_sample": {"sample": 0.0}}},
            scans,
            candidates,
            [component],
        )

        self.assertEqual(accepted, [])
        self.assertTrue(
            not region
            or all(
                row.get("sequence_inference_level")
                == "sequence_region_candidate"
                for row in region
            )
        )

    def test_sequence_inference_catalog_identifies_terminal_truncation(self) -> None:
        candidates = generate_sequence_inference_candidates(
            {"HC": "APEPTIDEK"},
            max_missed_cleavages=0,
            max_terminal_trim=2,
        )
        target = next(
            candidate for candidate in candidates
            if (
                not candidate.is_decoy
                and candidate.sequence == "PEPTIDEK"
                and candidate.proteolysis
                == "n_terminal_truncation_candidate"
            )
        )
        self.assertTrue(any(
            not candidate.is_decoy
            and candidate.proteolysis
            == "c_terminal_truncation_candidate"
            for candidate in candidates
        ))
        component_mass = target.neutral_mass
        fragment_map = {
            label: mz
            for label, mz, _, __ in theoretical_fragments(target)
            if "^" not in label
        }

        def truncation_scan(
            scan_id: str,
            charge: int,
            labels: list[str],
        ) -> LCMSSpectrumScan:
            mz_values = [fragment_map[label] for label in labels]
            return LCMSSpectrumScan(
                scan_id=scan_id,
                raw_file_id="raw",
                sample_id="sample",
                rt=12.0,
                ms_level=2,
                mz_array=mz_values,
                intensity_array=[1500.0] * len(mz_values),
                tic=1500.0 * len(mz_values),
                base_peak_mz=mz_values[0],
                base_peak_intensity=1500.0,
                precursor_mz=component_mass / charge + PROTON,
                precursor_charge=charge,
                precursor_intensity=8000.0,
                activation_method="HCD",
                collision_energy=25.0,
            )

        scans = [
            truncation_scan(
                "trunc=z2",
                2,
                ["b1", "b2", "b3", "b4", "y1", "y2"],
            ),
            truncation_scan(
                "trunc=z3",
                3,
                ["b5", "b6", "y3", "y4", "y5", "y6"],
            ),
        ]
        component = {
            "component_group_id": "MS1INFER_TRUNC",
            "inferred_component": True,
            "component_neutral_mass": component_mass,
            "component_representative_charge": 2,
            "members": [
                {
                    "feature_group_id": "TRUNC_Z2",
                    "true_peak_mz": component_mass / 2 + PROTON,
                    "representative_mz": component_mass / 2 + PROTON,
                    "representative_rt": 12.0,
                    "component_charge": 2,
                    "component_isotope_offset": 0,
                },
                {
                    "feature_group_id": "TRUNC_Z3",
                    "true_peak_mz": component_mass / 3 + PROTON,
                    "representative_mz": component_mass / 3 + PROTON,
                    "representative_rt": 12.0,
                    "component_charge": 3,
                    "component_isotope_offset": 0,
                },
            ],
        }

        accepted, _ = search_component_sequence_tag_scans(
            {"alignment": {"rt_shift_by_sample": {"sample": 0.0}}},
            scans,
            candidates,
            [component],
        )

        self.assertEqual(len(accepted), 1)
        self.assertEqual(accepted[0]["sequence"], "PEPTIDEK")
        self.assertEqual(
            accepted[0]["sequence_inference_level"],
            "truncation_sequence_supported",
        )
        self.assertIn(
            "N-terminal truncation candidate",
            accepted[0]["modification_text"],
        )

    def test_mzml_parser_reads_ms2_precursor_metadata(self) -> None:
        mz = encoded([100.0, 200.0])
        intensity = encoded([10.0, 20.0])
        xml = f"""<mzML><run><spectrumList count="2">
        <spectrum id="scan=1"><cvParam accession="MS:1000511" value="1"/><cvParam accession="MS:1000016" value="1.0"/>
        <binaryDataArrayList><binaryDataArray><cvParam accession="MS:1000514"/><cvParam accession="MS:1000523"/><binary>{mz}</binary></binaryDataArray><binaryDataArray><cvParam accession="MS:1000515"/><cvParam accession="MS:1000523"/><binary>{intensity}</binary></binaryDataArray></binaryDataArrayList></spectrum>
        <spectrum id="scan=2"><cvParam accession="MS:1000511" value="2"/><cvParam accession="MS:1000016" value="1.1"/><cvParam accession="MS:1000285" value="30"/>
        <precursorList><precursor spectrumRef="scan=1"><isolationWindow><cvParam accession="MS:1000828" value="0.8"/><cvParam accession="MS:1000829" value="0.8"/></isolationWindow><selectedIonList><selectedIon><cvParam accession="MS:1000744" value="500.25"/><cvParam accession="MS:1000041" value="2"/><cvParam accession="MS:1000042" value="12345"/></selectedIon></selectedIonList><activation><cvParam accession="MS:1000045" value="25"/><cvParam accession="MS:1000422" value=""/></activation></precursor></precursorList>
        <binaryDataArrayList><binaryDataArray><cvParam accession="MS:1000514"/><cvParam accession="MS:1000523"/><binary>{mz}</binary></binaryDataArray><binaryDataArray><cvParam accession="MS:1000515"/><cvParam accession="MS:1000523"/><binary>{intensity}</binary></binaryDataArray></binaryDataArrayList></spectrum>
        </spectrumList></run></mzML>"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.mzML"
            path.write_text(xml, encoding="utf-8")
            raw_file, scans = read_mzml(path, ms_levels=(2,))
            cached_raw_file, ms1_scans = read_mzml(path, ms_levels=(1,))
            self.assertEqual(cached_raw_file.parser_status, "mzML-cache")
            self.assertEqual(cached_raw_file.scan_count, 1)
            self.assertEqual(ms1_scans[0].ms_level, 1)

        self.assertEqual(raw_file.scan_count, 1)
        self.assertEqual(scans[0].precursor_scan_id, "scan=1")
        self.assertEqual(scans[0].precursor_charge, 2)
        self.assertAlmostEqual(scans[0].precursor_mz or 0.0, 500.25)
        self.assertEqual(scans[0].activation_method, "HCD")
        self.assertEqual(scans[0].collision_energy, 25.0)

    def test_fasta_rejects_empty_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "invalid.fasta"
            path.write_text(">heavy\n>light\nPEPTIDEK\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "no amino acid sequence"):
                read_fasta(path)

    def test_sequence_search_matches_target_fragments(self) -> None:
        candidates = generate_candidates({"HC": "PEPTIDEK"}, max_missed_cleavages=0)
        target = next(candidate for candidate in candidates if not candidate.is_decoy)
        fragments = theoretical_fragments(target)
        mz_array = sorted(ion[1] for ion in fragments[:12])
        scan = LCMSSpectrumScan(
            scan_id="scan=2",
            raw_file_id="synthetic",
            sample_id="reference",
            rt=10.0,
            ms_level=2,
            mz_array=mz_array,
            intensity_array=[1000.0 + index * 100.0 for index in range(len(mz_array))],
            tic=1.0,
            base_peak_mz=mz_array[-1],
            base_peak_intensity=2000.0,
            precursor_mz=target.neutral_mass / 2 + PROTON,
            precursor_charge=2,
            activation_method="HCD",
            collision_energy=25.0,
        )

        psms = search_scans([scan], candidates, precursor_tolerance_ppm=10.0, fragment_tolerance_ppm=10.0)

        self.assertEqual(len(psms), 1)
        self.assertFalse(psms[0]["is_decoy"])
        self.assertEqual(psms[0]["sequence"], "PEPTIDEK")
        self.assertGreaterEqual(psms[0]["matched_ion_count"], 6)
        self.assertGreater(psms[0]["score"], 40.0)

    def test_low_score_target_candidate_keeps_observed_spectrum(self) -> None:
        candidates = generate_candidates({"HC": "PEPTIDEK"}, max_missed_cleavages=0)
        target = next(candidate for candidate in candidates if not candidate.is_decoy)
        scan = LCMSSpectrumScan(
            scan_id="scan=low",
            raw_file_id="synthetic",
            sample_id="reference",
            rt=10.0,
            ms_level=2,
            mz_array=[50.0, 60.0, 70.0],
            intensity_array=[300.0, 200.0, 100.0],
            tic=600.0,
            base_peak_mz=50.0,
            base_peak_intensity=300.0,
            precursor_mz=target.neutral_mass / 2 + PROTON,
            precursor_charge=2,
        )

        psms = search_scans(
            [scan],
            candidates,
            precursor_tolerance_ppm=10.0,
            fragment_tolerance_ppm=10.0,
            min_score=20.0,
        )

        self.assertEqual(len(psms), 1)
        self.assertFalse(psms[0]["is_decoy"])
        self.assertLess(psms[0]["score"], 40.0)
        self.assertEqual(len(psms[0]["spectrum_peaks"]), 3)

    def test_d_level_sequence_is_candidate_not_identified_component(self) -> None:
        precursor_mz = 500.0
        feature = {
            "feature_group_id": "FG_LOW",
            "parent_tic_peak_id": "TICP_0001",
            "representative_mz": precursor_mz,
            "representative_rt": 10.0,
            "ranking_score": 2.0,
            "difference_type": "area_changed",
            "area_by_sample": {"sample": 100.0},
            "normalized_area_by_sample": {"sample": 100.0},
        }
        payload = {
            "alignment": {"rt_shift_by_sample": {"sample": 0.0}},
            "global_feature_groups": [feature],
            "peak_results": [],
        }
        annotation = {
            "feature_or_candidate_id": "FG_LOW",
            "candidate_id": "HC:1-8:PEPTIDEK:Unmodified",
            "base_peptide_id": "HC:1-8:PEPTIDEK",
            "neutral_mass": (precursor_mz - PROTON) * 2,
            "chain": "HC",
            "start": 1,
            "end": 8,
            "sequence": "PEPTIDEK",
            "modification": "Unmodified",
            "confidence": "D_low_evidence",
            "feature_link_type": "direct_precursor",
            "feature_isotope_offset": 0,
            "feature_charge": 2,
            "sample_evidence": {"sample": {"score": 39.0}},
        }
        psm = {
            "sample_id": "sample",
            "scan_id": "scan=low",
            "rt": 10.0,
            "precursor_mz": precursor_mz,
            "precursor_charge": 2,
            "sequence": "PEPTIDEK",
            "modification_text": "Unmodified",
            "score": 39.0,
            "q_value": 0.005,
            "matched_ion_count": 8,
            "fragment_coverage": 0.2,
            "spectrum_peaks": [{"mz": 100.0, "intensity": 50.0, "label": "b1"}],
            "feature_links": [{"feature_group_id": "FG_LOW"}],
        }
        scan = LCMSSpectrumScan(
            scan_id="scan=low",
            raw_file_id="synthetic",
            sample_id="sample",
            rt=10.0,
            ms_level=2,
            mz_array=[100.0],
            intensity_array=[50.0],
            tic=50.0,
            base_peak_mz=100.0,
            base_peak_intensity=50.0,
            precursor_mz=precursor_mz,
            precursor_charge=2,
            isolation_window_lower_offset=0.8,
            isolation_window_upper_offset=0.8,
        )

        rows, summary = build_feature_ms2_evidence(
            payload, [scan], [psm], [annotation]
        )
        self.assertEqual(rows[0]["ms2_status"], "low_evidence_sequence_candidate")
        self.assertEqual(
            rows[0]["unidentified_reason"],
            "low_evidence_sequence_candidate_below_identification_threshold",
        )
        self.assertEqual(summary["low_evidence_sequence_candidate"], 1)

        groups = build_ms1_component_groups(payload, [annotation])
        self.assertFalse(groups[0]["identified_component"])
        self.assertTrue(groups[0]["candidate_component"])
        self.assertEqual(groups[0]["component_candidate_sequence"], "PEPTIDEK")
        self.assertIn("candidate PEPTIDEK", groups[0]["component_label"])

    def test_feature_ms2_evidence_indexes_non_significant_features(self) -> None:
        payload = {
            "alignment": {"rt_shift_by_sample": {"sample": 0.0}},
            "global_feature_groups": [
                {
                    "feature_group_id": "FG_COMMON",
                    "representative_mz": 500.0,
                    "representative_rt": 10.0,
                    "difference_type": "common_feature",
                    "ranking_score": 0.1,
                }
            ],
        }

        rows, summary = build_feature_ms2_evidence(payload, [], [], [])

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["feature_group_id"], "FG_COMMON")
        self.assertEqual(rows[0]["ms2_status"], "no_ms2_acquired")
        self.assertEqual(rows[0]["unidentified_reason"], "no_ms2_scan")
        self.assertEqual(summary["no_ms2_acquired"], 1)

    def test_feature_open_mass_search_recovers_known_backbone_with_unknown_delta(self) -> None:
        candidates = generate_sequence_inference_candidates(
            {"HC": "PEPTIDEK"},
            max_missed_cleavages=0,
            max_terminal_trim=0,
        )
        target = next(candidate for candidate in candidates if not candidate.is_decoy)
        delta = 15.99491462
        component_mass = target.neutral_mass + delta
        fragment_map = {
            label: mz for label, mz, _, __ in theoretical_fragments(target)
            if "^" not in label
        }
        site = 3
        mz_values = []
        for label in ("b2", "b3", "b4", "b5", "y2", "y5", "y6"):
            mz = fragment_map[label]
            ordinal = int(label[1:])
            shifted = (
                label.startswith("b") and ordinal > site
            ) or (
                label.startswith("y") and ordinal >= len(target.sequence) - site
            )
            mz_values.append(mz + (delta if shifted else 0.0))
        scan = LCMSSpectrumScan(
            scan_id="scan=open",
            raw_file_id="raw",
            sample_id="sample",
            rt=10.0,
            ms_level=2,
            mz_array=sorted(mz_values),
            intensity_array=[1000.0] * len(mz_values),
            tic=5000.0,
            base_peak_mz=max(mz_values),
            base_peak_intensity=1000.0,
            precursor_mz=component_mass / 2 + PROTON,
            precursor_charge=2,
            precursor_intensity=10000.0,
            isolation_window_lower_offset=0.8,
            isolation_window_upper_offset=0.8,
        )
        payload = {
            "alignment": {"rt_shift_by_sample": {"sample": 0.0}},
            "global_feature_groups": [{
                "feature_group_id": "FG_OPEN",
                "representative_mz": component_mass / 2 + PROTON,
                "true_peak_mz": component_mass / 2 + PROTON,
                "representative_rt": 10.0,
                "difference_type": "area_changed",
                "ranking_score": 2.0,
            }],
        }
        matches = search_feature_open_mass_scans(
            payload,
            [scan],
            candidates,
            fragment_tolerance_ppm=20.0,
            precursor_tolerance_ppm=20.0,
        )
        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["sequence"], "PEPTIDEK")
        self.assertEqual(matches[0]["feature_link_type"], "feature_open_mass_search")
        self.assertAlmostEqual(float(matches[0]["mass_delta"]), delta, places=3)
        self.assertTrue(matches[0]["spectrum_peaks"])
        self.assertEqual(matches[0]["modification_text"], "Open ΔMass +15.9949 Da")

    def test_feature_targeted_glycopeptide_search_uses_diagnostic_and_y_ions(self) -> None:
        candidates = generate_candidates(
            {"HC": "EEQYNSTYR"},
            max_missed_cleavages=0,
            max_variable_modifications=0,
        )
        target = next(candidate for candidate in candidates if not candidate.is_decoy)
        glycan = next(
            row for row in TARGETED_N_GLYCANS
            if row["name"] == "G0F N-glycan"
        )
        glycan_mass = float(glycan["mass"])
        component_mass = target.neutral_mass + glycan_mass
        fragment_map = {
            label: mz for label, mz, _, __ in theoretical_fragments(target)
            if "^" not in label
        }
        pairs = [
            (138.0550, 2200.0),
            (204.086649, 6000.0),
            (366.139472, 4200.0),
            (target.neutral_mass + PROTON, 1800.0),
            (
                target.neutral_mass + 203.0793725330 + PROTON,
                2600.0,
            ),
            (
                target.neutral_mass + 2 * 203.0793725330 + PROTON,
                2300.0,
            ),
        ]
        for label in ("b2", "b3", "b4", "y2", "y3", "y4"):
            pairs.append((fragment_map[label], 1500.0))
        for label in ("b5", "b6", "y5"):
            pairs.append((fragment_map[label] + 203.0793725330, 1700.0))
        pairs.sort(key=lambda pair: pair[0])
        scan = LCMSSpectrumScan(
            scan_id="scan=glyco",
            raw_file_id="raw",
            sample_id="sample",
            rt=10.0,
            ms_level=2,
            mz_array=[pair[0] for pair in pairs],
            intensity_array=[pair[1] for pair in pairs],
            tic=sum(pair[1] for pair in pairs),
            base_peak_mz=204.086649,
            base_peak_intensity=6000.0,
            precursor_mz=component_mass / 2 + PROTON,
            precursor_charge=2,
            precursor_intensity=10000.0,
            isolation_window_lower_offset=0.8,
            isolation_window_upper_offset=0.8,
            activation_method="HCD",
        )
        payload = {
            "alignment": {"rt_shift_by_sample": {"sample": 0.0}},
            "global_feature_groups": [{
                "feature_group_id": "FG_GLYCO",
                "representative_mz": component_mass / 2 + PROTON,
                "true_peak_mz": component_mass / 2 + PROTON,
                "representative_rt": 10.0,
                "difference_type": "area_changed",
                "ranking_score": 2.0,
            }],
        }

        matches = search_feature_glycopeptide_scans(
            payload,
            [scan],
            candidates,
        )

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["sequence"], "EEQYNSTYR")
        self.assertEqual(matches[0]["glycan_name"], "G0F N-glycan")
        self.assertEqual(
            matches[0]["feature_link_type"],
            "feature_targeted_glycopeptide_search",
        )
        self.assertGreaterEqual(matches[0]["glycan_diagnostic_ion_count"], 3)
        self.assertGreaterEqual(matches[0]["glycan_core_y_ion_count"], 3)
        self.assertGreaterEqual(matches[0]["glycan_hexnac_fragment_count"], 1)

    def test_feature_open_mass_search_skips_features_without_ms2(self) -> None:
        payload = {
            "alignment": {"rt_shift_by_sample": {"sample": 0.0}},
            "global_feature_groups": [{
                "feature_group_id": "FG_NO_MS2",
                "representative_mz": 500.0,
                "representative_rt": 10.0,
                "difference_type": "area_changed",
            }],
        }
        candidates = generate_sequence_inference_candidates(
            {"HC": "PEPTIDEK"}, max_missed_cleavages=0, max_terminal_trim=0
        )
        self.assertEqual(
            search_feature_open_mass_scans(payload, [], candidates),
            [],
        )
        self.assertEqual(
            search_global_open_modification_scans([], candidates),
            ([], [], []),
        )

    def test_global_open_search_discovers_and_rescores_recurrent_oxidation(self) -> None:
        candidates = generate_sequence_inference_candidates(
            {"HC": "PEPMIDEK"},
            max_missed_cleavages=0,
            max_terminal_trim=0,
        )
        target = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and candidate.sequence == "PEPMIDEK"
        )
        delta = 15.99491462
        site = 3
        fragments = {
            label: mz for label, mz, _, __ in theoretical_fragments(target)
            if "^" not in label
        }

        def shifted_mz(label: str) -> float:
            ordinal = int(label[1:])
            contains_site = (
                (label.startswith("b") and ordinal > site)
                or (
                    label.startswith("y")
                    and ordinal >= len(target.sequence) - site
                )
            )
            return fragments[label] + (delta if contains_site else 0.0)

        labels = [
            *(f"b{index}" for index in range(1, len(target.sequence))),
            *(f"y{index}" for index in range(1, len(target.sequence))),
        ]
        mz_values = sorted(shifted_mz(label) for label in labels)
        component_mass = target.neutral_mass + delta
        scans = [
            LCMSSpectrumScan(
                scan_id=f"scan=global-open-{index}",
                raw_file_id="raw",
                sample_id="sample_a" if index < 2 else "sample_b",
                rt=10.0 + index * 0.25,
                ms_level=2,
                mz_array=mz_values,
                intensity_array=[1500.0 + peak * 75.0 for peak in range(len(mz_values))],
                tic=30000.0,
                base_peak_mz=mz_values[-1],
                base_peak_intensity=3000.0,
                precursor_mz=component_mass / 2 + PROTON,
                precursor_charge=2,
                precursor_intensity=20000.0,
                activation_method="HCD",
                collision_energy=25.0,
            )
            for index in range(4)
        ]

        accepted, tentative, clusters = search_global_open_modification_scans(
            scans,
            candidates,
            max_discovery_scans=10,
            max_delta_clusters=5,
        )
        self.assertTrue(accepted)
        self.assertFalse(tentative)
        self.assertTrue(clusters)
        self.assertEqual(accepted[0]["sequence"], "PEPMIDEK")
        self.assertAlmostEqual(float(accepted[0]["mass_delta"]), delta, places=3)
        self.assertEqual(
            accepted[0]["open_modification_confidence"],
            "B_open_modification_supported",
        )
        self.assertIn("Oxidation", accepted[0]["modification_text"])
        self.assertEqual(clusters[0]["classification"], "supported_open_modification")
        self.assertGreaterEqual(int(clusters[0]["target_psm_count"]), 2)
        known = build_global_known_modification_summary(
            accepted,
            tentative,
            clusters,
        )
        self.assertEqual(known[0]["modification"], "Oxidation")
        self.assertEqual(
            known[0]["classification"],
            "supported_known_modification",
        )
        self.assertEqual(known[0]["sample_count"], 2)
        self.assertIn("HC:M4", known[0]["localized_sites"])
        findings = build_modified_peptide_findings(accepted, [])
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0]["confidence"], "B_high_confidence_inferred")
        self.assertAlmostEqual(
            float(findings[0]["modification_mass_delta"]),
            delta,
            places=3,
        )

    def test_candidate_catalog_includes_antibody_quality_attributes(self) -> None:
        candidates = generate_candidates({"HC": "EEQYNSTYRSLSLSPGK"}, max_missed_cleavages=1)
        target_modifications = {
            candidate.modification_text
            for candidate in candidates
            if not candidate.is_decoy
        }

        self.assertTrue(any(text.startswith("G0F N-glycan@") for text in target_modifications))
        self.assertTrue(any(text.startswith("Glycation@") for text in target_modifications))
        self.assertTrue(any("C-terminal Lys clipping" in text for text in target_modifications))

    def test_expanded_catalog_includes_double_modifications_and_semitryptic_peptides(self) -> None:
        candidates = generate_candidates(
            {"HC": "PEPMIDEK"},
            max_missed_cleavages=0,
            max_variable_modifications=2,
            semitryptic_max_trim=1,
        )
        targets = [candidate for candidate in candidates if not candidate.is_decoy]

        self.assertTrue(any(len(candidate.modifications) == 2 for candidate in targets))
        self.assertTrue(any(candidate.proteolysis == "semi_tryptic" for candidate in targets))

    def test_modified_findings_include_unpaired_modified_psm(self) -> None:
        candidates = generate_candidates({"HC": "PEPMIDEK"}, max_missed_cleavages=0)
        modified = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and candidate.modification_text == "Oxidation@4"
        )
        psm = {
            **modified.__dict__,
            "sample_id": "test",
            "scan_id": "scan=10",
            "rt": 10.0,
            "precursor_mz": modified.neutral_mass / 2 + PROTON,
            "precursor_charge": 2,
            "precursor_intensity": 10000.0,
            "score": 70.0,
            "q_value": 0.001,
            "matched_ion_count": 8,
            "fragment_coverage": 0.4,
            "explained_intensity": 0.3,
            "modification_text": modified.modification_text,
            "feature_links": [{"feature_group_id": "FG_MOD"}],
        }
        payload = {
            "global_feature_groups": [{
                "feature_group_id": "FG_MOD",
                "difference_type": "area_changed",
                "ranking_score": 2.0,
                "max_fold_change": 4.0,
                "higher_abundance_sample": "test",
                "normalized_area_by_sample": {"test": 100.0},
            }],
        }

        findings = build_modified_peptide_findings([psm], [], payload, [])

        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0]["pairing_status"], "modified_identified_without_unmodified_ms2")
        self.assertEqual(findings[0]["evidence_scope"], "differential_feature_linked")
        self.assertEqual(findings[0]["confidence"], "B_high_confidence_inferred")

    def test_repeated_selected_ms2_scans_can_form_tentative_consensus_match(self) -> None:
        candidates = generate_candidates({"HC": "PEPTIDEK"}, max_missed_cleavages=0)
        target = next(candidate for candidate in candidates if not candidate.is_decoy and not candidate.modifications)
        fragment_mz = sorted(ion[1] for ion in theoretical_fragments(target))
        precursor_mz = target.neutral_mass / 2 + PROTON
        scans = [
            LCMSSpectrumScan(
                scan_id=f"scan={index}",
                raw_file_id="sample",
                sample_id="sample",
                rt=10.0 + index * 0.02,
                ms_level=2,
                mz_array=fragment_mz,
                intensity_array=[1000.0 + point * 10.0 for point in range(len(fragment_mz))],
                tic=1.0,
                base_peak_mz=fragment_mz[-1],
                base_peak_intensity=2000.0,
                precursor_mz=precursor_mz,
                precursor_charge=2,
                precursor_intensity=10000.0,
                isolation_window_lower_offset=0.8,
                isolation_window_upper_offset=0.8,
                activation_method="HCD",
                collision_energy=25.0,
            )
            for index in range(2)
        ]
        payload = {
            "alignment": {"rt_shift_by_sample": {"sample": 0.0}},
            "global_feature_groups": [{
                "feature_group_id": "FG_CONSENSUS",
                "representative_mz": precursor_mz,
                "representative_rt": 10.0,
                "difference_type": "area_changed",
                "ranking_score": 2.0,
            }],
        }

        matches = search_selected_feature_consensus_scans(payload, scans, candidates)

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["sequence"], "PEPTIDEK")
        self.assertEqual(matches[0]["feature_link_type"], "selected_precursor_consensus_search")
        self.assertEqual(matches[0]["consensus_source_scan_count"], 2)

    def test_component_search_combines_charge_state_fragment_evidence(self) -> None:
        candidates = generate_candidates({"HC": "PEPTIDEK"}, max_missed_cleavages=0)
        target = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and not candidate.modifications
        )
        fragments = {
            label: mz
            for label, mz, _, __ in theoretical_fragments(target)
        }
        z2_mz = target.neutral_mass / 2 + PROTON
        z3_mz = target.neutral_mass / 3 + PROTON

        def scan(scan_id: str, precursor_mz: float, charge: int, labels: list[str]) -> LCMSSpectrumScan:
            mz_values = sorted(fragments[label] for label in labels)
            return LCMSSpectrumScan(
                scan_id=scan_id,
                raw_file_id="sample",
                sample_id="sample",
                rt=10.0,
                ms_level=2,
                mz_array=mz_values,
                intensity_array=[1000.0 + index * 100.0 for index in range(len(mz_values))],
                tic=sum(1000.0 + index * 100.0 for index in range(len(mz_values))),
                base_peak_mz=mz_values[-1],
                base_peak_intensity=2000.0,
                precursor_mz=precursor_mz,
                precursor_charge=charge,
                precursor_intensity=10000.0,
                isolation_window_lower_offset=0.8,
                isolation_window_upper_offset=0.8,
                activation_method="HCD",
                collision_energy=25.0,
            )

        scans = [
            scan("scan=z2", z2_mz, 2, ["b1", "b2", "b3", "y1", "y2"]),
            scan("scan=z3", z3_mz, 3, ["b4", "b5", "b6", "y3", "y4"]),
        ]
        payload = {
            "alignment": {"rt_shift_by_sample": {"sample": 0.0}},
            "global_feature_groups": [],
        }
        component = {
            "component_group_id": "MS1INFER_0001",
            "inferred_component": True,
            "component_confidence": "MS1_high",
            "component_neutral_mass": target.neutral_mass,
            "component_representative_charge": 2,
            "members": [
                {
                    "feature_group_id": "FG_Z2",
                    "representative_mz": z2_mz,
                    "true_peak_mz": z2_mz,
                    "representative_rt": 10.0,
                    "component_charge": 2,
                    "component_isotope_offset": 0,
                    "difference_type": "area_changed",
                    "ranking_score": 2.0,
                },
                {
                    "feature_group_id": "FG_Z3",
                    "representative_mz": z3_mz,
                    "true_peak_mz": z3_mz,
                    "representative_rt": 10.0,
                    "component_charge": 3,
                    "component_isotope_offset": 0,
                    "difference_type": "area_changed",
                    "ranking_score": 1.8,
                },
            ],
        }

        matches = search_component_consensus_scans(
            payload,
            scans,
            candidates,
            [component],
            component_mass_tolerance_ppm=20.0,
            fragment_tolerance_ppm=20.0,
        )

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["sequence"], "PEPTIDEK")
        self.assertEqual(
            matches[0]["feature_link_type"],
            "component_charge_isotope_consensus_search",
        )
        self.assertEqual(matches[0]["component_member_support_count"], 2)
        self.assertEqual(matches[0]["component_charge_support_count"], 2)
        self.assertGreaterEqual(matches[0]["matched_ion_count"], 8)
        self.assertEqual(
            {link["feature_group_id"] for link in matches[0]["feature_links"]},
            {"FG_Z2", "FG_Z3"},
        )
        payload["global_feature_groups"] = [dict(member) for member in component["members"]]
        annotations = differential_annotations(matches)
        evidence_rows, status_summary = build_feature_ms2_evidence(
            payload,
            scans,
            matches,
            annotations,
            candidates=candidates,
        )
        self.assertEqual(
            {row["ms2_status"] for row in evidence_rows},
            {"tentative_component_consensus_identification"},
        )
        self.assertEqual(
            status_summary["tentative_component_consensus_identification"],
            2,
        )

    def test_isotope_envelope_links_one_psm_to_multiple_ms1_features(self) -> None:
        precursor_mz = 500.25
        isotope_mz = precursor_mz + ISOTOPE_MASS_DIFF / 2
        neutral_mass = (precursor_mz - PROTON) * 2
        charge_three_mz = neutral_mass / 3 + PROTON
        payload = {
            "alignment": {"rt_shift_by_sample": {"sample": 0.0}},
            "global_feature_groups": [
                {
                    "feature_group_id": "FG_MONO",
                    "representative_mz": precursor_mz,
                    "representative_rt": 10.0,
                    "difference_type": "area_changed",
                    "difference_score": 0.8,
                    "ranking_score": 2.0,
                    "area_by_sample": {"sample": 100.0},
                    "normalized_area_by_sample": {"sample": 100.0},
                },
                {
                    "feature_group_id": "FG_ISOTOPE",
                    "representative_mz": isotope_mz,
                    "representative_rt": 10.0,
                    "difference_type": "area_changed",
                    "difference_score": 0.7,
                    "ranking_score": 1.8,
                    "area_by_sample": {"sample": 60.0},
                    "normalized_area_by_sample": {"sample": 60.0},
                },
                {
                    "feature_group_id": "FG_CHARGE3",
                    "representative_mz": charge_three_mz,
                    "representative_rt": 10.0,
                    "difference_type": "moderate_difference",
                    "difference_score": 0.5,
                    "ranking_score": 1.5,
                    "area_by_sample": {"sample": 40.0},
                    "normalized_area_by_sample": {"sample": 40.0},
                },
            ],
            "peak_results": [],
        }
        psm = {
            "sample_id": "sample",
            "scan_id": "scan=2",
            "rt": 10.0,
            "precursor_mz": precursor_mz,
            "precursor_charge": 2,
            "candidate_id": "HC:1-8:PEPTIDEK:Unmodified",
            "base_peptide_id": "HC:1-8:PEPTIDEK",
            "chain": "HC",
            "start": 1,
            "end": 8,
            "sequence": "PEPTIDEK",
            "modification_text": "Unmodified",
            "score": 75.0,
            "q_value": 0.0,
            "matched_ion_count": 8,
            "fragment_coverage": 0.5,
        }

        annotate_feature_groups([psm], payload)
        annotations = differential_annotations([psm])

        self.assertEqual(
            {link["feature_group_id"] for link in psm["feature_links"]},
            {"FG_MONO", "FG_ISOTOPE", "FG_CHARGE3"},
        )
        isotope_annotation = next(row for row in annotations if row["feature_or_candidate_id"] == "FG_ISOTOPE")
        self.assertEqual(isotope_annotation["feature_link_type"], "isotope_envelope")
        self.assertIn("不是独立修饰", isotope_annotation["interpretation"])
        charge_annotation = next(row for row in annotations if row["feature_or_candidate_id"] == "FG_CHARGE3")
        self.assertEqual(charge_annotation["feature_link_type"], "charge_state_envelope")

        component_groups = build_ms1_component_groups(payload, annotations)
        self.assertEqual(len(component_groups), 1)
        self.assertTrue(component_groups[0]["identified_component"])
        self.assertEqual(component_groups[0]["component_member_count"], 3)
        self.assertEqual(component_groups[0]["component_primary_feature_id"], "FG_MONO")
        self.assertEqual(component_groups[0]["component_charge_states"], [2, 3])
        self.assertAlmostEqual(component_groups[0]["area_by_sample"]["sample"], 200.0)
        self.assertAlmostEqual(component_groups[0]["ranking_score"], 2.0)

        scan = LCMSSpectrumScan(
            scan_id="scan=2",
            raw_file_id="synthetic",
            sample_id="sample",
            rt=10.0,
            ms_level=2,
            mz_array=[100.0, 200.0],
            intensity_array=[10.0, 20.0],
            tic=30.0,
            base_peak_mz=200.0,
            base_peak_intensity=20.0,
            precursor_mz=precursor_mz,
            precursor_charge=2,
            isolation_window_lower_offset=0.8,
            isolation_window_upper_offset=0.8,
        )
        rows, summary = build_feature_ms2_evidence(payload, [scan], [psm], annotations)
        isotope_row = next(row for row in rows if row["feature_group_id"] == "FG_ISOTOPE")
        self.assertEqual(isotope_row["ms2_status"], "identified_isotope_envelope")
        self.assertEqual(summary["identified_isotope_envelope"], 1)
        charge_row = next(row for row in rows if row["feature_group_id"] == "FG_CHARGE3")
        self.assertEqual(charge_row["ms2_status"], "identified_charge_state_envelope")

    def test_unidentified_ms1_features_merge_by_isotope_charge_and_xic(self) -> None:
        mono_mz = 500.0
        neutral_mass = (mono_mz - PROTON) * 2
        isotope_mz = mono_mz + ISOTOPE_MASS_DIFF / 2
        charge_three_mz = neutral_mass / 3 + PROTON
        sample_ids = ["reference", "test"]

        def feature(feature_id: str, mz: float, reference_area: float, test_area: float) -> dict[str, object]:
            return {
                "feature_group_id": feature_id,
                "parent_tic_peak_id": "TICP_0001",
                "representative_mz": mz,
                "representative_rt": 10.0,
                "ranking_score": test_area / 100.0,
                "difference_type": "moderate_difference",
                "area_by_sample": {"reference": reference_area, "test": test_area},
                "normalized_area_by_sample": {"reference": reference_area, "test": test_area},
                "max_area": test_area,
            }

        payload = {
            "sample_ids": sample_ids,
            "global_feature_groups": [
                feature("FG_MONO_Z2", mono_mz, 100.0, 220.0),
                feature("FG_ISOTOPE_Z2", isotope_mz, 60.0, 132.0),
                feature("FG_MONO_Z3", charge_three_mz, 80.0, 176.0),
                feature("FG_UNRELATED", 700.0, 100.0, 20.0),
            ],
            "peak_results": [],
        }

        spectra_by_sample: dict[str, list[dict[str, object]]] = {}
        for sample_id in sample_ids:
            scans: list[dict[str, object]] = []
            scale = 2.2 if sample_id == "test" else 1.0
            for index in range(31):
                rt = 9.7 + index * 0.02
                shape = math.exp(-((rt - 10.0) / 0.09) ** 2)
                unrelated_shape = math.exp(-((rt - 10.25) / 0.05) ** 2)
                scans.append({
                    "aligned_rt": rt,
                    "rt": rt,
                    "mz": [charge_three_mz, mono_mz, isotope_mz, 700.0],
                    "intensity": [
                        800.0 * shape * scale,
                        1000.0 * shape * scale,
                        600.0 * shape * scale,
                        900.0 * unrelated_shape,
                    ],
                })
            spectra_by_sample[sample_id] = scans

        groups = build_ms1_component_groups(
            payload,
            [],
            spectra_loader=lambda sample_id: spectra_by_sample[sample_id],
        )

        inferred = [row for row in groups if row.get("inferred_component")]
        self.assertEqual(len(inferred), 1)
        self.assertEqual(inferred[0]["component_member_count"], 3)
        self.assertEqual(inferred[0]["component_charge_states"], [2, 3])
        self.assertEqual(
            {member["feature_group_id"] for member in inferred[0]["members"]},
            {"FG_MONO_Z2", "FG_ISOTOPE_Z2", "FG_MONO_Z3"},
        )
        self.assertAlmostEqual(inferred[0]["component_neutral_mass"], neutral_mass, places=3)
        self.assertTrue(any(row["component_primary_feature_id"] == "FG_UNRELATED" for row in groups))

    def test_isotope_envelope_refinement_merges_coarse_z2_z3_features(self) -> None:
        sample_ids = ["reference", "test"]
        z3_feature_mz = 627.9658
        z2_feature_mz = 941.2812
        z3_envelope = [627.687866, 628.021789, 628.355835]
        z2_envelope = [941.027954, 941.529114]

        def global_feature(feature_id: str, mz: float, rt: float) -> dict[str, object]:
            return {
                "feature_group_id": feature_id,
                "parent_tic_peak_id": "TICP_0074",
                "representative_mz": mz,
                "representative_rt": rt,
                "ranking_score": 8.0,
                "difference_type": "presence_absence",
                "area_by_sample": {"reference": 100.0, "test": 1_000_000.0},
                "normalized_area_by_sample": {"reference": 100.0, "test": 1_000_000.0},
                "max_area": 1_000_000.0,
                "source_feature_group_ids": [feature_id],
            }

        source_groups = []
        for feature_id, mz in (("FG_Z3", z3_feature_mz), ("FG_Z2", z2_feature_mz)):
            source_groups.append({
                "feature_group_id": feature_id,
                "features_by_sample": {
                    sample_id: {
                        "aligned_rt_apex": 10.0,
                        "rt_start": 9.85,
                        "rt_end": 10.15,
                        "normalized_area": 1_000_000.0 if sample_id == "test" else 100.0,
                    }
                    for sample_id in sample_ids
                },
            })
        payload = {
            "sample_ids": sample_ids,
            "global_feature_groups": [
                global_feature("FG_Z3", z3_feature_mz, 10.00),
                global_feature("FG_Z2", z2_feature_mz, 10.28),
            ],
            "peak_results": [{"feature_groups": source_groups}],
        }
        spectra_by_sample: dict[str, list[dict[str, object]]] = {}
        for sample_id in sample_ids:
            scans = []
            scale = 100.0 if sample_id == "test" else 1.0
            for index in range(31):
                rt = 9.7 + index * 0.02
                shape = math.exp(-((rt - 10.0) / 0.09) ** 2) * scale
                scans.append({
                    "scan_id": f"{sample_id}_{index}",
                    "aligned_rt": rt,
                    "rt": rt,
                    "mz": z3_envelope + z2_envelope,
                    "intensity": [
                        800.0 * shape,
                        1000.0 * shape,
                        600.0 * shape,
                        900.0 * shape,
                        850.0 * shape,
                    ],
                })
            spectra_by_sample[sample_id] = scans

        groups = build_ms1_component_groups(
            payload,
            [],
            spectra_loader=lambda sample_id: spectra_by_sample[sample_id],
        )

        inferred = next(row for row in groups if row.get("inferred_component"))
        self.assertEqual(inferred["component_member_count"], 2)
        self.assertEqual(inferred["component_charge_states"], [2, 3])
        self.assertIn("resolved_isotope_envelopes", inferred["component_mass_evidence"])
        self.assertIn("averagine_full_isotope_envelope", inferred["component_mass_evidence"])
        self.assertGreater(inferred["component_isotope_fit_score"], 0.85)
        self.assertAlmostEqual(inferred["component_neutral_mass"], 1880.0415, places=2)
        self.assertIsNotNone(inferred["component_representative_mz"])
        self.assertIn(inferred["component_representative_charge"], {2, 3})
        for member in inferred["members"]:
            self.assertIn("true_peak_mz", member)
            self.assertIn("envelope_representative_mz", member)
            self.assertIn("component_isotope_fit_score", member)
            self.assertIn("component_theoretical_isotope_intensity", member)

    def test_modified_peptide_pairs_unmodified_form_and_merges_charge_states(self) -> None:
        candidates = generate_candidates({"HC": "PEPMIDEK"}, max_missed_cleavages=0)
        modified = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and candidate.modification_text == "Oxidation@4"
        )
        unmodified = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and candidate.modification_text == "Unmodified"
        )
        dioxidized = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and candidate.modification_text == "Dioxidation@4"
        )
        sample_ids = ["reference", "test"]
        payload = {
            "sample_ids": sample_ids,
            "reference_sample": "reference",
            "alignment": {"rt_shift_by_sample": {sample: 0.0 for sample in sample_ids}},
            "feature_area_normalization": {"factor_by_sample": {sample: 1.0 for sample in sample_ids}},
            "global_feature_groups": [
                {"feature_group_id": "FG_MOD", "difference_type": "area_changed"},
                {"feature_group_id": "FG_DIOX", "difference_type": "common_feature"},
            ],
        }
        annotations = [{
            "feature_or_candidate_id": "FG_MOD",
            "candidate_id": modified.candidate_id,
            "feature_rt": 10.0,
            "modification": modified.modification_text,
            "confidence": "B_high_confidence_inferred",
        }]
        psms = [
            {"candidate_id": modified.candidate_id, "sample_id": sample, "rt": 10.0}
            for sample in sample_ids
        ] + [
            {"candidate_id": unmodified.candidate_id, "sample_id": sample, "rt": 12.0}
            for sample in sample_ids
        ] + [
            {
                "candidate_id": dioxidized.candidate_id,
                "sample_id": sample,
                "rt": 11.0,
                "precursor_charge": 2,
                "score": 80.0,
                "q_value": 0.0,
                "matched_ion_count": 8,
                "fragment_coverage": 0.5,
                "feature_links": [{
                    "feature_group_id": "FG_DIOX",
                    "feature_rt": 11.0,
                    "feature_ranking": 0.1,
                    "feature_difference_type": "common_feature",
                }],
            }
            for sample in sample_ids
        ]

        def spectra(sample_id: str) -> list[dict[str, object]]:
            scans: list[dict[str, object]] = []
            for index in range(121):
                rt = 8.0 + index * 0.05
                points: list[tuple[float, float]] = []
                for candidate, apex, amplitude in (
                    (modified, 10.0, 8000.0 if sample_id == "test" else 2000.0),
                    (dioxidized, 11.0, 3000.0),
                    (unmodified, 12.0, 8000.0),
                ):
                    peak = amplitude * (2.718281828 ** (-0.5 * ((rt - apex) / 0.12) ** 2))
                    for charge in (2, 3):
                        for isotope_offset in range(3):
                            mz = (candidate.neutral_mass + isotope_offset * ISOTOPE_MASS_DIFF) / charge + PROTON
                            points.append((mz, peak * (1.0 - isotope_offset * 0.2)))
                points.sort()
                scans.append({
                    "scan_id": f"scan={index}",
                    "rt": rt,
                    "aligned_rt": rt,
                    "mz": [mz for mz, _ in points],
                    "intensity": [intensity for _, intensity in points],
                })
            return scans

        rows = build_peptide_form_comparisons(
            payload,
            psms,
            annotations,
            candidates,
            spectra,
        )

        self.assertEqual(len(rows), 2)
        row = next(item for item in rows if item["modification"] == "Oxidation@4")
        recalled = next(item for item in rows if item["modification"] == "Dioxidation@4")
        self.assertIn("FG_DIOX", recalled["linked_feature_ids"])
        self.assertEqual(row["modified_form"]["status"], "ms2_confirmed")
        self.assertEqual(row["unmodified_form"]["status"], "ms2_confirmed")
        self.assertEqual(row["modified_form"]["charge_states_by_sample"]["reference"], [2, 3])
        comparison = row["pairwise_comparisons"][0]
        self.assertEqual(comparison["modification_direction"], "increased")
        self.assertEqual(comparison["modified_charge_trend"]["status"], "consistent")
        self.assertAlmostEqual(comparison["modified_unmodified_ratio_fold"], 4.0, places=3)

    def test_terminal_lys_family_is_triggered_by_retained_form_difference(self) -> None:
        candidates = generate_candidates({"HC": "PEPTIDEK"}, max_missed_cleavages=0)
        retained = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and candidate.sequence == "PEPTIDEK" and not candidate.modifications
        )
        clipped = next(
            candidate for candidate in candidates
            if not candidate.is_decoy and "C-terminal Lys clipping" in candidate.modification_text
        )
        sample_ids = ["reference", "test"]
        payload = {
            "sample_ids": sample_ids,
            "reference_sample": "reference",
            "alignment": {"rt_shift_by_sample": {sample: 0.0 for sample in sample_ids}},
            "feature_area_normalization": {"factor_by_sample": {sample: 1.0 for sample in sample_ids}},
            "global_feature_groups": [
                {"feature_group_id": "FG_RETAINED", "difference_type": "area_changed"},
                {"feature_group_id": "FG_CLIPPED", "difference_type": "common_feature"},
            ],
        }
        retained_psms = [{
            "candidate_id": retained.candidate_id,
            "base_peptide_id": retained.base_peptide_id,
            "sample_id": sample,
            "scan_id": f"{sample}-retained",
            "rt": 10.0,
            "precursor_charge": 2,
            "chain": retained.chain,
            "start": retained.start,
            "end": retained.end,
            "sequence": retained.sequence,
            "neutral_mass": retained.neutral_mass,
            "modifications": [],
            "modification_text": "Unmodified",
            "proteolysis": "fully_tryptic",
            "score": 90.0,
            "q_value": 0.0,
            "matched_ion_count": 12,
            "fragment_coverage": 0.7,
            "feature_links": [{
                "feature_group_id": "FG_RETAINED",
                "feature_rt": 10.0,
                "feature_ranking": 2.0,
                "feature_difference_type": "area_changed",
            }],
        } for sample in sample_ids]
        clipped_psms = [{
            "candidate_id": "HC:1-7:PEPTIDE:Unmodified",
            "base_peptide_id": "HC:1-7:PEPTIDE",
            "sample_id": sample,
            "scan_id": f"{sample}-clipped",
            "rt": 12.0,
            "precursor_charge": 2,
            "chain": clipped.chain,
            "start": clipped.start,
            "end": clipped.end,
            "sequence": clipped.sequence,
            "neutral_mass": clipped.neutral_mass,
            "modifications": [],
            "modification_text": "Unmodified",
            "proteolysis": "semi_tryptic",
            "score": 88.0,
            "q_value": 0.0,
            "matched_ion_count": 10,
            "fragment_coverage": 0.65,
            "feature_links": [{
                "feature_group_id": "FG_CLIPPED",
                "feature_rt": 12.0,
                "feature_ranking": 0.1,
                "feature_difference_type": "common_feature",
            }],
        } for sample in sample_ids]
        psms = retained_psms + clipped_psms
        annotations = differential_annotations(psms)

        def spectra(sample_id: str) -> list[dict[str, object]]:
            scans: list[dict[str, object]] = []
            for index in range(121):
                rt = 8.0 + index * 0.05
                points: list[tuple[float, float]] = []
                for candidate, apex, amplitude in (
                    (retained, 10.0, 8000.0 if sample_id == "test" else 2000.0),
                    (clipped, 12.0, 8000.0),
                ):
                    peak = amplitude * math.exp(-0.5 * ((rt - apex) / 0.12) ** 2)
                    for charge in (2, 3):
                        for isotope_offset in range(3):
                            mz = (candidate.neutral_mass + isotope_offset * ISOTOPE_MASS_DIFF) / charge + PROTON
                            points.append((mz, peak * (1.0 - isotope_offset * 0.2)))
                points.sort()
                scans.append({
                    "scan_id": f"scan={index}",
                    "rt": rt,
                    "aligned_rt": rt,
                    "mz": [mz for mz, _ in points],
                    "intensity": [intensity for _, intensity in points],
                })
            return scans

        pairs = build_peptide_form_comparisons(payload, psms, annotations, candidates, spectra)
        terminal_pair = next(row for row in pairs if "C-terminal Lys clipping" in row["modification"])
        self.assertEqual(terminal_pair["sequence"], "PEPTIDEK")
        self.assertIn("FG_RETAINED", terminal_pair["linked_feature_ids"])
        self.assertIn("FG_CLIPPED", terminal_pair["linked_feature_ids"])
        self.assertIn("FG_RETAINED", terminal_pair["unmodified_linked_feature_ids"])
        self.assertIn("FG_CLIPPED", terminal_pair["modified_linked_feature_ids"])

        levels = build_modification_level_quantitation(pairs)
        terminal_level = next(row for row in levels if row["event_type"] == "c_terminal_lys_processing")
        self.assertEqual(terminal_level["quantitation_status"], "formal_relative_quantitation")
        retained_form = next(form for form in terminal_level["forms"] if form["label"] == "C端 Lys 保留")
        clipped_form = next(form for form in terminal_level["forms"] if form["label"] == "C端 Lys 已剪切")
        self.assertIn("FG_RETAINED", retained_form["linked_feature_ids"])
        self.assertIn("FG_CLIPPED", clipped_form["linked_feature_ids"])
        self.assertIsNotNone(retained_form["best_psm"])
        self.assertIsNotNone(clipped_form["best_psm"])
        self.assertGreater(float(retained_form["neutral_mass"]), 0.0)
        self.assertAlmostEqual(retained_form["relative_level_by_sample"]["reference"], 0.2, places=2)
        self.assertAlmostEqual(retained_form["relative_level_by_sample"]["test"], 0.5, places=2)
        self.assertTrue(terminal_level["high_value_candidate"])

    def test_nested_peptides_create_site_consensus_and_digestion_warning(self) -> None:
        def row(
            quantitation_id: str,
            sequence: str,
            end: int,
            reference_total: float,
            test_total: float,
            reference_level: float,
            test_level: float,
        ) -> dict[str, object]:
            form_id = f"modified:{sequence}:PyroGlu-Q@1"
            delta = test_level - reference_level
            return {
                "event_type": "site_modification_distribution",
                "quantifiable": True,
                "quantitation_id": quantitation_id,
                "rank": 1,
                "chain": "HC",
                "start": 1,
                "end": end,
                "sequence": sequence,
                "reference_sample": "QL2519",
                "quantitation_status": "formal_relative_quantitation",
                "total_area_by_sample": {"QL2519": reference_total, "Reference": test_total},
                "forms": [{
                    "form_id": form_id,
                    "label": "PyroGlu-Q@1",
                    "modification": "PyroGlu-Q@1",
                    "is_unmodified": False,
                    "included_in_denominator": True,
                    "ms2_confidence": "B_high_confidence_inferred",
                    "status": "ms2_confirmed",
                }],
                "pairwise_comparisons": [{
                    "reference_sample": "QL2519",
                    "test_sample": "Reference",
                    "form_changes": [{
                        "form_id": form_id,
                        "reference_level": reference_level,
                        "test_level": test_level,
                        "percentage_point_difference": delta,
                        "direction": "increased" if delta > 0.005 else "stable",
                    }],
                }],
            }

        summaries = build_proteolytic_modification_families([
            row("MODLEVEL_0002", "QVQLVQSGAEVK", 12, 124.91, 92.65, 0.8498, 0.9815),
            row("MODLEVEL_0009", "QVQLVQSGAEVKKPGASVK", 19, 6.06, 23.53, 0.7727, 0.9557),
        ])
        self.assertEqual(len(summaries), 1)
        summary = summaries[0]
        self.assertTrue(summary["high_value_consensus"])
        site = summary["site_conclusions"][0]
        self.assertEqual(site["event_label"], "PyroGlu-Q@1")
        self.assertTrue(site["pairwise_comparisons"][0]["consistent"])
        digestion = summary["digestion_comparisons"][0]
        self.assertTrue(digestion["suspected_digestion_redistribution"])
        self.assertAlmostEqual(
            digestion["uncorrected_combined_test_over_reference_fold"],
            (92.65 + 23.53) / (124.91 + 6.06),
        )
        glyco_rows = [
            row("MODLEVEL_0010", "QVQLVQSGAEVK", 12, 10.0, 8.0, 0.4, 0.6),
            row("MODLEVEL_0011", "QVQLVQSGAEVKKPGASVK", 19, 2.0, 6.0, 0.3, 0.7),
        ]
        for glyco_row in glyco_rows:
            glyco_row["event_type"] = "glycoform_distribution"
        self.assertEqual(len(build_proteolytic_modification_families(glyco_rows)), 1)


if __name__ == "__main__":
    unittest.main()
