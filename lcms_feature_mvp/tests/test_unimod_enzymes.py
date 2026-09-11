from copy import deepcopy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.lcms_enzymes import digest_enzyme, validate_enzyme
from core.lcms_msms import _candidate, digest_trypsin, generate_candidates, theoretical_fragments, PROTON
from core.lcms_models import LCMSSpectrumScan
from core.lcms_unimod import applicable_sites, load_catalog, search_unimod_rescue


class EnzymeTests(unittest.TestCase):
    def peptides(self, seq, enzyme, missed=0):
        return [p[3] for p in digest_enzyme({"p": seq}, missed, 1, 60, enzyme)]

    def test_default_and_proline(self):
        self.assertEqual(self.peptides("AKPRK", "trypsin"), ["AKPR", "K"])
        self.assertEqual(self.peptides("AKPRK", "trypsin_p"), ["AK", "PR", "K"])
        self.assertEqual(digest_trypsin({"p": "AKPRKPEPTIDER"}, 2, 1, 60),
                         digest_enzyme({"p": "AKPRKPEPTIDER"}, 2, 1, 60))
        self.assertEqual(validate_enzyme(None), "trypsin")
        with self.assertRaises(ValueError):
            validate_enzyme("unknown")

    def test_n_terminal_and_gluc_variants(self):
        self.assertEqual(self.peptides("ADAD", "asp_n"), ["A", "DA", "D"])
        self.assertEqual(self.peptides("AKAK", "lys_n"), ["A", "KA", "K"])
        self.assertEqual(self.peptides("ADEAK", "glu_c"), ["ADE", "AK"])
        self.assertEqual(self.peptides("ADEAK", "glu_c_de"), ["AD", "E", "AK"])
        self.assertEqual(self.peptides("AKEPK", "lys_c"), ["AK", "EPK"])
        self.assertIn("ADEAK", self.peptides("ADEAK", "glu_c", 1))

    def test_candidates_and_decoys_use_selected_enzyme(self):
        candidates = generate_candidates({"p": "KPEPTIDEKAACDE"}, min_length=3,
                                         max_missed_cleavages=0, enzyme="lys_n",
                                         max_variable_modifications=0)
        targets = [p for p in candidates if not p.is_decoy]
        decoys = [p for p in candidates if p.is_decoy]
        self.assertEqual(len(targets), len(decoys))
        for target, decoy in zip(targets, decoys):
            self.assertTrue(target.sequence.startswith("K"))
            self.assertTrue(decoy.sequence.startswith("K"))
            self.assertAlmostEqual(target.neutral_mass, decoy.neutral_mass)
            self.assertEqual(target.proteolysis, "fully_enzymatic:lys_n")


class UnimodTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = load_catalog()

    def record(self, site="S", position="Anywhere", classification="Post-translational"):
        return {"rules": [{"site": site, "position": position, "classification": classification}]}

    def test_full_snapshot_mass_sign_and_precision(self):
        self.assertGreater(len(self.catalog.records), 1500)
        self.assertIn("Phospho", [r["title"] for r in self.catalog.match_mass(79.966331, 0.001)])
        self.assertIn("Amidated", [r["title"] for r in self.catalog.match_mass(-0.984016, 0.001)])
        self.assertEqual(self.catalog.match_mass(float("nan"), 0.01), [])
        self.assertEqual(self.catalog.match_mass(79.970, 0.0001), [])

    def test_residue_peptide_and_protein_terminal_rules(self):
        internal = _candidate("p", 3, 8, "SACDES")
        self.assertEqual([i for i, _ in applicable_sites(self.record(), internal)], [0, 5])
        self.assertEqual([i for i, _ in applicable_sites(self.record("S", "Any N-term"), internal)], [0])
        self.assertEqual(applicable_sites(self.record("N-term", "Protein N-term"), internal), [])
        self.assertEqual(applicable_sites(self.record("C-term", "Protein C-term"), internal), [])
        self.assertEqual(len(applicable_sites(self.record("C-term", "Protein C-term"), internal, 8)), 1)
        self.assertEqual(applicable_sites(self.record("C-term", "Protein C-term"), internal, 9), [])
        self.assertEqual(applicable_sites(self.record("C"), internal), [])
        self.assertEqual(applicable_sites(self.record(classification="AA substitution"), internal), [])
        self.assertEqual(applicable_sites(self.record(position="Unrecognized"), internal), [])
        self.assertEqual(applicable_sites({"title": "Xlink:DSS[156]", **self.record()}, internal), [])
        self.assertEqual(applicable_sites({"title": "b-type-ion", **self.record("C-term")}, internal), [])
        self.assertEqual(applicable_sites({"title": "Lys-loss", **self.record("K")}, _candidate("p", 1, 4, "AKAK")), [])

    def example(self, relation="selected_precursor"):
        sequence = "ASPEPTIDEK"
        backbone = _candidate("p", 1, len(sequence), sequence)
        modified = _candidate("p", 1, len(sequence), sequence, ((1, "Phospho", 79.966331),))
        peaks = sorted(mz for _, mz, _, _ in theoretical_fragments(modified))
        scan = LCMSSpectrumScan("scan=1", "raw", "sample-A", 1.0, 2,
                               peaks, [100.0]*len(peaks), 100.0*len(peaks), peaks[0], 100.0,
                               precursor_mz=modified.neutral_mass/2+PROTON, precursor_charge=2)
        feature = {"feature_group_id": "FG1", "difference_type": "area_changed",
                   "confidence": None, "sequence": None, "ms2_status": "selected_precursor_unidentified",
                   "coverage_scan_by_sample": {"sample-A": {"sample_id": "sample-A", "scan_id": "scan=1",
                       "relation_type": relation, "precursor_charge": 2}}}
        return feature, scan, backbone

    def test_second_pass_scores_new_modification_without_changing_identification(self):
        feature, scan, backbone = self.example()
        original = deepcopy(feature)
        summary = search_unimod_rescue([feature], [scan], [backbone], {"p": backbone.sequence})
        self.assertEqual(summary["features_with_candidates"], 1)
        for key, value in original.items():
            self.assertEqual(feature[key], value)
        hits = feature["unimod_rescue_candidates"]
        phospho = next(r for r in hits if r["unimod_title"] == "Phospho")
        self.assertIsNone(phospho["q_value"])
        self.assertTrue(phospho["exploratory_only"])
        self.assertEqual(phospho["sample_id"], "sample-A")
        self.assertEqual(phospho["display_site"], 2)
        self.assertEqual(phospho["localization_status"], "not_validated")
        self.assertGreater(phospho["fragment_coverage"], 0.9)
        self.assertTrue(any(p["label"].startswith("b") for p in phospho["spectrum_peaks"]))

    def test_no_rescue_for_coisolation_or_existing_formal_identification(self):
        for confidence, relation in [(None, "isolation_window_only"), ("B_high_confidence_inferred", "selected_precursor"),
                                     ("C_tentative", "selected_precursor")]:
            feature, scan, backbone = self.example(relation)
            feature["confidence"] = confidence
            summary = search_unimod_rescue([feature], [scan], [backbone], {"p": backbone.sequence})
            self.assertEqual(summary["scored_site_candidates"], 0)

    def test_noise_and_unknown_charge_do_not_produce_identifications(self):
        feature, scan, backbone = self.example()
        scan.mz_array = [10.1, 11.2, 12.3]
        scan.intensity_array = [1, 1, 1]
        summary = search_unimod_rescue([feature], [scan], [backbone], {"p": backbone.sequence})
        self.assertEqual(summary["reported_candidates"], 0)
        feature["coverage_scan_by_sample"]["sample-A"]["precursor_charge"] = None
        summary = search_unimod_rescue([feature], [scan], [backbone], {"p": backbone.sequence})
        self.assertEqual(summary["selected_scans"], 0)


if __name__ == "__main__":
    unittest.main()
