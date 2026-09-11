from copy import deepcopy
from pathlib import Path
import contextlib
import io
import json
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.lcms_sample_prep import ALL_PREP_FIELDS, REAGENT_IDS, normalize_sample_prep, filter_reagent_rules, sample_prep_model, filter_linear_candidates
from core.lcms_unimod import load_catalog, search_unimod_rescue, applicable_sites
from core.lcms_msms import _candidate, theoretical_fragments, PROTON, CARBAMIDOMETHYL, generate_candidates, generate_sequence_inference_candidates, search_scans
from core.lcms_models import LCMSSpectrumScan
from run_msms_compare import parse_args
import run_msms_compare as runner


class SamplePrepTests(unittest.TestCase):
    def test_cli_pipeline_persists_and_uses_mass_model(self):
        # Exercise the real runner/search/report writer; only the input reader
        # is replaced by a known synthetic spectrum, not the search itself.
        sequence = "ACDEFGHICK"
        free = _candidate("p", 1, len(sequence), sequence, carbamidomethyl_cys=False)
        peaks = sorted(mz for _, mz, _, _ in theoretical_fragments(free))
        scan = LCMSSpectrumScan("scan=1", "raw", "A", 1, 2, peaks, [100.]*len(peaks),
                               100.*len(peaks), peaks[0], 100., precursor_mz=free.neutral_mass/2+PROTON, precursor_charge=2)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "test.fasta"
            fasta.write_text(">p\n"+sequence+"\n", encoding="utf-8")
            for reduction, alkylation in (("reduced", "none"), ("reduced", "cam"), ("none", "none")):
                out = root / (reduction+"_"+alkylation)
                argv = ["run", "--mzml", str(root/"a.mzML"), "--fasta", str(fasta), "--output-dir", str(out),
                        "--prep-reduction", reduction, "--prep-alkylation", alkylation, "--ms2-workers", "1"]
                with patch.object(sys, "argv", argv), patch.object(runner, "read_mzml", return_value=(SimpleNamespace(sample_id="A",file_name="a.mzML"), [scan])), contextlib.redirect_stdout(io.StringIO()):
                    runner.main()
                report = json.loads((out/"lcms_msms_identifications.json").read_text(encoding="utf-8"))
                self.assertEqual(report["parameters"]["sample_prep"]["alkylation"], alkylation)
                self.assertEqual(report["parameters"]["fixed_modification"], "None" if alkylation=="none" else "Carbamidomethyl@C")
                matching = [row for row in report["psms"] if row["sequence"] == sequence]
                if reduction == "reduced" and alkylation == "none":
                    self.assertTrue(matching)
                    self.assertFalse(matching[0]["carbamidomethyl_cys"])
                else:
                    self.assertFalse(matching)
                if reduction == "none":
                    self.assertTrue(all("C" not in row["sequence"] for row in report["psms"]))

    def test_validation_and_cli(self):
        self.assertEqual(normalize_sample_prep(), {key: "unknown" for key in ALL_PREP_FIELDS})
        for bad in ({"biotin": False}, {"biotin": "maybe"}, {"unknown_key": "no"}, []):
            with self.assertRaises(ValueError):
                normalize_sample_prep(bad)
        with patch.object(sys, "argv", ["run", "--mzml", "a", "--fasta", "a", "--output-dir", "out",
                                       "--prep-biotin", "no", "--prep-tmt", "yes"]):
            args = parse_args()
        self.assertEqual((args.prep_biotin, args.prep_tmt, args.prep_itraq), ("no", "yes", "unknown"))

    def test_alkylation_changes_precursor_fragments_and_decoys(self):
        chains = {"p": "ACDEFGHICK"}
        generated = {}
        for reagent in ("unknown", "cam", "none"):
            model = sample_prep_model({"alkylation": reagent})
            generated[reagent] = generate_candidates(chains, max_variable_modifications=0,
                                carbamidomethyl_cys=model["carbamidomethyl_cys"])
        self.assertEqual(generated["unknown"], generated["cam"])
        for cam, free in zip(generated["cam"], generated["none"]):
            self.assertAlmostEqual(cam.neutral_mass-free.neutral_mass, 2*CARBAMIDOMETHYL)
            free_frag = {label: mz for label, mz, _, _ in theoretical_fragments(free)}
            for label, mz, series, ordinal in theoretical_fragments(cam):
                sequence = cam.sequence[:ordinal] if series == "b" else cam.sequence[-ordinal:]
                charge = int(label.split("^")[1]) if "^" in label else 1
                self.assertAlmostEqual(mz-free_frag[label], sequence.count("C")*CARBAMIDOMETHYL/charge)
        free = next(c for c in generated["none"] if not c.is_decoy)
        peaks = sorted(mz for _, mz, _, _ in theoretical_fragments(free))
        scan = LCMSSpectrumScan("s", "raw", "A", 1, 2, peaks, [100.]*len(peaks),
                               100.*len(peaks), peaks[0], 100., precursor_mz=free.neutral_mass/2+PROTON, precursor_charge=2)
        self.assertEqual(len(search_scans([scan], generated["none"])), 1)
        self.assertEqual(search_scans([scan], generated["cam"]), [])

    def test_reduction_restricts_cys_in_both_candidate_pools(self):
        chains = {"p": "ACDEFGHIKPEPTIDER"}
        for generate in (generate_candidates, generate_sequence_inference_candidates):
            candidates = generate(chains)
            self.assertTrue(any("C" in c.sequence for c in candidates))
            for reagent in ("unknown", "reduced"):
                self.assertEqual(filter_linear_candidates(candidates, {"reduction": reagent}), candidates)
            filtered = filter_linear_candidates(candidates, {"reduction": "none"})
            self.assertTrue(all("C" not in c.sequence for c in filtered))
        self.assertTrue(sample_prep_model({"reduction": "none"})["warnings"])

    def test_no_alkylation_excludes_cam_rescue_and_legacy_flag(self):
        record = next(r for r in load_catalog().records if r["id"] == 4)
        self.assertEqual(filter_reagent_rules(record, normalize_sample_prep({"alkylation": "none"})), [])
        argv = ["run", "--mzml", "a", "--fasta", "a", "--output-dir", "out", "--no-carbamidomethyl-cys"]
        with patch.object(sys, "argv", argv):
            self.assertEqual(parse_args().sample_prep["alkylation"], "none")
        with patch.object(sys, "argv", argv+["--prep-alkylation", "cam"]), patch("sys.stderr"):
            with self.assertRaises(SystemExit):
                parse_args()

    def test_curated_ids_and_no_collateral_filtering(self):
        records = {r["id"]: r for r in load_catalog().records}
        for ids in REAGENT_IDS.values():
            self.assertTrue(ids <= records.keys())
        no = normalize_sample_prep({key: "no" for key in REAGENT_IDS})
        for key, ids in REAGENT_IDS.items():
            for record_id in ids - {3}:
                self.assertEqual(filter_reagent_rules(records[record_id], no), [], (key, record_id))
        # Common PTMs, unlabelled HNE, glycans and unrelated artefacts survive.
        for record_id in (4, 7, 21, 24, 35, 53, 2119):
            self.assertEqual(filter_reagent_rules(records[record_id], no), records[record_id]["rules"])
        for state in ("unknown", "yes"):
            prep = normalize_sample_prep({key: state for key in REAGENT_IDS})
            for record in records.values():
                self.assertEqual(filter_reagent_rules(record, prep), record["rules"])

    def test_natural_biotin_is_not_artificial_labelling(self):
        record = next(r for r in load_catalog().records if r["id"] == 3)
        filtered = {**record, "rules": filter_reagent_rules(record, normalize_sample_prep({"biotin": "no"}))}
        self.assertEqual([r["classification"] for r in filtered["rules"]], ["Post-translational"])
        backbone = _candidate("p", 1, 10, "ASPEPTIDEK")
        self.assertEqual([i for i, rule in applicable_sites(filtered, backbone)], [9])

    def test_filter_before_scoring_does_not_change_formal_evidence(self):
        backbone = _candidate("p", 1, 10, "ASPEPTIDEK")
        record = next(r for r in load_catalog().records if r["id"] == 289)
        index = applicable_sites(record, backbone)[0][0]
        modified = _candidate("p", 1, 10, backbone.sequence, ((index, record["title"], record["mass"]),))
        peaks = sorted(mz for _, mz, _, _ in theoretical_fragments(modified))
        scan = LCMSSpectrumScan("scan=1", "raw", "A", 1, 2, peaks, [100.]*len(peaks),
                               100.*len(peaks), peaks[0], 100., precursor_mz=modified.neutral_mass/2+PROTON, precursor_charge=2)
        original = {"feature_group_id": "FG", "difference_type": "area_changed", "sequence": None,
                    "coverage_scan_by_sample": {"A": {"sample_id": "A", "scan_id": "scan=1",
                        "precursor_charge": 2, "relation_type": "selected_precursor"}}}
        outputs = []
        for prep in (None, {"biotin": "yes"}, {"biotin": "no"}):
            f = deepcopy(original)
            summary = search_unimod_rescue([f], [scan], [backbone], {"p": backbone.sequence}, sample_prep=prep)
            for key, value in original.items():
                self.assertEqual(f[key], value)
            outputs.append(f["unimod_rescue_candidates"])
        self.assertEqual(outputs[0], outputs[1])
        self.assertIn(289, {r["unimod_id"] for r in outputs[0]})
        self.assertNotIn(289, {r["unimod_id"] for r in outputs[2]})
        self.assertIn(289, summary["prep_excluded_record_ids"])
        self.assertGreater(summary["prep_excluded_rule_count"], 0)


if __name__ == "__main__":
    unittest.main()
