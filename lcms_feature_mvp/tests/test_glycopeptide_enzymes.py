import ast
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import lcms_msms as msms
from core.lcms_enzymes import ENZYMES
from core.lcms_models import LCMSSpectrumScan
import run_msms_compare


class GlycopeptideEnzymeTests(unittest.TestCase):
    sequences = {
        "trypsin": "ACNYSTK", "trypsin_p": "ACNYSTK",
        "lys_c": "ACNYSTK", "lys_c_p": "ACNYSTK",
        "arg_c": "ACNYSTR", "glu_c": "ACNYSTE", "glu_c_de": "ACNYSTD",
        "asp_n": "DACNYSTQ", "lys_n": "KACNYSTQ", "chymotrypsin": "ACNYSTF",
    }

    def test_all_enzymes_preserve_terminal_and_modification_mass(self):
        self.assertEqual(set(self.sequences), set(ENZYMES))
        for enzyme, sequence in self.sequences.items():
            for fixed_cam in (False, True):
                with self.subTest(enzyme=enzyme, fixed_cam=fixed_cam):
                    backbone = msms._candidate("p", 1, len(sequence), sequence, carbamidomethyl_cys=fixed_cam)
                    candidates = msms._targeted_n_glycan_candidates([backbone], enzyme=enzyme)
                    targets = [c for c in candidates if not c.is_decoy]
                    decoys = [c for c in candidates if c.is_decoy]
                    self.assertGreater(len(targets), 0)
                    self.assertEqual(len(targets), len(decoys))
                    for target in targets:
                        expected = msms._decoy(target, enzyme)
                        self.assertIn(expected, decoys)
                        self.assertAlmostEqual(expected.neutral_mass, target.neutral_mass)
                        self.assertEqual(expected.carbamidomethyl_cys, fixed_cam)
                        endpoint = 0 if ENZYMES[enzyme].terminus == "N" else -1
                        self.assertEqual(expected.sequence[endpoint], sequence[endpoint])
                        for position, _, _ in expected.modifications:
                            self.assertEqual(expected.sequence[position], "N")
                        if enzyme in {"asp_n", "lys_n", "glu_c", "glu_c_de", "chymotrypsin"}:
                            self.assertNotEqual(expected.sequence, msms._decoy(target).sequence)

    def test_public_search_forwards_enzyme_to_glycan_generator(self):
        with patch.object(msms, "_targeted_n_glycan_candidates", wraps=msms._targeted_n_glycan_candidates) as generate:
            self.assertEqual(msms.search_feature_glycopeptide_scans({}, [], [], enzyme="lys_n"), [])
            generate.assert_called_once_with([], enzyme="lys_n")
        with self.assertRaises(ValueError):
            msms.search_feature_glycopeptide_scans({}, [], [], enzyme="invalid")

    def test_non_trypsin_glycopeptide_spectra_are_identified(self):
        for enzyme in ("asp_n", "lys_n", "glu_c"):
            with self.subTest(enzyme=enzyme):
                sequence = self.sequences[enzyme]
                backbones = msms.generate_candidates({"p": sequence}, enzyme=enzyme,
                            max_missed_cleavages=0, max_variable_modifications=0)
                target = next(c for c in backbones if not c.is_decoy and c.sequence == sequence)
                glycan = next(g for g in msms.TARGETED_N_GLYCANS if g["name"] == "G0F N-glycan")
                mass = target.neutral_mass + glycan["mass"]
                hexnac = 203.0793725330
                pairs = [(138.0550,2200.), (204.086649,6000.), (366.139472,4200.)]
                pairs.extend((target.neutral_mass+n*hexnac+msms.PROTON, 2600.) for n in (0,1,2))
                for label, mz, series, ordinal in msms.theoretical_fragments(target):
                    if "^" in label:
                        continue
                    pairs.append((mz,1500.))
                    fragment = sequence[:ordinal] if series == "b" else sequence[-ordinal:]
                    if "N" in fragment:
                        pairs.append((mz+hexnac,1700.))
                pairs.sort()
                scan = LCMSSpectrumScan("scan=glyco", "raw", "sample", 10., 2,
                            [p[0] for p in pairs], [p[1] for p in pairs], sum(p[1] for p in pairs),
                            204.086649, 6000., precursor_mz=mass/2+msms.PROTON, precursor_charge=2,
                            isolation_window_lower_offset=.8, isolation_window_upper_offset=.8,
                            activation_method="HCD")
                payload = {"global_feature_groups": [{"feature_group_id":"FG_GLYCO",
                            "representative_mz":scan.precursor_mz, "true_peak_mz":scan.precursor_mz,
                            "representative_rt":10., "difference_type":"area_changed", "ranking_score":2.}]}
                matches = msms.search_feature_glycopeptide_scans(payload,[scan],backbones,enzyme=enzyme)
                self.assertEqual(len(matches),1)
                self.assertEqual(matches[0]["sequence"],sequence)
                self.assertEqual(matches[0]["glycan_name"],"G0F N-glycan")

    def test_runner_passes_selected_enzyme_not_default(self):
        tree = ast.parse(Path(run_msms_compare.__file__).read_text(encoding="utf-8"))
        calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
                 and isinstance(node.func, ast.Name) and node.func.id == "search_feature_glycopeptide_scans"]
        self.assertEqual(len(calls), 1)
        selected = next(keyword.value for keyword in calls[0].keywords if keyword.arg == "enzyme")
        self.assertEqual(ast.unparse(selected), "args.enzyme")

    def test_legacy_trypsin_default_is_unchanged(self):
        sequence = self.sequences["trypsin"]
        backbone = msms._candidate("p", 1, len(sequence), sequence)
        self.assertEqual(msms._targeted_n_glycan_candidates([backbone]),
                         msms._targeted_n_glycan_candidates([backbone], enzyme="trypsin"))


if __name__ == "__main__":
    unittest.main()
