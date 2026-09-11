import re
import shutil
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_peak_first_compare import PEAK_FIRST_TEMPLATE
from run_msms_compare import parse_args


class UnimodUITests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("node"), "Node is needed to validate report JavaScript")
    def test_embedded_javascript_syntax(self):
        for script in re.findall(r"<script[^>]*>(.*?)</script>", PEAK_FIRST_TEMPLATE, re.S):
            completed = subprocess.run(["node", "--check", "-"], input=script, encoding="utf-8", capture_output=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)

    @unittest.skipUnless(shutil.which("node"), "Node is needed to validate report JavaScript")
    def test_unimod_sample_filter_never_falls_back_to_other_sample(self):
        start = PEAK_FIRST_TEMPLATE.index("function unimodEvidenceCandidates(")
        end = PEAK_FIRST_TEMPLATE.index("function renderUnimodEvidenceSelector", start)
        script = PEAK_FIRST_TEMPLATE[start:end] + '''
const assert = require('node:assert/strict');
const feature={unimod_rescue_candidates:[{sample_id:'A'},{sample_id:'B'}]};
assert.equal(unimodEvidenceCandidates(feature).length,2);
assert.deepEqual(unimodEvidenceCandidates(feature,'B'),[{sample_id:'B'}]);
assert.deepEqual(unimodEvidenceCandidates(feature,'missing'),[]);
assert.deepEqual(unimodEvidenceCandidates({}),[]);
'''
        completed = subprocess.run(["node", "-e", script], capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_cli_enzyme_and_controlled_comparison_flag(self):
        with patch.object(sys, "argv", ["run", "--mzml", "a.mzML", "--fasta", "a.fasta", "--output-dir", "out",
                                       "--enzyme", "glu_c_de", "--disable-unimod-rescue"]):
            args = parse_args()
        self.assertEqual(args.enzyme, "glu_c_de")
        self.assertTrue(args.disable_unimod_rescue)


if __name__ == "__main__":
    unittest.main()
