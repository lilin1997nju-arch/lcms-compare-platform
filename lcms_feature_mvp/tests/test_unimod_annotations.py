import re
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.lcms_unimod_annotations import COMMON_NAMES, modification_explanation, unimod_name_notes
from core.lcms_unimod import load_catalog
from run_peak_first_compare import PEAK_FIRST_TEMPLATE


class UnimodAnnotationTests(unittest.TestCase):
    def test_common_names_have_no_note(self):
        for name in COMMON_NAMES:
            self.assertEqual(modification_explanation({"id": 0, "title": name}), "")
        for record in load_catalog().records:
            if record["title"] in COMMON_NAMES:
                self.assertNotIn(str(record["id"]), unimod_name_notes())

    def test_obscure_notes_include_original_definition_and_caveat(self):
        records = {str(row["id"]): row for row in load_catalog().records}
        for record_id, keyword in [("2119", "VG"), ("53", "HNE"), ("1913", "AGE"), ("1301", "转肽")]:
            note = unimod_name_notes()[record_id]
            self.assertIn(keyword, note)
            self.assertIn(records[record_id]["full_name"], note)
            self.assertIn("不代表本样品已确认", note)

    def test_special_families_and_unknown_names(self):
        for title, keyword in [("Hex(3)HexNAc(2)", "糖组成"), ("Cation:Fe[III]", "价态")]:
            self.assertIn(keyword, modification_explanation({"id": -1, "title": title}))
        self.assertEqual(modification_explanation({"id": -1, "title": "Unrecognized"}), "")
        self.assertEqual(modification_explanation({"id": -1, "title": "Hex(3)unknown"}), "")

    @unittest.skipUnless(shutil.which("node"), "Node is needed for report rendering tests")
    def test_html_escaping_old_records_and_selector_reset(self):
        start = PEAK_FIRST_TEMPLATE.index("const UNIMOD_NAME_NOTES=")
        end = PEAK_FIRST_TEMPLATE.index("function drawFeatureMs2", start)
        # Use the actual report escape helper, including quote escaping for title attributes.
        helper = re.search(r"function escapeHtml\([^\n]+", PEAK_FIRST_TEMPLATE).group(0)
        script = helper + "\n" + PEAK_FIRST_TEMPLATE[start:end] + r'''
const assert=require('node:assert/strict');
assert.equal(unimodModificationHtml({unimod_id:35},'Oxidation'),'Oxidation');
assert.equal(unimodModificationHtml({unimod_id:999999},'<unknown>'),'&lt;unknown&gt;');
assert.equal(unimodNameNote({unimod_id:'__proto__'}),'');
const html=unimodModificationHtml({unimod_id:2119},'ValGly <test> "quote"');
assert.ok(html.includes('class="unimod-explained"'));
assert.ok(html.includes('title="'));
assert.ok(html.includes('UNIMOD:2119'));
assert.ok(html.includes('&lt;test&gt; &quot;quote&quot;'));
assert.ok(!html.includes('<test>'));
// No new result fields are required: old unimod_id-only records work.
assert.equal(unimodNameNote({unimod_id:'2119'}),unimodNameNote({unimod_id:2119}));
const select={dataset:{context:'F|'},value:'0',options:[{value:''},{value:'0'},{value:'1'}],removeAttribute(key){delete this[key];}};
const note={};
const $=id=>id==='featureMs2Evidence'?select:note;
const state={msms:{}};
const featureMs2SampleDisplayName=id=>id||'';
const nice=value=>String(value||0);
const feature={feature_group_id:'F',unimod_rescue_candidates:[{unimod_id:2119},{unimod_id:35}]};
renderUnimodEvidenceSelector(feature);
assert.ok(select.title.includes('UNIMOD:2119'));
select.value='1';
renderUnimodEvidenceSelector(feature);
assert.equal(select.title,undefined);
select.value='0';
renderUnimodEvidenceSelector(feature);
renderUnimodEvidenceSelector({...feature,feature_group_id:'other'});
assert.equal(select.title,undefined);
assert.equal(select.value,'');
'''
        result = subprocess.run(["node", "-"], input=script, encoding="utf-8", capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
