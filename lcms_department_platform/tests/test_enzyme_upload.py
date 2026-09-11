import json
from http.server import ThreadingHTTPServer
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from urllib.request import Request, urlopen
from urllib.error import HTTPError

from lcms_department_platform.server import PortalConfig, TaskStore, make_handler
from core.lcms_sample_prep import ALL_PREP_FIELDS, normalize_sample_prep, sample_prep_cli_args


class EnzymeUploadTests(unittest.TestCase):
    def test_upload_persists_enzyme_and_preserves_fixed_parameters(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = PortalConfig("127.0.0.1", 0, Path(temporary), Path("unused-parser"), sys.executable)
            store = TaskStore(config)
            server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(config, store))
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            def submit(enzyme, prep=None):
                boundary = "lcms-enzyme-test-boundary"
                parts = []
                fields = {"task_name": "enzyme-test", "sample_names": '["A","B"]'}
                fields.update(prep or {})
                if enzyme is not None:
                    fields["enzyme"] = enzyme
                for name, value in fields.items():
                    parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="{name}"\r\n\r\n{value}\r\n')
                for name in ("a.mzML", "b.mzML"):
                    parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="files"; filename="{name}"\r\nContent-Type: application/xml\r\n\r\n<mzML/>\r\n')
                body = ("".join(parts)+f"--{boundary}--\r\n").encode()
                request = Request(f"http://127.0.0.1:{server.server_port}/api/tasks", data=body,
                                  headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
                with urlopen(request, timeout=10) as response:
                    return json.load(response)
            try:
                with urlopen(f"http://127.0.0.1:{server.server_port}/", timeout=10) as response:
                    page = response.read().decode("utf-8")
                self.assertNotIn("__PREP_OPTIONS__", page)
                for key in ALL_PREP_FIELDS:
                    self.assertIn(f'name="prep_{key}"', page)
                submit("glu_c_de", {"prep_biotin": "no", "prep_tmt": "yes", "prep_itraq": "no", "prep_reduction": "reduced", "prep_alkylation": "none"})
                submit(None)
                tasks = store.load()
                self.assertEqual({task["params"]["enzyme"] for task in tasks}, {"glu_c_de", "trypsin"})
                for task in tasks:
                    prep = task["params"]["sample_prep"]
                    expected = normalize_sample_prep({"biotin": "no", "tmt": "yes", "itraq": "no", "reduction": "reduced", "alkylation": "none"}) if task["params"]["enzyme"] == "glu_c_de" else normalize_sample_prep()
                    self.assertEqual(prep, expected)
                    self.assertEqual(task["params"]["top_n_mz"], 200)
                    self.assertEqual(task["params"]["candidate_max_per_tic"], 200)
                    self.assertEqual(task["params"]["top_n_peaks"], 0)
                with self.assertRaises(HTTPError) as caught:
                    submit("not-an-enzyme")
                self.assertEqual(caught.exception.code, 400)
                self.assertEqual(len(store.load()), 2)
                with self.assertRaises(HTTPError) as caught:
                    submit(None, {"prep_alkylation": "nem"})
                self.assertEqual(caught.exception.code, 400)
                self.assertEqual(len(store.load()), 2)
                with self.assertRaises(HTTPError) as caught:
                    submit(None, {"prep_biotin": "false"})
                self.assertEqual(caught.exception.code, 400)
                self.assertEqual(len(store.load()), 2)
            finally:
                server.shutdown()
                server.server_close()
                thread.join(timeout=2)

    def test_worker_cli_preserves_unknown_for_old_tasks(self):
        self.assertEqual(sample_prep_cli_args(), [part for key in ALL_PREP_FIELDS for part in (f"--prep-{key}", "unknown")])
        command = sample_prep_cli_args({"biotin": "no", "alkylation": "none"})
        self.assertEqual(command[command.index("--prep-biotin")+1], "no")
        self.assertEqual(command[command.index("--prep-alkylation")+1], "none")


if __name__ == "__main__":
    unittest.main()
