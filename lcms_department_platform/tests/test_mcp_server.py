import json
import sqlite3
import tempfile
import unittest
from contextlib import closing
from pathlib import Path
from types import SimpleNamespace

from lcms_department_platform.mcp_server import MCPApplication


class FakeStore:
    def __init__(self, task):
        self.task = task

    def load(self):
        return [self.task]

    def get(self, task_id):
        return self.task if task_id == self.task["task_id"] else None


class MCPServerTests(unittest.TestCase):
    def test_read_only_tools_and_ms2_policy(self):
        with tempfile.TemporaryDirectory() as temporary:
            jobs_dir = Path(temporary) / "jobs"
            task_id = "task_001"
            db_path = jobs_dir / task_id / "output" / "lcms_peak_first_compare.sqlite"
            db_path.parent.mkdir(parents=True)
            with closing(sqlite3.connect(db_path)) as connection:
                connection.execute(
                    "CREATE TABLE peak_first_artifacts (artifact_key TEXT PRIMARY KEY, payload_json TEXT NOT NULL, updated_at TEXT)"
                )
                bootstrap = {"project_id": task_id, "sample_ids": ["A", "B"], "ms1_component_groups": []}
                msms = {
                    "total_ms2_scans": 2,
                    "total_accepted_psms": 1,
                    "feature_evidence": [
                        {"feature_group_id": "F_MS2", "ms2_status": "identified_direct_precursor", "sequence": "PEPTIDE", "coverage_scan": {"scan_id": "S1", "spectrum_peaks": [{"mz": 100.0, "intensity": 10.0}]}},
                        {"feature_group_id": "F_NO_MS2", "ms2_status": "no_ms2_acquired", "sequence": None},
                    ],
                    "ms1_component_groups": {"component_groups": []},
                }
                connection.executemany(
                    "INSERT INTO peak_first_artifacts VALUES (?, ?, ?)",
                    [("bootstrap", json.dumps(bootstrap), "now"), ("msms_identifications", json.dumps(msms), "now")],
                )
                connection.commit()

            task = {"task_id": task_id, "task_name": "test", "status": "finished", "progress": 100}
            config = SimpleNamespace(jobs_dir=jobs_dir)
            app = MCPApplication(config, FakeStore(task))

            tools = app.dispatch({"jsonrpc": "2.0", "id": 1, "method": "tools/list", "params": {}})
            self.assertEqual(tools["result"]["tools"][-1]["name"], "get_msms_spectrum")

            response = app.dispatch({"jsonrpc": "2.0", "id": 2, "method": "tools/call", "params": {"name": "list_features", "arguments": {"task_id": task_id, "limit": 10}}})
            payload = json.loads(response["result"]["content"][0]["text"])
            self.assertEqual(payload["total_matching"], 1)
            self.assertEqual(payload["features"][0]["feature_group_id"], "F_MS2")

            response = app.dispatch({"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": {"name": "get_feature_evidence", "arguments": {"task_id": task_id, "feature_id": "F_MS2"}}})
            payload = json.loads(response["result"]["content"][0]["text"])
            self.assertTrue(payload["ms2_available"])
            self.assertTrue(payload["analysis_eligible"])

            response = app.dispatch({"jsonrpc": "2.0", "id": 4, "method": "tools/call", "params": {"name": "get_feature_evidence", "arguments": {"task_id": task_id, "feature_id": "F_NO_MS2"}}})
            payload = json.loads(response["result"]["content"][0]["text"])
            self.assertFalse(payload["ms2_available"])
            self.assertFalse(payload["analysis_eligible"])

    def test_task_path_is_not_traversable(self):
        with tempfile.TemporaryDirectory() as temporary:
            task = {"task_id": "safe", "status": "finished"}
            app = MCPApplication(SimpleNamespace(jobs_dir=Path(temporary)), FakeStore(task))
            response = app.dispatch({"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": "get_task_status", "arguments": {"task_id": "../safe"}}})
            self.assertTrue(response["result"]["isError"])


if __name__ == "__main__":
    unittest.main()
