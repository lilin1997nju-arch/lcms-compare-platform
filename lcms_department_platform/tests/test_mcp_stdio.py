import io
import json
import tempfile
import unittest
from pathlib import Path

from lcms_department_platform.mcp_server import MCPApplication
from lcms_department_platform.mcp_stdio import MCPConfig, ReadOnlyTaskStore, serve_stdio


class MCPStdioTests(unittest.TestCase):
    def test_initialize_and_tool_list_over_json_lines(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            state = root / "state"
            state.mkdir()
            (state / "tasks.json").write_text("[]", encoding="utf-8")
            config = MCPConfig(root)
            app = MCPApplication(config, ReadOnlyTaskStore(config))
            requests = b"\n".join(
                [
                    json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": "2025-06-18"}}).encode(),
                    json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized", "params": {}}).encode(),
                    json.dumps({"jsonrpc": "2.0", "id": 2, "method": "tools/list", "params": {}}).encode(),
                ]
            ) + b"\n"
            output = io.BytesIO()
            serve_stdio(app, io.BytesIO(requests), output)
            responses = [json.loads(line) for line in output.getvalue().splitlines()]
            self.assertEqual([item["id"] for item in responses], [1, 2])
            self.assertEqual(responses[0]["result"]["serverInfo"]["name"], "lcms-department-platform")
            tools = responses[1]["result"]["tools"]
            self.assertTrue(next(item for item in tools if item["name"] == "list_tasks")["annotations"]["readOnlyHint"])
            self.assertFalse(next(item for item in tools if item["name"] == "save_agent_annotation")["annotations"]["readOnlyHint"])

    def test_parse_error_does_not_stop_server(self):
        config = MCPConfig(Path("unused"))
        app = MCPApplication(config, ReadOnlyTaskStore(config))
        output = io.BytesIO()
        serve_stdio(app, io.BytesIO(b"not-json\n{\"jsonrpc\":\"2.0\",\"id\":3,\"method\":\"ping\"}\n"), output)
        responses = [json.loads(line) for line in output.getvalue().splitlines()]
        self.assertEqual(responses[0]["error"]["code"], -32700)
        self.assertEqual(responses[1]["id"], 3)


if __name__ == "__main__":
    unittest.main()
