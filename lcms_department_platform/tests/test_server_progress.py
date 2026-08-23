from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


from lcms_department_platform.server import TaskWorker


class ProgressStore:
    def __init__(self, jobs_dir: Path) -> None:
        self.config = SimpleNamespace(jobs_dir=jobs_dir)
        self.updates: list[dict[str, object]] = []

    def get(self, task_id: str) -> dict[str, object]:
        return {"task_id": task_id, "status": "running", "cancel_requested": False}

    def update(self, task_id: str, **updates: object) -> dict[str, object]:
        self.updates.append(dict(updates))
        return {"task_id": task_id, **updates}


class ServerProgressTests(unittest.TestCase):
    def test_command_progress_markers_update_overall_range(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            store = ProgressStore(Path(tmp) / "jobs")
            worker = TaskWorker(store)  # type: ignore[arg-type]
            worker.run_command(
                "task_progress",
                [
                    sys.executable,
                    "-c",
                    "print('LCMS_PROGRESS\\t0.5\\tmiddle_stage', flush=True)",
                ],
                Path.cwd(),
                progress_range=(20, 80),
                progress_stage="测试阶段",
            )

        progress_updates = [update for update in store.updates if "progress" in update]
        self.assertIn(50, [int(update["progress"]) for update in progress_updates])
        self.assertEqual(progress_updates[-1]["progress"], 80)
        self.assertIn("中间阶段", [str(update.get("stage")) for update in store.updates])


if __name__ == "__main__":
    unittest.main()
