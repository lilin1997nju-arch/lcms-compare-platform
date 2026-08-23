"""Build the self-contained Windows STDIO MCP executable with PyInstaller."""

from __future__ import annotations

import shutil
import subprocess
import sys
import os
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[1]
OUTPUT = WORKSPACE / "release" / "mcp"
WORK = WORKSPACE / "release" / "mcp-build"
BUILD_TOOLS = WORKSPACE / "release" / "build-tools"


def main() -> None:
    if BUILD_TOOLS.exists():
        sys.path.insert(0, str(BUILD_TOOLS))
    try:
        import PyInstaller  # noqa: F401
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "PyInstaller is required: python -m pip install --target release/build-tools -r desktop/requirements-build.txt"
        ) from exc

    if WORK.exists():
        shutil.rmtree(WORK)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--noconfirm",
        "--clean",
        "--onefile",
        "--console",
        "--name",
        "LCMS-MCP",
        "--distpath",
        str(OUTPUT),
        "--workpath",
        str(WORK / "work"),
        "--specpath",
        str(WORK),
        "--paths",
        str(WORKSPACE / "lcms_department_platform"),
        "--paths",
        str(WORKSPACE / "lcms_feature_mvp"),
        str(WORKSPACE / "lcms_department_platform" / "mcp_stdio.py"),
    ]
    environment = os.environ.copy()
    existing_python_path = environment.get("PYTHONPATH", "")
    environment["PYTHONPATH"] = os.pathsep.join(part for part in (str(BUILD_TOOLS), existing_python_path) if part)
    subprocess.run(command, cwd=WORKSPACE, env=environment, check=True)


if __name__ == "__main__":
    main()
