from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
import zipfile
from pathlib import Path


EXCLUDE_DIRS = {
    ".git",
    ".venv",
    "venv",
    "__pycache__",
    "outputs",
    "jobs",
    "logs",
    "state",
    ".local-tools",
    "handoff_packages",
}
EXCLUDE_SUFFIXES = {".raw", ".mzml", ".sqlite", ".db", ".pyc"}
EXCLUDE_NAMES = {
    "python-bg-err.txt",
    "python-bg-out.txt",
    "LC-MSMS平台_MS2功能展示.pptx.inspect.ndjson",
}
INCLUDE_DIRS = {"data", "lcms_feature_mvp", "lcms_realdata_platform", "lcms_department_platform", "docs", "examples", "tools"}
INCLUDE_FILES = {
    "START_HERE.md",
    "VERSION.json",
    "requirements.txt",
    ".gitignore",
    "LC-MS-MS差异组分鉴定推测功能设计.md",
    "LC-MSMS平台_MS2功能展示.pptx",
    "一键启动LC-MSMS平台.cmd",
    "一键启动部门平台.cmd",
    "停止部门平台.cmd",
    "start_department_platform.cmd",
    "stop_department_platform.cmd",
    "build_migration_package.cmd",
    "打包部门迁移包.cmd",
}


def should_copy(path: Path, root: Path) -> bool:
    relative = path.relative_to(root)
    if any(part in EXCLUDE_DIRS for part in relative.parts):
        return False
    if any(part == "outputs" or part.startswith("outputs_") for part in relative.parts):
        return False
    if path.name in EXCLUDE_NAMES or path.suffix.lower() in EXCLUDE_SUFFIXES:
        return False
    if path.name.startswith("lcms_department_platform_") and path.suffix == ".pid":
        return False
    return True


def copy_tree(root: Path, staging: Path) -> list[str]:
    copied: list[str] = []
    for top in sorted(root.iterdir(), key=lambda item: item.name.lower()):
        if top.is_dir() and top.name in INCLUDE_DIRS:
            for source in top.rglob("*"):
                if source.is_file() and should_copy(source, root):
                    target = staging / source.relative_to(root)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, target)
                    copied.append(str(target.relative_to(staging)))
        elif top.is_file() and top.name in INCLUDE_FILES and should_copy(top, root):
            target = staging / top.name
            shutil.copy2(top, target)
            copied.append(top.name)
    return copied


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a code-and-runtime migration zip for the LC-MS platform.")
    parser.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--output", default="")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    output = Path(args.output).resolve() if args.output else root / "handoff_packages" / f"LC-MSMS_department_migration_{__import__('time').strftime('%Y%m%d')}_with_ThermoRawFileParser.zip"
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="lcms_migration_") as temp:
        staging = Path(temp)
        copied = copy_tree(root, staging)
        readme = staging / "MIGRATION_README.md"
        readme.write_text(
            """# LC-MS 部门迁移包

这是可复制到部门其他 Windows 电脑的代码与运行包。

## 启动

1. 安装 Python 3.10 或更高版本。
2. 迁移包已包含 ThermoRawFileParser 的完整自包含运行目录，可直接处理 Thermo `.raw`；如需替换为其他授权版本，可在启动脚本中传入 `-ParserPath`。仅使用 `.mzML` 时不需要转换器。
3. 双击 `start_department_platform.cmd`（中文名称入口也一并提供）。
4. 浏览器打开 `http://127.0.0.1:8770/`。

## 验证

运行 `powershell -ExecutionPolicy Bypass -File .\\tools\\verify_migration.ps1`。

## 数据与结果

迁移包刻意不包含 RAW/mzML、任务 jobs、SQLite 结果和运行日志。到新电脑后，从页面上传数据；任务结果会写入 `lcms_department_platform\\jobs\\<task_id>`。
该包不包含权限管理，适合部门内网或单机使用。
""",
            encoding="utf-8",
        )
        manifest = {
            "package_type": "lcms_department_migration",
            "source_root": str(root),
            "files": sorted(copied + ["MIGRATION_README.md", "MIGRATION_MANIFEST.json"]),
            "excluded_runtime_data": ["RAW", "mzML", "SQLite results", "jobs", "logs", "state", "cache"],
        }
        (staging / "MIGRATION_MANIFEST.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
        if output.exists():
            output.unlink()
        with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for source in sorted(staging.rglob("*")):
                if source.is_file():
                    archive.write(source, source.relative_to(staging).as_posix())
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(output.suffix + ".sha256").write_text(f"{digest}  {output.name}\n", encoding="ascii")
    print(json.dumps({"output": str(output), "bytes": output.stat().st_size, "sha256": digest}, ensure_ascii=False))


if __name__ == "__main__":
    main()
