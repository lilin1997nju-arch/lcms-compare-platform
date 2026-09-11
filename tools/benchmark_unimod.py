"""Sequential A/B on a task's cached mzML and copied SQLite; never edits the task.

Reports MS2 workflow elapsed time, not RAW conversion or MS1 computation.
"""
import argparse
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Use a new output directory to preserve previous comparisons")
    args.output.mkdir(parents=True)
    summary = {"task": str(args.task.resolve()), "scope": "MS2-only, existing cached mzML; same MS1 bootstrap",
               "order": ["baseline", "unimod"], "runs": {}}
    reports = {}
    for name, extra in [("baseline", ["--disable-unimod-rescue"]), ("unimod", [])]:
        destination = args.output / name
        destination.mkdir()
        db_path = destination / "lcms_peak_first_compare.sqlite"
        source_path = args.task / "output/lcms_peak_first_compare.sqlite"
        with sqlite3.connect(source_path.as_uri() + "?mode=ro", uri=True) as source, sqlite3.connect(db_path) as target:
            source.backup(target)
        command = [sys.executable, str(ROOT / "lcms_feature_mvp/run_msms_compare.py"),
                   "--fasta", str(args.task / "references/sequence.fasta"),
                   "--output-dir", str(destination), "--feature-sqlite", str(db_path),
                   "--project-id", args.task.name, "--enzyme", "trypsin", *extra]
        for path in sorted((args.task / "mzML").glob("*.mzML")):
            command.extend(["--mzml", str(path)])
        started = time.perf_counter()
        print(f"Starting {name}: {destination}", flush=True)
        with (destination / "worker.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, cwd=ROOT)
        elapsed = time.perf_counter() - started
        report = json.loads((destination / "lcms_msms_identifications.json").read_text(encoding="utf-8"))
        reports[name] = report
        timings = {}
        for line in (destination / "worker.log").read_text(encoding="utf-8").splitlines():
            if line.startswith("LCMS_TIMING\t"):
                _, stage, value = line.split("\t")
                timings[stage] = float(value)
        summary["runs"][name] = {"wall_seconds": elapsed, "timings": timings,
                                 "formal_accepted_psms": report["total_accepted_psms"],
                                 "unimod_rescue": report["unimod_rescue"]}
        (args.output / "comparison.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"Finished {name}: {elapsed:.2f}s", flush=True)
    def without_rescue(report):
        report = json.loads(json.dumps(report))
        report.pop("unimod_rescue", None)
        for row in report["feature_evidence"]:
            row.pop("unimod_rescue_candidates", None)
            row.pop("unimod_rescue_note", None)
        return report
    summary["first_pass_and_quantitation_identical"] = without_rescue(reports["baseline"]) == without_rescue(reports["unimod"])
    baseline = summary["runs"]["baseline"]["wall_seconds"]
    enhanced = summary["runs"]["unimod"]["wall_seconds"]
    summary["added_seconds"] = enhanced - baseline
    summary["added_percent"] = 100 * (enhanced / baseline - 1)
    (args.output / "comparison.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    if not summary["first_pass_and_quantitation_identical"]:
        raise SystemExit("First-pass outputs differ; investigate before releasing")


if __name__ == "__main__":
    main()
