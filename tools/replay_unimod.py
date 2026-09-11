"""Repeat just the added search on saved first-pass results (no task writes)."""
import argparse
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lcms_feature_mvp"))
from core.lcms_msms import generate_sequence_inference_candidates
from core.lcms_parser import read_mzml
from core.lcms_unimod import search_unimod_rescue
from run_peak_first_compare import PEAK_FIRST_TEMPLATE


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--task", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads((args.comparison / "baseline/lcms_msms_identifications.json").read_text(encoding="utf-8"))
    scans = []
    for path in sorted((args.task / "mzML").glob("*.mzML")):
        scans.extend(read_mzml(path, source["project_id"], ms_levels=(2,))[1])
    params = source["parameters"]
    candidates = generate_sequence_inference_candidates(
        source["chains"], max_missed_cleavages=params["sequence_inference_max_missed_cleavages"],
        max_terminal_trim=params["sequence_inference_max_terminal_trim"], enzyme=params["enzyme"],
        carbamidomethyl_cys=params["fixed_modification"] != "None")
    runs = []
    first_hits = None
    for repeat in range(3):
        features = deepcopy(source["feature_evidence"])
        summary = search_unimod_rescue(features, scans, candidates, source["chains"],
                                       precursor_ppm=params["precursor_ppm"] * 2,
                                       fragment_ppm=params["fragment_ppm"])
        hits = {f["feature_group_id"]: f["unimod_rescue_candidates"] for f in features if f.get("unimod_rescue_candidates")}
        if first_hits is None:
            first_hits = hits
        elif hits != first_hits:
            raise RuntimeError("Repeated second-pass candidates changed")
        runs.append(summary)
        print(f"Replay {repeat+1}: {summary['elapsed_seconds']:.3f}s, {len(hits)} Features", flush=True)
    ab_report = json.loads((args.comparison / "unimod/lcms_msms_identifications.json").read_text(encoding="utf-8"))
    ab_hits = {f["feature_group_id"]: f["unimod_rescue_candidates"] for f in ab_report["feature_evidence"] if f.get("unimod_rescue_candidates")}
    counts = Counter(r["unimod_title"] for rows in first_hits.values() for r in rows)
    result = {"runs": runs, "median_seconds": statistics.median(r["elapsed_seconds"] for r in runs),
              "identical_to_ab_candidates": first_hits == ab_hits,
              "candidate_names": dict(counts.most_common()),
              "example_feature": next(iter(first_hits), None)}
    (args.comparison / "replay.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    for name in ("baseline", "unimod"):
        (args.comparison / name / "lcms_peak_first_compare.html").write_text(PEAK_FIRST_TEMPLATE, encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "runs"}, ensure_ascii=True, indent=2))
    if not result["identical_to_ab_candidates"]:
        raise SystemExit("Updated code differs from A/B; rerun comparison")


if __name__ == "__main__":
    main()
