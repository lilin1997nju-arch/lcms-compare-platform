#!/usr/bin/env python3
"""Run the LC-MS feature comparison MVP on local RAW/CSV inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from core.lcms_difference import build_feature_groups
from core.lcms_export import (
    payload_for_report,
    write_candidates_csv,
    write_feature_matrix_csv,
    write_json,
)
from core.lcms_feature_matching import match_features
from core.lcms_models import CandidateMz, LCMSFeature, XICTrace
from core.lcms_parser import load_lcms_directory, write_mock_centroid_csv
from core.lcms_peak_detection import detect_xic_peaks
from core.lcms_xic import extract_xic, screen_candidate_mz


REPORT_TEMPLATE = """<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>LC-MS Feature MVP</title>
  <style>
    :root { color-scheme: light; --line: #d8dde8; --ink: #172033; --muted: #667085; --accent: #2563eb; }
    body { margin: 0; font-family: Arial, "Microsoft YaHei", sans-serif; color: var(--ink); background: #f7f8fb; }
    header { padding: 16px 20px; background: #ffffff; border-bottom: 1px solid var(--line); }
    h1 { margin: 0 0 4px; font-size: 20px; }
    main { display: grid; grid-template-columns: 1.1fr 0.9fr; gap: 12px; padding: 12px; }
    section { background: #fff; border: 1px solid var(--line); border-radius: 8px; padding: 12px; min-width: 0; }
    h2 { margin: 0 0 10px; font-size: 15px; }
    canvas { width: 100%; height: 260px; border: 1px solid #eef1f6; border-radius: 6px; }
    table { width: 100%; border-collapse: collapse; font-size: 12px; }
    th, td { padding: 6px 7px; border-bottom: 1px solid #edf0f5; text-align: right; white-space: nowrap; }
    th:first-child, td:first-child { text-align: left; }
    tr:hover { background: #f3f6ff; cursor: pointer; }
    .muted { color: var(--muted); }
    .full { grid-column: 1 / -1; }
    .scroll { overflow: auto; max-height: 320px; }
    .pill { display: inline-block; padding: 2px 7px; border-radius: 999px; background: #eaf1ff; color: #1d4ed8; }
  </style>
</head>
<body>
  <header>
    <h1>LC-MS RT-m/z Feature 精确比对 MVP</h1>
    <div class="muted">四区视图：TIC/BPC 总览、Feature 列表、XIC 叠加、Feature Matrix。RAW 文件如未转换，解析状态会显示为 mock_from_vendor_raw。</div>
  </header>
  <main>
    <section>
      <h2>区域 1：TIC/BPC 总览</h2>
      <canvas id="ticCanvas" width="900" height="300"></canvas>
      <canvas id="bpcCanvas" width="900" height="300" style="margin-top:8px"></canvas>
      <div id="rawStatus" class="muted"></div>
    </section>
    <section>
      <h2>区域 2：候选 m/z / Feature 列表</h2>
      <div class="scroll"><table id="featureTable"></table></div>
    </section>
    <section>
      <h2>区域 3：XIC 叠加图 <span id="selectedFeature" class="pill"></span></h2>
      <canvas id="xicCanvas" width="900" height="300"></canvas>
    </section>
    <section>
      <h2>区域 4：Feature Matrix</h2>
      <div class="scroll"><table id="matrixTable"></table></div>
    </section>
    <section class="full">
      <h2>数据来源</h2>
      <div id="sourceTable" class="scroll"></div>
    </section>
  </main>
  <script>
    const DATA = __PAYLOAD__;
    const colors = ["#2563eb", "#dc2626", "#059669", "#9333ea", "#ea580c", "#0891b2"];
    function drawLines(canvas, series, yLabel) {
      const ctx = canvas.getContext("2d");
      const w = canvas.width, h = canvas.height;
      ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, w, h);
      const xs = series.flatMap(s => s.x);
      const ys = series.flatMap(s => s.y);
      const xmin = Math.min(...xs), xmax = Math.max(...xs);
      const ymax = Math.max(...ys, 1);
      const left = 54, right = 14, top = 12, bottom = 34;
      ctx.strokeStyle = "#d8dde8"; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(left, top); ctx.lineTo(left, h-bottom); ctx.lineTo(w-right, h-bottom); ctx.stroke();
      ctx.fillStyle = "#667085"; ctx.font = "12px Arial"; ctx.fillText("RT (min)", w/2-24, h-8); ctx.fillText(yLabel, 6, 16);
      series.forEach((s, idx) => {
        ctx.strokeStyle = colors[idx % colors.length]; ctx.lineWidth = 1.6; ctx.beginPath();
        s.x.forEach((x, i) => {
          const px = left + (x - xmin) / Math.max(xmax - xmin, 1e-9) * (w - left - right);
          const py = h - bottom - (s.y[i] / ymax) * (h - top - bottom);
          if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
        });
        ctx.stroke();
        ctx.fillStyle = colors[idx % colors.length]; ctx.fillText(s.name, left + 8, top + 15 + idx * 15);
      });
    }
    function ticSeries() {
      return DATA.sample_ids.map((sample, idx) => {
        const scans = DATA.tic_by_sample[sample] || [];
        return { name: sample, x: scans.map(p => p[0]), y: scans.map(p => p[1]) };
      });
    }
    function bpcSeries() {
      return DATA.sample_ids.map((sample, idx) => {
        const scans = DATA.bpc_by_sample[sample] || [];
        return { name: sample, x: scans.map(p => p[0]), y: scans.map(p => p[1]) };
      });
    }
    function renderFeatureTable() {
      const rows = DATA.feature_groups;
      const html = ["<tr><th>feature_group_id</th><th>m/z</th><th>RT</th><th>samples</th><th>max_area</th><th>fold</th><th>type</th><th>confidence</th></tr>"];
      rows.forEach(g => {
        const maxArea = Math.max(...Object.values(g.area_by_sample));
        html.push(`<tr data-group="${g.feature_group_id}"><td>${g.feature_group_id}</td><td>${g.representative_mz.toFixed(5)}</td><td>${g.representative_rt.toFixed(3)}</td><td>${g.sample_count}</td><td>${maxArea.toFixed(1)}</td><td>${g.max_fold_change ? g.max_fold_change.toFixed(2) : ""}</td><td>${g.difference_type}</td><td>${g.confidence_score.toFixed(2)}</td></tr>`);
      });
      document.getElementById("featureTable").innerHTML = html.join("");
      document.querySelectorAll("#featureTable tr[data-group]").forEach(row => {
        row.addEventListener("click", () => renderXic(row.dataset.group));
      });
    }
    function renderMatrix() {
      const header = ["<tr><th>feature</th><th>type</th>", ...DATA.sample_ids.map(s => `<th>${s}</th>`), "</tr>"].join("");
      const rows = DATA.feature_groups.map(g => `<tr><td>${g.feature_group_id}</td><td>${g.difference_type}</td>${DATA.sample_ids.map(s => `<td>${(g.area_by_sample[s] || 0).toFixed(1)}</td>`).join("")}</tr>`);
      document.getElementById("matrixTable").innerHTML = header + rows.join("");
    }
    function renderXic(groupId) {
      const traces = DATA.xics_by_group[groupId] || [];
      document.getElementById("selectedFeature").textContent = groupId || "";
      drawLines(document.getElementById("xicCanvas"), traces.map(t => ({ name: t.sample_id, x: t.rt_array, y: t.intensity_array })), "XIC");
    }
    function renderSources() {
      document.getElementById("rawStatus").textContent = DATA.raw_files.map(f => `${f.sample_id}: ${f.parser_status}`).join(" | ");
      const rows = DATA.raw_files.map(f => `<tr><td>${f.sample_id}</td><td>${f.file_name}</td><td>${f.data_format}</td><td>${f.scan_count}</td><td>${f.rt_min.toFixed(2)}-${f.rt_max.toFixed(2)}</td><td>${f.mz_min.toFixed(1)}-${f.mz_max.toFixed(1)}</td><td>${f.parser_status}</td></tr>`);
      document.getElementById("sourceTable").innerHTML = `<table><tr><th>sample</th><th>file</th><th>format</th><th>scans</th><th>RT</th><th>m/z</th><th>parser</th></tr>${rows.join("")}</table>`;
    }
    renderSources();
    drawLines(document.getElementById("ticCanvas"), ticSeries(), "TIC");
    drawLines(document.getElementById("bpcCanvas"), bpcSeries(), "BPC");
    renderFeatureTable();
    renderMatrix();
    renderXic(DATA.feature_groups[0]?.feature_group_id);
  </script>
</body>
</html>
"""


def traces_for_tic(scans_by_sample: dict[str, list[object]]) -> dict[str, list[list[float]]]:
    return {
        sample_id: [[scan.rt, scan.tic] for scan in scans]
        for sample_id, scans in scans_by_sample.items()
    }


def traces_for_bpc(scans_by_sample: dict[str, list[object]]) -> dict[str, list[list[float]]]:
    return {
        sample_id: [[scan.rt, scan.base_peak_intensity] for scan in scans]
        for sample_id, scans in scans_by_sample.items()
    }


def run_pipeline(args: argparse.Namespace) -> dict[str, object]:
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    raw_files, scans_by_sample = load_lcms_directory(input_dir, project_id=args.project_id)
    sample_ids = [raw_file.sample_id for raw_file in raw_files]
    reference_samples = {sample for sample in sample_ids if sample.startswith(args.reference_prefix)}
    target_samples = {sample for sample in sample_ids if sample.startswith(args.target_prefix)}
    if not reference_samples or not target_samples:
        raise ValueError("Could not infer both reference and target sample groups")

    candidates_by_sample: dict[str, list[CandidateMz]] = {}
    all_features: list[LCMSFeature] = []
    xic_lookup: dict[tuple[str, float], XICTrace] = {}
    for sample_id, scans in scans_by_sample.items():
        candidates = screen_candidate_mz(
            scans,
            args.rt_start,
            args.rt_end,
            mz_min=args.mz_min,
            mz_max=args.mz_max,
            mz_tolerance_ppm=args.mz_tolerance_ppm,
            intensity_threshold=args.intensity_threshold,
            min_scan_count=args.min_scan_count,
            top_n_mz=args.top_n_mz,
        )
        candidates_by_sample[sample_id] = candidates
        write_candidates_csv(output_dir / f"candidates_{sample_id}.csv", candidates)
        for candidate in candidates:
            xic = extract_xic(scans, candidate.candidate_mz, args.mz_tolerance_ppm)
            xic_lookup[(sample_id, round(candidate.candidate_mz, 4))] = xic
            peaks = detect_xic_peaks(
                xic,
                smoothing_window=args.smoothing_window,
                min_peak_height=args.min_peak_height,
                min_peak_area=args.min_peak_area,
                min_snr=args.min_snr,
                min_peak_width=args.min_peak_width,
                max_peak_width=args.max_peak_width,
            )
            all_features.extend(peaks[:1])

    matched = match_features(
        all_features,
        sample_ids,
        mz_tolerance_ppm=args.mz_tolerance_ppm,
        rt_tolerance_min=args.rt_tolerance_min,
    )
    groups = build_feature_groups(
        matched,
        sample_ids,
        reference_samples=reference_samples,
        target_samples=target_samples,
        rt_shift_threshold_min=args.rt_shift_threshold_min,
        fold_change_threshold=args.fold_change_threshold,
    )
    groups.sort(key=lambda group: (group.difference_type == "common_feature", -group.confidence_score, group.representative_rt))

    xics_by_group: dict[str, list[XICTrace]] = {}
    for group in groups:
        traces: list[XICTrace] = []
        for sample_id in sample_ids:
            scans = scans_by_sample[sample_id]
            traces.append(extract_xic(scans, group.representative_mz, args.mz_tolerance_ppm))
        xics_by_group[group.feature_group_id] = traces

    payload = payload_for_report(raw_files, candidates_by_sample, all_features, groups, xics_by_group, sample_ids)
    payload["tic_by_sample"] = traces_for_tic(scans_by_sample)
    payload["bpc_by_sample"] = traces_for_bpc(scans_by_sample)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "lcms_mvp_results.json", payload)
    write_feature_matrix_csv(output_dir / "feature_matrix.csv", groups, sample_ids)
    raw_paths = [Path(raw_file.file_path) for raw_file in raw_files if raw_file.data_format == "raw"]
    if raw_paths:
        write_mock_centroid_csv(output_dir / "mock_centroid_scans.csv", raw_paths)
    report_html = REPORT_TEMPLATE.replace("__PAYLOAD__", Path(output_dir / "lcms_mvp_results.json").read_text(encoding="utf-8"))
    (output_dir / "lcms_mvp_report.html").write_text(report_html, encoding="utf-8")
    return {
        "sample_count": len(sample_ids),
        "feature_count": len(all_features),
        "feature_group_count": len(groups),
        "output_dir": str(output_dir.resolve()),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run LC-MS feature comparison MVP.")
    parser.add_argument("--input-dir", default="lcms_feature_mvp/data/raw/zenodo_5005513")
    parser.add_argument("--output-dir", default="lcms_feature_mvp/outputs")
    parser.add_argument("--project-id", default="zenodo_5005513")
    parser.add_argument("--reference-prefix", default="MabThera")
    parser.add_argument("--target-prefix", default="Reditux")
    parser.add_argument("--rt-start", type=float, default=4.0)
    parser.add_argument("--rt-end", type=float, default=8.8)
    parser.add_argument("--mz-min", type=float, default=200.0)
    parser.add_argument("--mz-max", type=float, default=2000.0)
    parser.add_argument("--mz-tolerance-ppm", type=float, default=20.0)
    parser.add_argument("--intensity-threshold", type=float, default=1500.0)
    parser.add_argument("--min-scan-count", type=int, default=4)
    parser.add_argument("--top-n-mz", type=int, default=30)
    parser.add_argument("--smoothing-window", type=int, default=5)
    parser.add_argument("--min-peak-height", type=float, default=2500.0)
    parser.add_argument("--min-peak-area", type=float, default=200.0)
    parser.add_argument("--min-snr", type=float, default=3.0)
    parser.add_argument("--min-peak-width", type=float, default=0.04)
    parser.add_argument("--max-peak-width", type=float, default=0.8)
    parser.add_argument("--rt-tolerance-min", type=float, default=0.2)
    parser.add_argument("--rt-shift-threshold-min", type=float, default=0.2)
    parser.add_argument("--fold-change-threshold", type=float, default=1.5)
    return parser.parse_args()


def main() -> int:
    result = run_pipeline(parse_args())
    print(f"Samples: {result['sample_count']}")
    print(f"Features: {result['feature_count']}")
    print(f"Feature groups: {result['feature_group_count']}")
    print(f"Outputs: {result['output_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
