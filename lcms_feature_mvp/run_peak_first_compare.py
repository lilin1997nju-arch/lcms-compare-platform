#!/usr/bin/env python3
"""Generate the standalone Peak-first LC-MS compare workbench."""

from __future__ import annotations

import argparse
import json
import sqlite3
from contextlib import closing
from pathlib import Path

from core.lcms_parser import load_lcms_directory
from core.lcms_msms import build_ms1_component_groups
from core.lcms_peak_first import PeakFirstParams, prepare_peak_first_payload
from core.lcms_workbench import spectrum_payload


PEAK_FIRST_TEMPLATE = r"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>LC-MS 峰优先差异分析 / Peak-first Compare</title>
  <style>
    :root {
      --ink:#17253d;
      --muted:#718096;
      --line:#dce6f2;
      --canvas:#f3f7fb;
      --card:#fff;
      --navy:#0a1730;
      --blue:#2563eb;
      --cyan:#0ea5b7;
      --green:#15966b;
      --amber:#c77710;
      --red:#c2414d;
      --shadow:0 18px 45px rgba(31,58,92,.08);
    }
    * { box-sizing:border-box; }
    html { min-width:320px; background:var(--canvas); }
    body { margin:0; min-height:100vh; font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI","Microsoft YaHei",sans-serif; color:var(--ink); background:var(--canvas); line-height:1.5; overflow-x:hidden; }
    body::before { content:""; position:fixed; z-index:-1; width:620px; height:620px; right:-230px; top:96px; border-radius:50%; background:radial-gradient(circle,rgba(38,99,235,.08),rgba(38,99,235,0) 68%); pointer-events:none; }
    header { position:sticky; top:0; z-index:20; min-height:70px; padding:14px clamp(18px,3.5vw,54px); display:grid; grid-template-columns:minmax(0,1fr) auto; gap:18px; align-items:center; color:#fff; background:rgba(10,23,48,.97); border-bottom:1px solid rgba(148,191,255,.18); box-shadow:0 8px 24px rgba(10,23,48,.16); backdrop-filter:blur(14px); }
    header h1 { margin:0; color:#f6faff; font-size:18px; font-weight:780; letter-spacing:-.01em; line-height:1.25; }
    header > a { padding:8px 12px; color:#d8e8ff; font-size:12px; font-weight:700; text-decoration:none; border:1px solid rgba(153,192,237,.2); border-radius:9px; background:rgba(255,255,255,.07); white-space:nowrap; }
    header > a:hover { color:#fff; background:rgba(255,255,255,.12); }
    main { width:min(1540px,100%); margin:0 auto; display:grid; grid-template-columns:minmax(0,1.75fr) minmax(360px,.75fr); gap:20px; padding:28px clamp(14px,3.5vw,52px) 58px; }
    section { position:relative; min-width:0; padding:22px; overflow:hidden; border:1px solid var(--line); border-radius:18px; background:var(--card); box-shadow:var(--shadow); }
    section::before { content:""; position:absolute; left:0; top:0; width:100%; height:3px; background:linear-gradient(90deg,rgba(37,99,235,.8),rgba(14,165,183,.32),transparent 80%); opacity:.72; }
    h2 { margin:0 0 14px; color:#182945; font-size:16px; font-weight:780; letter-spacing:-.015em; }
    h2::first-letter { color:var(--blue); }
    canvas { display:block; width:100%; border:1px solid #dfe9f3; border-radius:12px; background:#fff; box-shadow:inset 0 1px 0 rgba(255,255,255,.7); }
    .wide { grid-column:1 / -1; }
    .left { grid-column:1; }
    .right { grid-column:2; }
    .controls { display:flex; gap:8px; flex-wrap:wrap; align-items:center; margin-bottom:11px; padding:9px 11px; border:1px solid #e3ebf4; border-radius:12px; background:#f8fafd; }
    .chrom-legend { display:flex; gap:14px; flex-wrap:wrap; align-items:center; margin:0 0 10px; padding:8px 11px; color:#536a86; font-size:11px; border:1px solid #e3ebf4; border-radius:10px; background:#fbfdff; }
    .chrom-legend-item { display:inline-flex; gap:6px; align-items:center; white-space:nowrap; }
    .chrom-legend-swatch { width:24px; height:3px; border-radius:3px; display:inline-block; }
    .chrom-legend-swatch.reference { background:#1f77b4; }
    .chrom-legend-swatch.test { background:#e74c3c; }
    [hidden] { display:none !important; }
    .report-settings { overflow:visible; }
    .report-settings details { margin:0; }
    .report-settings summary { display:flex; gap:12px; align-items:center; justify-content:space-between; cursor:pointer; list-style:none; color:#294866; font-size:15px; font-weight:780; }
    .report-settings summary::-webkit-details-marker { display:none; }
    .report-settings summary::after { content:"⌄"; color:#6b86a4; font-size:18px; transition:transform .18s ease; }
    .report-settings details:not([open]) summary::after { transform:rotate(-90deg); }
    .settings-intro { margin:5px 0 14px; }
    .settings-layout { display:grid; grid-template-columns:minmax(180px,.85fr) minmax(250px,1fr) minmax(340px,1.45fr); gap:13px; }
    .settings-group { min-width:0; margin:0; padding:12px; border:1px solid #e1eaf3; border-radius:12px; background:#f8fafd; }
    .settings-group h3 { margin:0 0 9px; color:#294866; font-size:12px; }
    .settings-options { display:grid; gap:7px; }
    .settings-check, .settings-select { display:flex; gap:7px; align-items:center; color:#536a86; font-size:11px; line-height:1.35; }
    .settings-check input { accent-color:#2563eb; }
    .settings-select select { min-width:0; flex:1; padding:5px 7px; font-size:11px; }
    .settings-color-list { display:grid; gap:6px; margin-top:11px; padding-top:10px; border-top:1px solid #e1eaf3; }
    .settings-color-row { display:flex; gap:7px; align-items:center; min-width:0; color:#536a86; font-size:11px; }
    .settings-color-row input[type="color"] { width:28px; height:24px; padding:2px; flex:0 0 auto; cursor:pointer; }
    .settings-color-row span { overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
    .settings-table-group { padding:7px 0; border-top:1px solid #e1eaf3; }
    .settings-table-group:first-child { padding-top:0; border-top:0; }
    .settings-table-title { margin:0 0 5px; color:#536a86; font-size:11px; font-weight:750; }
    .settings-column-grid { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:4px 8px; }
    .settings-column-grid label { display:flex; min-width:0; gap:5px; align-items:center; color:#718096; font-size:10px; }
    .settings-column-grid span { overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
    .settings-footer { display:flex; gap:10px; align-items:center; justify-content:flex-end; margin-top:12px; }
    .settings-footer .note { margin-right:auto; }
    .settings-reset { padding:6px 9px; color:#536a86; font-size:11px; }
    .report-submodule { min-width:0; }
    #mzXicSection.single-module .xic-grid { grid-template-columns:1fr; }
    #globalFeatureTable th.global-component-column, #globalFeatureTable td.global-component-column { width:350px; min-width:350px; max-width:350px; white-space:normal; overflow-wrap:anywhere; word-break:break-word; vertical-align:top; }
    .controls label { display:inline-flex; gap:5px; align-items:center; color:#536a86; font-size:11px; font-weight:650; }
    button, select, input { border:1px solid #cfdcea; border-radius:9px; padding:8px 10px; color:var(--ink); font:inherit; font-size:12px; background:#fff; outline:none; transition:border-color .18s ease,box-shadow .18s ease,background .18s ease,transform .18s ease; }
    button { cursor:pointer; font-weight:700; }
    button:hover:not(:disabled) { transform:translateY(-1px); border-color:#a9c5e5; background:#f5f9ff; }
    button:disabled, input:disabled { cursor:default; opacity:.45; }
    select:focus, input:focus { border-color:#6ca5f5; box-shadow:0 0 0 3px rgba(37,99,235,.12); }
    table { border-collapse:separate; border-spacing:0; width:100%; font-size:11px; }
    th, td { padding:9px 9px; white-space:nowrap; text-align:right; border-bottom:1px solid #e9eff6; }
    th { position:sticky; top:0; z-index:2; color:#7b8da5; font-size:9px; font-weight:850; letter-spacing:.08em; text-transform:uppercase; background:#f8fafd; }
    th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) { text-align:left; }
    tr:last-child td { border-bottom:0; }
    tr[data-peak], tr[data-mz], tr[data-global-feature], tr[data-saved-feature] { cursor:pointer; }
    tr:hover { background:#f7fbff; }
    tr.selected { background:#eef7ff; box-shadow:inset 3px 0 0 #3182ce; }
    .scroll { overflow:auto; max-height:360px; border:1px solid #e5edf5; border-radius:12px; }
    .tall-scroll { max-height:420px; }
    .note { color:var(--muted); font-size:11px; line-height:1.55; }
    .pan { width:220px; }
    .xic-grid { display:grid; grid-template-columns:minmax(260px,.56fr) minmax(0,1.44fr); gap:18px; align-items:start; }
    .xic-grid canvas { height:auto; }
    .feature-map { height:auto; }
    .feature-analysis-grid { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1fr); gap:18px; align-items:start; }
    .feature-analysis-pane { min-width:0; }
    .feature-ms2-heading { display:flex; gap:10px; align-items:center; justify-content:space-between; margin-bottom:8px; }
    .feature-ms2-heading h3 { margin:0; color:#294866; font-size:13px; }
    .feature-ms2-heading a { color:#205ab6; font-size:11px; }
    .feature-ms2-detail { min-height:38px; margin-bottom:8px; padding:9px 11px; color:#536a86; font-size:11px; line-height:1.5; border:1px solid #e3ebf4; border-radius:10px; background:#f8fafd; }
    .feature-ms2-status { display:inline-block; margin-right:6px; padding:2px 7px; border-radius:10px; background:#eef2ff; color:#344054; font-weight:700; }
    .feature-ms2-status.identified { background:#e7f8f0; color:#087b55; }
    .feature-ms2-status.unresolved { background:#fff5e8; color:#a45d0c; }
    .feature-ms2-status.missing { background:#f1f4f7; color:#667085; }
    #featureMs2Annotations { margin-top:9px; padding:9px 11px; color:#6a4b0b; border:1px solid #f0d790; border-radius:10px; background:#fffbeb; }
    .agent-annotation { margin-top:7px; padding-top:7px; border-top:1px solid #f5e7b8; }
    .mod-quant-summary { display:flex; gap:8px; flex-wrap:wrap; margin:8px 0 11px; }
    .digestion-warning { color:#9a3412; font-weight:750; }
    .mod-view-switch { display:inline-flex; overflow:hidden; border:1px solid #b9cbe4; border-radius:8px; background:#f8fbff; }
    .mod-view-switch button { min-width:82px; border:0; border-radius:0; color:#536a86; background:transparent; }
    .mod-view-switch button.active { color:#fff; background:#2563eb; }
    .mod-quant-badge { display:inline-flex; padding:3px 8px; border-radius:999px; color:#175cd3; background:#eff8ff; font-size:11px; font-weight:750; }
    .mod-quant-badge.formal { color:#067647; background:#ecfdf3; }
    .mod-quant-badge.high-value { color:#9a3412; background:#fff2df; }
    .mod-quant-table { min-width:1050px; table-layout:fixed; }
    .mod-quant-table td { vertical-align:top; }
    .mod-quant-table td.wrap { white-space:normal; text-align:left; }
    .mod-quant-table th:nth-child(1), .mod-quant-table td:nth-child(1) { width:38px; }
    .mod-quant-table th:nth-child(2), .mod-quant-table td:nth-child(2) { width:48px; }
    .mod-quant-table th:nth-child(3), .mod-quant-table td:nth-child(3) { width:190px; }
    .mod-quant-table th:nth-child(4), .mod-quant-table td:nth-child(4) { width:230px; }
    .mod-quant-table th:nth-child(6), .mod-quant-table td:nth-child(6) { width:135px; }
    .mod-quant-table th:nth-child(7), .mod-quant-table td:nth-child(7) { width:70px; }
    .mod-quant-event { color:#172033; font-weight:780; }
    .mod-quant-site { margin-top:2px; color:#667085; font-size:10px; }
    .mod-quant-toggle { width:28px; height:26px; padding:0; border-radius:7px; color:#175cd3; background:#eff8ff; }
    .mod-quant-conclusion { white-space:normal; color:#344054; text-align:left; }
    .mod-quant-conclusion.high-value { color:#9a3412; font-weight:750; }
    .mod-quant-detail-row > td { padding:0 !important; background:#f8fbff; }
    .mod-proteolytic-group-row > td { background:#eef5ff; border-top:2px solid #b8d1f5; }
    .mod-proteolytic-group-row.high-value > td { background:#fff8e8; border-top-color:#efb65b; }
    .mod-proteolytic-member-row > td { background:#fbfdff; }
    .mod-proteolytic-member-row td:nth-child(3) { padding-left:20px; }
    .mod-group-sequences { margin-top:4px; color:#667085; font-size:10px; line-height:1.45; }
    .mod-quant-detail { padding:14px 16px 16px 54px; border-top:1px solid #dbeafe; }
    .mod-quant-detail-grid { display:grid; grid-template-columns:minmax(300px,.8fr) minmax(520px,1.5fr); gap:16px; align-items:start; }
    .mod-composition-panel, .mod-form-panel { padding:11px; border:1px solid #dbe3ef; border-radius:10px; background:#fff; }
    .mod-detail-title { margin-bottom:8px; color:#344054; font-size:11px; font-weight:780; }
    .mod-composition-row { display:grid; grid-template-columns:86px 1fr; gap:8px; align-items:center; margin:8px 0; }
    .mod-composition-sample { overflow:hidden; color:#475467; font-size:10px; font-weight:750; text-overflow:ellipsis; }
    .mod-composition-bar { display:flex; height:24px; overflow:hidden; border:1px solid #d0d5dd; border-radius:6px; background:#f2f4f7; }
    .mod-composition-scale { display:flex; width:0; height:100%; overflow:hidden; transition:width .58s cubic-bezier(.22,.8,.3,1); }
    .mod-composition-segment { min-width:0; height:100%; border-right:1px solid rgba(255,255,255,.7); opacity:.9; transition:opacity .35s ease; }
    .mod-composition-scale.animated .mod-composition-segment { opacity:1; }
    .mod-form-matrix { width:100%; min-width:0; font-size:10px; }
    .mod-form-matrix th, .mod-form-matrix td { padding:6px 7px; text-align:right; white-space:nowrap; }
    .mod-form-matrix th:first-child, .mod-form-matrix td:first-child { width:auto; text-align:left; white-space:normal; }
    .mod-form-label { display:flex; gap:6px; align-items:flex-start; }
    .mod-form-swatch { flex:0 0 9px; width:9px; height:9px; margin-top:3px; border-radius:2px; }
    .mod-form-link { padding:0; border:0; color:#175cd3; background:transparent; font:inherit; text-align:left; text-decoration:underline; text-decoration-style:dotted; text-underline-offset:3px; cursor:pointer; }
    .mod-form-link:hover { color:#0b4cac; }
    .mod-form-link-count { flex:0 0 auto; margin-left:3px; color:#667085; font-size:9px; white-space:nowrap; }
    .mod-form-matrix tr.mod-form-selected td { background:#eef6ff; box-shadow:inset 0 1px #b9d4fb,inset 0 -1px #b9d4fb; }
    .mod-delta-up { color:#b42318; font-weight:750; }
    .mod-delta-down { color:#175cd3; font-weight:750; }
    .mod-delta-stable { color:#667085; }
    .mod-feature-details { margin-top:10px; }
    .mod-feature-details summary { cursor:pointer; color:#475467; font-size:10px; font-weight:700; }
    .mod-quant-feature { margin:2px 3px 2px 0; padding:2px 5px; color:#175cd3; font-size:10px; background:#eff8ff; }
    @media(max-width:1100px){ .mod-quant-detail-grid{grid-template-columns:1fr;} }
    @media(prefers-reduced-motion:reduce){ .mod-composition-scale,.mod-composition-segment{transition:none;} }
    @media (max-width:1100px) { main { grid-template-columns:1fr; } .left,.right { grid-column:1; } .feature-analysis-grid { grid-template-columns:1fr; } .settings-layout { grid-template-columns:repeat(2,minmax(0,1fr)); } .settings-layout .settings-group:last-child { grid-column:1 / -1; } }
    .tooltip { position:fixed; display:none; z-index:20; max-width:360px; padding:10px 12px; color:#284360; border:1px solid #cfdcea; border-radius:10px; background:rgba(255,255,255,.97); box-shadow:0 12px 28px rgba(15,23,42,.15); font-size:11px; line-height:1.5; pointer-events:none; }
    .tooltip-title { display:flex; gap:7px; align-items:center; margin-bottom:5px; color:#1d3554; font-weight:800; }
    .tooltip-muted { color:#8a9bb0; font-size:10px; }
    .ion-badge { display:inline-flex; margin:3px 0 5px; padding:2px 7px; color:#2456a6; background:#eef5ff; border:1px solid #d4e4fa; border-radius:999px; font-weight:750; }
    .mini-table { margin:7px 0 11px; max-height:150px; overflow:auto; border:1px solid #e5edf5; border-radius:9px; }
    .status-high_consistency { color:#15966b; font-weight:750; }
    .status-chromatogram_consistent_spectrum_changed,
    .status-need_local_rt_mz_analysis,
    .status-chromatogram_changed,
    .status-possible_rt_shift,
    .status-local_alignment_rejected_no_clear_peak { color:#c2414d; font-weight:750; }
    .local-align-skipped_no_clear_peak { color:#c2414d; font-weight:750; }
    .local-align-capped_possible_rt_shift { color:#b26b13; font-weight:750; }
    .local-align-applied { color:#205ab6; }
    .local-align-not_needed { color:#718096; }
    th.sortable { cursor:pointer; user-select:none; color:#607893; }
    th.sortable:hover { color:#205ab6; background:#eef5ff; }
    .section-title-row { display:flex; align-items:center; justify-content:space-between; gap:12px; }
    .section-title-row h2 { margin-bottom:14px; }
    .export-button { padding:8px 12px; color:#205ab6; border-color:#c9dcf5; background:#eff6ff; }
    .export-button::before { content:"↓"; margin-right:6px; }
    .feature-map-pane { display:flex; min-width:0; flex-direction:column; padding-top:48px; }
    .feature-map-legend { order:0; display:flex; min-height:30px; align-items:center; justify-content:center; gap:8px; padding:0 0 10px; color:#718096; font-size:10px; }
    .feature-map-legend .legend-gradient { display:inline-block; width:min(360px,48%); height:12px; border:1px solid rgba(23,37,61,.12); border-radius:4px; }
    .feature-map-legend .legend-gradient.direction { background:linear-gradient(90deg,#1e40af 0%,#93c5fd 28%,#e5e7eb 50%,#fca5a5 72%,#991b1b 100%); }
    .feature-map-legend .legend-gradient.magnitude { background:linear-gradient(90deg,#16a34a 0%,#eab308 30%,#f97316 62%,#dc2626 82%,#581c87 100%); }
    .feature-map-legend .legend-label { white-space:nowrap; }
    .feature-map-legend .legend-caption { margin-left:8px; color:#8a9bb0; white-space:nowrap; }
    .feature-map-pane canvas { order:1; }
    .feature-ms2-heading { min-height:28px; }
    .sequence-overview-pane { min-width:0; }
    .sample-order { display:inline-flex; gap:4px; align-items:center; justify-content:flex-end; min-width:42px; }
    .sample-dot { width:11px; height:11px; border-radius:50%; display:inline-block; border:1px solid rgba(15,23,42,.25); box-shadow:0 0 0 1px rgba(255,255,255,.75) inset; }
    .direction-up { color:#b42318; font-weight:600; }
    .direction-down { color:#175cd3; font-weight:600; }
    .direction-equal { color:#667085; }
    .pair-legend { display:inline-flex; gap:10px; flex-wrap:wrap; align-items:center; }
    .pair-legend > span { display:inline-flex; gap:4px; align-items:center; }
    .component-row { background:#f8fafc; font-weight:600; }
    .component-row.selected, .component-child.selected { background:#fff3d8; }
    .component-child { background:#fff; color:#475467; }
    .component-child td:first-child { color:#98a2b3; }
    .component-cell { text-align:left !important; }
    #globalFeatureTable { min-width:2200px; }
    #globalFeatureTable th:nth-child(2), #globalFeatureTable td:nth-child(2) { width:350px; min-width:350px; max-width:350px; white-space:normal; overflow-wrap:anywhere; word-break:break-word; vertical-align:top; }
    #globalFeatureTable .component-label, #globalFeatureTable .component-meta { white-space:normal; overflow-wrap:anywhere; word-break:break-word; }
    .component-toggle { width:24px; min-width:24px; padding:1px 4px; margin-right:5px; border-color:#d0d5dd; line-height:18px; }
    .component-toggle-placeholder { display:inline-block; width:29px; }
    .component-label { color:#172033; }
    .component-meta { margin-left:7px; color:#667085; font-size:11px; font-weight:400; }
    .component-tree { display:inline-block; width:24px; color:#98a2b3; text-align:center; }
    .identification-badge { display:inline-block; margin-left:6px; padding:1px 5px; border-radius:10px; background:#ecfdf3; color:#027a48; font-size:10px; font-weight:600; }
    .inference-badge { display:inline-block; margin-left:6px; padding:1px 5px; border-radius:10px; background:#eff8ff; color:#175cd3; font-size:10px; font-weight:600; }
    .candidate-badge { display:inline-block; margin-left:6px; padding:1px 5px; border-radius:10px; background:#fff5e8; color:#a45d0c; font-size:10px; font-weight:600; }
    .structure-grid { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1fr); gap:18px; align-items:start; }
    .structure-pane { min-width:0; padding:15px; border:1px solid #e2eaf3; border-radius:13px; background:#f8fafd; }
    .structure-pane h3 { margin:0 0 10px; color:#294866; font-size:13px; }
    .structure-toolbar { display:flex; gap:8px; flex-wrap:wrap; align-items:center; margin-bottom:10px; padding:9px 11px; border:1px solid #e3ebf4; border-radius:11px; background:#fff; }
    .structure-toolbar input[type="text"] { width:90px; text-transform:uppercase; }
    .structure-toolbar input[type="file"] { max-width:240px; }
    .structure-upload-button { padding:5px 8px; }
    .structure-status { min-height:34px; margin:7px 0; line-height:1.5; }
    .structure-badge { display:inline-block; margin-right:6px; padding:2px 7px; border-radius:10px; font-size:11px; font-weight:600; background:#eff8ff; color:#175cd3; }
    .structure-badge.mapped { background:#ecfdf3; color:#027a48; }
    .structure-badge.candidate { background:#fff4e5; color:#b54708; }
    .structure-badge.missing { background:#f2f4f7; color:#667085; }
    .global-sequence-overview { margin:15px 0 17px; padding:14px; border:1px solid #dfe8f2; border-radius:13px; background:linear-gradient(135deg,#f8fafd,#f2f8fb); }
    .global-sequence-overview h3 { margin:0 0 8px; color:#294866; font-size:14px; }
    .global-sequence-tracks { display:grid; grid-template-columns:1fr; gap:10px; margin-top:8px; }
    .global-chain-card { min-width:0; padding:11px; border:1px solid #e2eaf3; border-radius:10px; background:#fff; }
    .global-chain-heading { display:flex; gap:10px; flex-wrap:wrap; align-items:baseline; margin-bottom:6px; }
    .global-chain-heading b { color:#101828; }
    .global-chain-summary { color:#667085; font-size:12px; }
    .global-chain-track { max-height:230px; overflow:auto; font:12.5px/1.72 Consolas,"Courier New",monospace; }
    .sequence-track { max-height:410px; overflow:auto; padding:11px; border:1px solid #e2eaf3; border-radius:10px; background:white; font:13px/1.75 Consolas,"Courier New",monospace; }
    .sequence-block { display:flex; align-items:flex-start; gap:8px; margin-bottom:2px; white-space:nowrap; }
    .sequence-position { width:42px; color:#98a2b3; text-align:right; user-select:none; }
    .sequence-residues { letter-spacing:1.5px; }
    .sequence-residue { display:inline-block; min-width:10px; text-align:center; border-radius:2px; }
    .sequence-residue.global-difference { color:#fff; cursor:pointer; font-weight:700; }
    .sequence-residue.diff-higher { background:#dc2626; }
    .sequence-residue.diff-lower { background:#2563eb; }
    .sequence-residue.diff-conflict { background:linear-gradient(135deg,#dc2626 0 46%,#fff 46% 54%,#2563eb 54% 100%); box-shadow:0 0 0 1px #6d28d9 inset; }
    .sequence-residue.diff-overlap-selected { border-top:2px solid #111827; border-bottom:2px solid #111827; border-radius:0; opacity:1; position:relative; z-index:2; }
    .sequence-residue.diff-overlap-start { border-left:2px solid #111827; border-radius:4px 0 0 4px; }
    .sequence-residue.diff-overlap-end { border-right:2px solid #111827; border-radius:0 4px 4px 0; }
    .sequence-residue.diff-tentative { opacity:.62; }
    .sequence-residue.diff-selected { outline:2px solid #111827; outline-offset:1px; opacity:1; position:relative; z-index:1; }
    .sequence-residue.peptide { color:#fff; font-weight:700; }
    .sequence-residue.modified { box-shadow:0 0 0 2px #f59e0b inset; }
    .sequence-legend { display:flex; gap:12px; flex-wrap:wrap; margin:7px 0 0; }
    .sequence-legend span { display:inline-flex; align-items:center; gap:5px; }
    .legend-swatch { width:12px; height:12px; border-radius:3px; background:#dc2626; display:inline-block; }
    .legend-swatch.higher { background:#dc2626; }
    .legend-swatch.lower { background:#2563eb; }
    .legend-swatch.conflict { background:linear-gradient(135deg,#dc2626 0 46%,#fff 46% 54%,#2563eb 54% 100%); box-shadow:0 0 0 1px #6d28d9 inset; }
    .legend-swatch.tentative { background:#64748b; opacity:.62; }
    .legend-swatch-site { width:12px; height:12px; border-radius:50%; background:#f59e0b; display:inline-block; }
    .sequence-overlap-details { margin-top:12px; padding:12px; border:1px solid #c4b5fd; border-radius:11px; background:#faf8ff; }
    .sequence-overlap-heading { display:flex; gap:10px; align-items:center; justify-content:space-between; margin-bottom:7px; }
    .sequence-overlap-heading b { color:#4c1d95; }
    .sequence-overlap-heading button { padding:3px 8px; }
    .sequence-overlap-note { margin-bottom:9px; color:#667085; font-size:12px; }
    .sequence-overlap-scroll { overflow:auto; padding:9px; border:1px solid #e4def7; border-radius:9px; background:#fff; }
    .sequence-overlap-reference { display:grid; grid-template-columns:268px minmax(390px,1fr); column-gap:10px; align-items:center; min-height:28px; padding-bottom:5px; border-bottom:1px solid #edf0f5; color:#475467; }
    .sequence-overlap-row { display:grid; grid-template-columns:108px 150px minmax(390px,1fr); column-gap:10px; align-items:center; min-height:28px; }
    .sequence-overlap-row { padding:5px 0; border-bottom:1px dashed #edf0f5; }
    .sequence-overlap-row:last-child { border-bottom:0; }
    .sequence-overlap-direction { display:flex; flex-direction:column; align-items:flex-start; gap:2px; font-size:12px; white-space:nowrap; }
    .sequence-overlap-direction-label { display:inline-flex; align-items:center; gap:5px; font-weight:700; }
    .sequence-overlap-direction-label::before { content:""; width:9px; height:9px; border-radius:3px; background:#dc2626; }
    .sequence-overlap-direction.lower .sequence-overlap-direction-label::before { background:#2563eb; }
    .sequence-overlap-fold { padding-left:14px; color:#344054; font-size:13px; font-weight:700; }
    .sequence-overlap-evidence { display:inline-flex; width:max-content; max-width:150px; padding:3px 6px; border:1px solid #d0d5dd; border-radius:8px; color:#475467; background:#f8fafc; font-size:10px; font-weight:750; line-height:1.3; white-space:normal; }
    .sequence-overlap-evidence.confirmed { color:#087b55; border-color:#a7e6c7; background:#ecfdf3; }
    .sequence-overlap-evidence.tentative { color:#a45d0c; border-color:#f5cf7c; background:#fff7e8; }
    .sequence-overlap-evidence.low { color:#9a6700; border-color:#efd99b; background:#fffaeb; }
    .sequence-overlap-evidence.mass-only { color:#667085; border-color:#d0d5dd; background:#f2f4f7; }
    .sequence-overlap-track { white-space:pre; font:12.5px/1.7 Consolas,"Courier New",monospace; letter-spacing:1.5px; color:#344054; }
    .sequence-overlap-peptide { display:inline; padding:2px 1px; border:0; border-radius:3px; background:#dc2626; color:#fff; font:inherit; letter-spacing:inherit; cursor:pointer; }
    .sequence-overlap-peptide.lower { background:#2563eb; }
    .sequence-overlap-modified-site { padding:1px 0; border-radius:2px; background:#facc15; color:#111827; font-weight:800; box-shadow:0 0 0 1px #d97706 inset; }
    .sequence-overlap-modifications { display:block; min-height:18px; padding:3px 0 1px; color:#92400e; font:11.5px/1.4 Arial,sans-serif; letter-spacing:0; white-space:normal; }
    .sequence-overlap-modification { display:inline-block; margin-right:6px; padding:1px 5px; border-radius:8px; background:#fef3c7; border:1px solid #f59e0b; }
    #structureViewer { position:relative; width:100%; height:500px; overflow:hidden; border:1px solid #dbe6f2; border-radius:10px; background:#fff; }
    .structure-source-link { color:#175cd3; font-size:12px; }
    @media (max-width:1100px) { .structure-grid { grid-template-columns:1fr; } #structureViewer { height:430px; } }
    @media (max-width:650px) { header { display:flex; flex-wrap:wrap; gap:10px; padding:13px 16px; } header h1 { flex:1 1 180px; font-size:15px; } main { display:grid; gap:15px; padding:22px 11px 42px; } section { padding:17px 14px; border-radius:14px; } .controls { align-items:flex-start; flex-direction:column; } .controls label { width:100%; justify-content:space-between; } .controls select,.controls input:not([type="checkbox"]) { max-width:100%; } .pan { width:100%; } .xic-grid { grid-template-columns:1fr; } .settings-layout { grid-template-columns:1fr; } .settings-layout .settings-group:last-child { grid-column:auto; } .settings-column-grid { grid-template-columns:repeat(2,minmax(0,1fr)); } .report-settings summary { align-items:flex-start; flex-direction:column; gap:4px; } .settings-footer { align-items:stretch; flex-direction:column; } .settings-footer .note { margin-right:0; } .structure-toolbar label { width:100%; } .structure-toolbar input[type="file"] { max-width:none; width:100%; } #structureViewer { height:360px; } }
  </style>
  <script src="/assets/3Dmol-min.js"></script>
</head>
<body>
<header>
  <h1 id="reportTitle">分析报告</h1>
  <a id="homeLink" href="/">返回主页</a>
</header>
<main>
  <section class="wide report-settings" id="reportSettings">
    <details id="settingsDetails">
      <summary><span>报告显示设置</span><span class="note">模块、表格列、图表颜色与 XIC 设置会自动保存</span></summary>
      <p class="note settings-intro">设置仅影响当前报告在本浏览器中的显示，不会改变已经完成的分析结果。</p>
      <div class="settings-layout">
        <div class="settings-group">
          <h3>报告模块</h3>
          <div id="moduleSettings" class="settings-options"></div>
        </div>
        <div class="settings-group">
          <h3>图表与 XIC</h3>
          <div id="chartSettings" class="settings-options"></div>
          <div id="sampleColorSettings" class="settings-color-list"></div>
        </div>
        <div class="settings-group">
          <h3>表格列</h3>
          <div id="tableSettings"></div>
        </div>
      </div>
      <div class="settings-footer"><span id="settingsStatus" class="note">修改后会立即应用</span><button id="resetReportSettings" class="settings-reset" type="button">恢复默认设置</button></div>
    </details>
  </section>

  <section class="wide" id="chromModule">
    <h2>1. Aligned TIC/BPC overlay with detected TIC peaks</h2>
    <div class="controls">
      <label>Signal <select id="chromMode"><option value="tic">TIC</option><option value="bpc">BPC</option></select></label>
      <label><input type="checkbox" id="showAligned" checked> aligned RT</label>
      <label><input type="checkbox" id="showLocalPeakAligned" checked> local peak RT correction</label>
      <label>Vertical offset <input id="offset" type="range" min="0" max="1.2" step="0.05" value="0"></label>
      <button id="resetChrom">Reset zoom</button>
      <label>Pan <input id="chromPan" class="pan" type="range" min="0" max="1" step="0.001" value="0" disabled></label>
    </div>
    <div id="chromLegend" class="chrom-legend" aria-label="TIC 图例"></div>
    <canvas id="chromCanvas" width="1500" height="360"></canvas>
    <div id="chromInfo" class="note"></div>
  </section>

  <section class="left" id="spectrumModule">
    <h2>2. Selected TIC peak summed spectrum</h2>
    <div class="controls">
      <button id="resetDetail">Reset spectrum zoom</button>
      <label>Pan <input id="detailPan" class="pan" type="range" min="0" max="1" step="0.001" value="0" disabled></label>
      <span id="featureInfo" class="note"></span>
    </div>
    <div id="spectrumLegend" class="chrom-legend" aria-label="MS/MS 图例"></div>
    <canvas id="detailCanvas" width="1050" height="340"></canvas>
    <div id="detailInfo" class="note"></div>
    <div id="detailTooltip" class="tooltip"></div>
  </section>

  <section class="right" id="peakTableModule">
    <h2>3. TIC Peak consistency table</h2>
    <div class="scroll tall-scroll"><table id="peakTable"></table></div>
  </section>

  <section class="wide" id="mzXicSection">
    <div class="xic-grid">
      <div id="mzModule" class="report-submodule">
        <h2>5. Top changed m/z in selected TIC peak</h2>
        <div class="scroll"><table id="mzTable"></table></div>
      </div>
      <div id="xicModule" class="report-submodule">
        <h2>6. XIC confirmation</h2>
        <div class="controls">
          <button id="resetXic">Reset XIC view</button>
          <button id="toggleXicFull">Show full XIC</button>
          <label>Pan <input id="xicPan" class="pan" type="range" min="0" max="1" step="0.001" value="0" disabled></label>
      </div>
      <div id="xicLegend" class="chrom-legend" aria-label="XIC 图例"></div>
      <canvas id="xicCanvas" width="1200" height="380"></canvas>
      <div id="xicInfo" class="note"></div>
      </div>
    </div>
  </section>

  <section class="wide" id="globalFeatureModule">
    <div class="section-title-row">
      <h2>7. Overall low-similarity feature groups</h2>
      <button id="exportGlobalFeature" type="button" class="export-button">导出 Excel</button>
    </div>
    <div class="note">MS2 identified peptide forms and strict MS1-inferred unknown components are each ranked once. MS1 inference requires neutral-mass agreement, coelution, XIC shape correlation and consistent sample trend; it is not a peptide identification. Click ▸ to inspect member Features.</div>
    <div class="scroll"><table id="globalFeatureTable"></table></div>
  </section>

  <section class="wide" id="featureMapModule">
    <h2>8. Feature-level RT-m/z heatmap</h2>
    <div class="controls">
      <label>Reference <select id="featureMapReference"></select></label>
      <label>Test <select id="featureMapTest"></select></label>
      <button id="swapFeatureMapSamples">Swap</button>
      <label>Color mode <select id="featureMapMode"><option value="direction">Pair direction</option><option value="magnitude">Fold magnitude</option></select></label>
      <button id="resetFeatureMap">Reset heatmap zoom</button>
      <label>Pan <input id="featureMapPan" class="pan" type="range" min="0" max="1" step="0.001" value="0" disabled></label>
      <span id="featureMapPairLegend" class="pair-legend note"></span>
      <span id="featureMapInfo" class="note"></span>
    </div>
    <div class="feature-analysis-grid">
      <div class="feature-analysis-pane feature-map-pane">
        <div id="featureMapLegend" class="feature-map-legend" aria-label="Heatmap color scale"></div>
        <canvas id="featureMapCanvas" class="feature-map" width="900" height="430"></canvas>
      </div>
      <div class="feature-analysis-pane">
        <div class="feature-ms2-heading">
          <h3>Selected component MS/MS spectrum</h3>
        </div>
        <div id="featureMs2Detail" class="feature-ms2-detail">Select a Feature from the heatmap or Feature list.</div>
        <canvas id="featureMs2Canvas" width="900" height="430"></canvas>
        <div id="featureMs2Info" class="note"></div>
        <div id="featureMs2Annotations" class="note" hidden></div>
      </div>
    </div>
  </section>

  <section class="wide" id="modificationQuantitationModule">
    <h2>9. 修饰水平与蛋白形式差异定量</h2>
    <div class="note">将同一肽段骨架的未修饰、已知修饰及末端加工形式归并，合并电荷态和同位素包络后计算各样本的相对 MS 响应构成。B 级进入正式相对定量，C 级仅作暂定结果；仅质量预测不进入分母。</div>
    <div class="controls">
      <label>筛选 <input id="modQuantFilter" type="text" placeholder="序列、修饰、事件或 Feature"></label>
      <label>结果 <select id="modQuantStatus"><option value="quantifiable">可定量</option><option value="formal">仅正式</option><option value="all">全部</option></select></label>
      <span class="mod-view-switch" aria-label="定量图表显示方式"><button id="modViewComposition" type="button">组成比例</button><button id="modViewXic" type="button">XIC面积</button></span>
    </div>
    <div id="modQuantSummary" class="mod-quant-summary"></div>
    <div class="scroll tall-scroll"><table id="modQuantTable" class="mod-quant-table"></table></div>
    <div class="note" style="margin-top:8px">“构成变化”表示该形式占同肽段合格形式总响应的比例变化；“形式响应倍数”表示该形式自身归一化面积的组间倍数，两者含义不同。包含关系的长短肽会先形成“酶切形态族”：最终结论采用各肽段内部修饰比例的一致性，不直接相加不同序列面积；族内合计响应只用于判断酶切分配。比例未经标准品或相对响应因子校正，不等同于申报级绝对占有率。</div>
  </section>

  <section class="wide" id="structureMappingModule">
    <h2>10. 差异组分的序列与三级结构定位</h2>
    <div class="note">全局视图把所有已获得肽段候选的差异组分定位到轻链（LC）或重链（HC）。红色和蓝色分别表示当前所选两个样本中哪一个丰度更高，具体样本名称见下方动态图例；点击有颜色的残基可反选对应Feature。三级结构可从RCSB PDB自动获取或上传PDB/mmCIF。结构定位只显示证据覆盖，不会提高原MS2鉴定置信等级。</div>
    <div class="structure-grid">
      <div class="sequence-overview-pane">
        <div class="global-sequence-overview">
          <h3>LC/HC 全序列差异覆盖总览</h3>
          <div id="globalSequenceStatus" class="structure-status note">正在载入MS2序列证据…</div>
          <div id="globalSequenceTracks" class="global-sequence-tracks"></div>
          <div id="sequenceOverlapDetails" class="sequence-overlap-details" hidden></div>
          <div class="sequence-legend note">
            <span id="sequenceHigherLegend"><i class="legend-swatch higher"></i>红色：当前第二个样本更高</span>
            <span id="sequenceLowerLegend"><i class="legend-swatch lower"></i>蓝色：当前第一个样本更高</span>
            <span><i class="legend-swatch conflict"></i>红蓝重叠：不同 Feature 方向相反（点击查看）</span>
            <span><i class="legend-swatch tentative"></i>候选证据</span>
            <span><i class="legend-swatch-site"></i>修饰位点</span>
          </div>
        </div>
      </div>
      <div class="structure-pane">
        <h3>三级结构（全部差异区域及选中Feature联动）</h3>
        <div class="structure-toolbar">
          <label>PDB ID <input id="pdbIdInput" type="text" maxlength="12" placeholder="例如 4HHB"></label>
          <button id="loadPdbStructure">从RCSB PDB获取</button>
          <label>导入结构（自动保存） <input id="structureFileInput" type="file" accept=".pdb,.ent,.cif,.mmcif,text/plain,chemical/x-pdb,chemical/x-mmcif"></label>
          <button id="saveStructureFile" class="structure-upload-button" type="button">重新保存并映射</button>
          <label>结构链 <select id="structureChainSelect"><option value="">自动匹配</option></select></label>
          <button id="resetStructureView">重置视角</button>
          <a id="structureSourceLink" class="structure-source-link" href="#" target="_blank" rel="noopener" hidden>在RCSB查看</a>
        </div>
        <div id="structureMappingStatus" class="structure-status note">请上传结构文件或输入PDB ID。</div>
        <div id="structureViewer"></div>
      </div>
    </div>
  </section>
</main>
<script>
let DATA = null;
let COMPARISONS = [];
let state = { comparison:null, selectedPeakId:null, selectedMz:null, selectedFeatureRt:null, selectedFeatureGroupId:null, selectedModificationFormKey:null, selectedModificationFormContext:null, scrollGlobalSelectionIntoView:false, chromZoom:null, detailZoom:null, xicZoom:null, xicFull:true, xic:null, xicLoading:false, xicRequestId:0, featureMapZoom:null, pairReference:null, pairTest:null, featureMapMode:"direction", globalSort:{key:"ranking_score",dir:"desc"}, mzSort:{key:"ranking_score",dir:"desc"}, expandedComponents:new Set(), expandedModificationFamilies:new Set(), expandedProteolyticFamilies:new Set(), modQuantViewMode:"composition", features:[], agentAnnotations:[], msms:null, msmsLoading:false, msmsLoaded:false, msmsError:null, msmsPromise:null, sequenceOverlapSelection:null, structureViewer:null, structureModel:null, structureText:null, structureFormat:null, structureLabel:null, structureSourceMeta:null, structureChains:[], structureSourceChainMap:{}, structureSelectedChain:"", structureRenderSignature:"", structureRequestId:0, structureLoading:false, structureAutoAttemptedComparison:null };
const colors = ["#1f77b4","#e74c3c","#2ecc71","#9b59b6","#f39c12","#00bcd4","#795548","#e91e63"];
const REPORT_MODULES = [
  {key:"chromModule", label:"TIC / BPC 图"},
  {key:"spectrumModule", label:"选中峰二级质谱图"},
  {key:"peakTableModule", label:"TIC 峰一致性表"},
  {key:"mzModule", label:"变化 m/z 表"},
  {key:"xicModule", label:"XIC 图"},
  {key:"globalFeatureModule", label:"低相似度 Feature 组"},
  {key:"featureMapModule", label:"Feature 热图与 MS/MS"},
  {key:"modificationQuantitationModule", label:"修饰水平与蛋白形式定量"},
  {key:"structureMappingModule", label:"序列与三级结构"}
];
const TABLE_COLUMN_CONFIG = {
  peakTable:[
    {key:"peak_id",label:"peak_id"},{key:"status",label:"status"},{key:"local_rt",label:"local RT"},{key:"rt_start",label:"rt_start"},{key:"rt_apex",label:"rt_apex"},{key:"rt_end",label:"rt_end"},{key:"consistency",label:"consistency"},{key:"spectrum",label:"spectrum"},{key:"chrom",label:"chrom"},{key:"presence",label:"presence"},{key:"area",label:"area"},{key:"snr",label:"S/N"}
  ],
  mzTable:[
    {key:"mz",label:"true peak m/z",sortKey:"mz"},{key:"pair_type",label:"pair type",sortKey:"pair_type"},{key:"direction",label:"test direction",sortKey:"pair_log_ratio"},{key:"abundance",label:"abundance order",sortKey:"pair_log_ratio"},{key:"confidence",label:"confidence",sortKey:"quantitation_confidence"},{key:"ranking",label:"ranking",sortKey:"ranking_score"},{key:"cv",label:"CV",sortKey:"cv"},{key:"fold",label:"pair fold",sortKey:"pair_fold"}
  ],
  globalFeatureTable:[
    {key:"rank",label:"rank",sortKey:"ranking_score"},{key:"component",label:"component / Feature",sortKey:"component_label"},{key:"tic_peak",label:"TIC peak",sortKey:"parent_tic_peak_id"},{key:"rt",label:"RT",sortKey:"representative_rt"},{key:"true_peak_mz",label:"true peak m/z",sortKey:"true_peak_mz"},{key:"envelope_mz",label:"envelope m/z",sortKey:"envelope_representative_mz"},{key:"neutral_mass",label:"neutral mass",sortKey:"component_neutral_mass"},{key:"pair_type",label:"pair type",sortKey:"pair_type"},{key:"direction",label:"test direction",sortKey:"pair_log_ratio"},{key:"abundance",label:"abundance order",sortKey:"pair_log_ratio"},{key:"pair_fold",label:"pair fold",sortKey:"pair_fold"},{key:"ln_ratio",label:"ln(test/ref)",sortKey:"pair_log_ratio"},{key:"ranking",label:"ranking",sortKey:"ranking_score"},{key:"max_area",label:"max raw area",sortKey:"max_area"},{key:"confidence",label:"confidence",sortKey:"quantitation_confidence"},{key:"ions",label:"ions",sortKey:"component_member_count"},{key:"merged",label:"merged",sortKey:"merged_feature_count"},{key:"cohort_type",label:"cohort type",sortKey:"difference_type"}
  ]
};
function defaultReportSettings(){
  return {
    modules:Object.fromEntries(REPORT_MODULES.map(item=>[item.key,true])),
    legends:{chrom:true,spectrum:true,xic:true},
    xicFull:true,
    xicIntegrationBand:true,
    chartColors:{},
    tableColumns:Object.fromEntries(Object.entries(TABLE_COLUMN_CONFIG).map(([scope,columns])=>[scope,columns.map(column=>column.key)]))
  };
}
let reportSettings = defaultReportSettings();
const MAX_FEATURE_FOLD = 1000;
const MAX_FEATURE_LN_FOLD = Math.log(MAX_FEATURE_FOLD);
const PROJECT_STRUCTURE_DEFAULTS = Object.freeze({
  vedolizumab_ql2519_vs_originator:{pdbId:"3V4P",kind:"homologous",label:"ACT-1 Fab–α4β7同源模板",note:"RCSB当前没有与vedolizumab可变区精确一致的实验结构；3V4P为其亲本抗体ACT-1与α4β7复合物，HC/LC可变区约85%一致。"},
  denosumab_boyoupei_vs_mailishu:{pdbId:"5I1C",kind:"homologous",label:"人源IGHV3-23/IGKV3-20 Fab同源模板",note:"RCSB当前没有与denosumab可变区精确一致的实验结构；5I1C是可变区较接近的Fab模板，HC/LC可变区约90%/96%一致。"}
});
const $ = id => document.getElementById(id);
function nice(v,d=3){ if(v===null||v===undefined||v==="")return ""; const n=Number(v); return Number.isFinite(n)?n.toFixed(d):""; }
function color(i){ const sample=DATA?.sample_ids?.[i]; return chartColor(sample,i); }
function sampleColor(sample){ const index=DATA?.sample_ids?.indexOf(sample)??-1; return chartColor(sample,Math.max(0,index)); }
function sampleShort(sample){ const text=String(sample||""); if(text.includes("QL2519"))return "QL2519"; if(text.toLowerCase().includes("vedolizumab"))return "VDZ"; return text.length>24?`${text.slice(0,21)}...`:text; }
function escapeHtml(value){ return String(value??"").replace(/[&<>"']/g,char=>({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"})[char]); }
async function fetchJson(url, options){ const r=await fetch(url,{cache:"no-store",...(options||{})}); if(!r.ok) throw new Error(`${url} HTTP ${r.status}`); return await r.json(); }
function reportSettingsKey(){ return `lcms_report_settings_v3_${encodeURIComponent(state.comparison||"default")}`; }
function validHexColor(value){ return /^#[0-9a-f]{6}$/i.test(String(value||"")); }
function loadReportSettings(){
  const defaults=defaultReportSettings();
  try{
    const saved=JSON.parse(localStorage.getItem(reportSettingsKey())||"null");
    if(!saved||typeof saved!=="object")return defaults;
    Object.keys(defaults.modules).forEach(key=>{ if(typeof saved.modules?.[key]==="boolean")defaults.modules[key]=saved.modules[key]; });
    Object.keys(defaults.legends).forEach(key=>{ if(typeof saved.legends?.[key]==="boolean")defaults.legends[key]=saved.legends[key]; });
    if(typeof saved.xicFull==="boolean")defaults.xicFull=saved.xicFull;
    if(typeof saved.xicIntegrationBand==="boolean")defaults.xicIntegrationBand=saved.xicIntegrationBand;
    Object.entries(saved.chartColors||{}).forEach(([sample,value])=>{ if(validHexColor(value))defaults.chartColors[sample]=value; });
    Object.entries(defaults.tableColumns).forEach(([scope,columns])=>{
      const savedColumns=Array.isArray(saved.tableColumns?.[scope])?saved.tableColumns[scope].filter(key=>columns.includes(key)):[];
      if(savedColumns.length)defaults.tableColumns[scope]=savedColumns;
    });
  }catch(_error){}
  return defaults;
}
function saveReportSettings(message="已自动保存"){
  try{ localStorage.setItem(reportSettingsKey(),JSON.stringify(reportSettings)); }catch(_error){}
  const status=$("settingsStatus"); if(status)status.textContent=message;
  applyReportSettings();
}
function chartColor(sample,index=0){ return reportSettings.chartColors?.[sample]||colors[index%colors.length]; }
function tableColumns(scope){
  const defaults=(TABLE_COLUMN_CONFIG[scope]||[]).map(column=>column.key);
  const selected=reportSettings.tableColumns?.[scope];
  return defaults.filter(key=>Array.isArray(selected)&&selected.includes(key));
}
function tableColumnVisible(scope,key){ return tableColumns(scope).includes(key); }
function applyTableColumnVisibility(scope,tableId){
  const table=$(tableId); if(!table)return;
  const columns=(TABLE_COLUMN_CONFIG[scope]||[]).map(column=>column.key);
  table.querySelectorAll("tr").forEach(row=>Array.from(row.children).forEach((cell,index)=>{ if(columns[index])cell.hidden=!tableColumnVisible(scope,columns[index]); }));
}
function renderReportSettings(){
  const moduleContainer=$("moduleSettings"), chartContainer=$("chartSettings"), colorContainer=$("sampleColorSettings"), tableContainer=$("tableSettings");
  if(!moduleContainer||!chartContainer||!colorContainer||!tableContainer)return;
  moduleContainer.innerHTML=REPORT_MODULES.map(item=>`<label class="settings-check"><input type="checkbox" data-setting-module="${item.key}" ${reportSettings.modules[item.key]?"checked":""}><span>${item.label}</span></label>`).join("");
  chartContainer.innerHTML=`<label class="settings-check"><input type="checkbox" data-setting-legend="chrom" ${reportSettings.legends.chrom?"checked":""}><span>TIC / BPC 图例</span></label><label class="settings-check"><input type="checkbox" data-setting-legend="spectrum" ${reportSettings.legends.spectrum?"checked":""}><span>二级质谱图例</span></label><label class="settings-check"><input type="checkbox" data-setting-legend="xic" ${reportSettings.legends.xic?"checked":""}><span>XIC 图例</span></label><label class="settings-select"><span>默认 XIC</span><select id="settingXicMode"><option value="full">整体 XIC</option><option value="local">局部 XIC</option></select></label><label class="settings-check"><input id="settingXicBand" type="checkbox" ${reportSettings.xicIntegrationBand?"checked":""}><span>显示 XIC 积分区间</span></label>`;
  $("settingXicMode").value=reportSettings.xicFull?"full":"local";
  const samples=DATA?.sample_ids||[];
  colorContainer.innerHTML=samples.length?`<div class="settings-table-title">样品颜色</div>${samples.map((sample,index)=>`<label class="settings-color-row"><input type="color" data-setting-color="${escapeHtml(sample)}" value="${chartColor(sample,index)}"><span title="${escapeHtml(sample)}">${escapeHtml(sample)}</span></label>`).join("")}`:`<div class="note">载入样品后可设置颜色</div>`;
  const tableLabels={peakTable:"TIC 峰一致性表",mzTable:"变化 m/z 表",globalFeatureTable:"Feature 组表"};
  tableContainer.innerHTML=Object.entries(TABLE_COLUMN_CONFIG).map(([scope,columns])=>`<div class="settings-table-group"><div class="settings-table-title">${tableLabels[scope]||scope}</div><div class="settings-column-grid">${columns.map(column=>`<label title="${escapeHtml(column.label)}"><input type="checkbox" data-setting-column-scope="${scope}" data-setting-column="${column.key}" ${tableColumnVisible(scope,column.key)?"checked":""}><span>${escapeHtml(column.label)}</span></label>`).join("")}</div></div>`).join("");
}
function bindReportSettings(){
  document.querySelectorAll("[data-setting-module]").forEach(input=>input.onchange=()=>{ reportSettings.modules[input.dataset.settingModule]=input.checked; saveReportSettings(); });
  document.querySelectorAll("[data-setting-legend]").forEach(input=>input.onchange=()=>{ reportSettings.legends[input.dataset.settingLegend]=input.checked; saveReportSettings(); drawAll(); });
  document.querySelectorAll("[data-setting-color]").forEach(input=>input.oninput=()=>{ reportSettings.chartColors[input.dataset.settingColor]=input.value; saveReportSettings(); drawAll(); });
  document.querySelectorAll("[data-setting-column]").forEach(input=>input.onchange=()=>{
    const scope=input.dataset.settingColumnScope, key=input.dataset.settingColumn, current=tableColumns(scope);
    if(!input.checked&&current.length<=1){ input.checked=true; const status=$("settingsStatus"); if(status)status.textContent="每张表至少保留一列"; return; }
    reportSettings.tableColumns[scope]=input.checked?[...current.filter(item=>item!==key),key]:(current.filter(item=>item!==key));
    saveReportSettings();
    const tableId=scope==="peakTable"?"peakTable":scope==="mzTable"?"mzTable":"globalFeatureTable";
    applyTableColumnVisibility(scope,tableId);
  });
  $("settingXicMode").onchange=async()=>{
    state.xicFull=$("settingXicMode").value==="full"; reportSettings.xicFull=state.xicFull; saveReportSettings(); updateXicToggleLabel();
    if(state.selectedMz)await selectSpectrumMz(state.selectedMz,state.selectedFeatureRt); else drawAll();
  };
  $("settingXicBand").onchange=()=>{ reportSettings.xicIntegrationBand=$("settingXicBand").checked; saveReportSettings(); drawAll(); };
  $("resetReportSettings").onclick=()=>{ reportSettings=defaultReportSettings(); state.xicFull=true; saveReportSettings("已恢复默认设置"); renderReportSettings(); bindReportSettings(); drawAll(); };
}
function setupReportSettings(){ renderReportSettings(); bindReportSettings(); applyReportSettings(); }
function applyReportSettings(){
  REPORT_MODULES.forEach(item=>{ const element=$(item.key); if(element)element.hidden=reportSettings.modules[item.key]===false; });
  const parent=$("mzXicSection"), mzVisible=reportSettings.modules.mzModule!==false, xicVisible=reportSettings.modules.xicModule!==false;
  if(parent){ parent.hidden=!mzVisible&&!xicVisible; parent.classList.toggle("single-module",mzVisible!==xicVisible); }
  const chromLegend=$("chromLegend"), spectrumLegend=$("spectrumLegend"), xicLegend=$("xicLegend");
  if(chromLegend)chromLegend.hidden=!reportSettings.legends.chrom;
  if(spectrumLegend)spectrumLegend.hidden=!reportSettings.legends.spectrum;
  if(xicLegend)xicLegend.hidden=!reportSettings.legends.xic;
  updateXicToggleLabel();
}
function updateXicToggleLabel(){ const button=$("toggleXicFull"); if(button)button.textContent=state.xicFull?"显示局部 XIC":"显示整体 XIC"; }
function apiUrl(path){ const sep=path.includes("?")?"&":"?"; return `${path}${sep}comparison=${encodeURIComponent(state.comparison||"")}`; }
function selectedPair(){
  const samples=DATA?.sample_ids||[];
  let reference=samples.includes(state.pairReference)?state.pairReference:(DATA?.reference_sample||samples[0]);
  let test=samples.includes(state.pairTest)&&state.pairTest!==reference?state.pairTest:samples.find(sample=>sample!==reference);
  return {reference,test};
}
function featureValueMap(item){ return item?.normalized_area_by_sample||item?.normalized_intensity_by_sample||item?.area_by_sample||item?.raw_intensity_by_sample||{}; }
function truePeakMz(item){ const value=Number(item?.true_peak_mz??item?.representative_mz??item?.mz); return Number.isFinite(value)?value:null; }
function envelopeRepresentativeMz(item){ const value=Number(item?.component_representative_mz??item?.envelope_representative_mz??item?.component_observed_first_isotope_mz??truePeakMz(item)); return Number.isFinite(value)?value:null; }
function featureGroupIds(item){
  const values=[item?.feature_group_id,item?.component_primary_feature_id,...(item?.source_feature_group_ids||[]),...(item?.merged_feature_group_ids||[])];
  return [...new Set(values.filter(Boolean).map(String))];
}
function componentContainsFeature(item,featureGroupId){
  const target=String(featureGroupId||""); if(!target)return false;
  return featureGroupIds(item).includes(target)||(item?.members||[]).some(member=>featureGroupIds(member).includes(target));
}
function currentComponentRows(){
  const updated=state.msms?.ms1_component_groups;
  return Array.isArray(updated)&&updated.length?updated:(DATA?.ms1_component_groups||[]);
}
const ms2StatusText={
  identified_direct_precursor:"direct precursor identified",
  identified_isotope_envelope:"identified through isotope envelope",
  identified_charge_state_envelope:"identified through charge-state envelope",
  tentative_component_consensus_identification:"component charge/isotope joint b/y candidate",
  tentative_backbone_sequence_support:"mass-offset backbone sequence candidate",
  tentative_truncation_sequence_support:"truncation backbone sequence candidate",
  tentative_sequence_region_candidate:"partial b/y sequence-region candidate",
  tentative_feature_consensus_identification:"single-Feature consensus candidate",
  tentative_feature_guided_identification:"MS1-guided tentative identification",
  tentative_feature_open_mass_identification:"single-Feature open-mass candidate",
  tentative_feature_glycopeptide_identification:"targeted glycopeptide candidate",
  tentative_global_open_modification_identification:"global open-modification search",
  low_evidence_sequence_candidate:"low-evidence sequence candidate",
  selected_precursor_unidentified:"selected precursor, unidentified",
  coisolated_ms2_unresolved:"isolation window covered, unresolved",
  no_ms2_acquired:"no MS2 scan acquired"
};
function selectedComponentFeatureIds(){
  const target=String(state.selectedFeatureGroupId||"");
  const ids=new Set(target?[target]:[]);
  const component=currentComponentRows().find(item=>componentContainsFeature(item,target));
  if(component){
    featureGroupIds(component).forEach(id=>ids.add(id));
    (component.members||[]).forEach(member=>featureGroupIds(member).forEach(id=>ids.add(id)));
  }
  return ids;
}
function selectedMs2Evidence(){
  const target=String(state.selectedFeatureGroupId||"");
  const evidence=state.msms?.feature_evidence||[];
  const exact=evidence.find(item=>String(item.feature_group_id||"")===target);
  if(exact)return {feature:exact,relation:"exact Feature"};
  const componentIds=selectedComponentFeatureIds();
  const related=evidence.filter(item=>componentIds.has(String(item.feature_group_id||"")));
  related.sort((a,b)=>Number(Boolean(b.best_psm))-Number(Boolean(a.best_psm)) || Number(a.rank||1e9)-Number(b.rank||1e9));
  return related.length?{feature:related[0],relation:"linked component member"}:null;
}
function componentCandidateEvidence(item){
  if(item?.component_candidate_sequence)return {sequence:item.component_candidate_sequence,modification:item.component_candidate_modification,confidence:item.component_candidate_confidence||"D_low_evidence"};
  const ids=new Set(featureGroupIds(item));
  (item?.members||[]).forEach(member=>featureGroupIds(member).forEach(id=>ids.add(id)));
  const candidates=(state.msms?.feature_evidence||[]).filter(evidence=>ids.has(String(evidence.feature_group_id||""))&&evidence.sequence&&String(evidence.confidence||"")==="D_low_evidence");
  candidates.sort((a,b)=>Number(b.best_psm?.score||b.candidate_psm?.score||0)-Number(a.best_psm?.score||a.candidate_psm?.score||0));
  return candidates[0]||null;
}
function preferredComponentEvidenceFeatureId(item){
  const ids=new Set(featureGroupIds(item));
  (item?.members||[]).forEach(member=>featureGroupIds(member).forEach(id=>ids.add(id)));
  const statusPriority={
    identified_direct_precursor:900,
    tentative_feature_glycopeptide_identification:850,
    identified_isotope_envelope:800,
    identified_charge_state_envelope:790,
    tentative_feature_open_mass_identification:760,
    tentative_feature_guided_identification:750,
    tentative_feature_consensus_identification:740,
    tentative_component_consensus_identification:730,
    tentative_backbone_sequence_support:720,
    tentative_truncation_sequence_support:710,
    tentative_sequence_region_candidate:700,
    low_evidence_sequence_candidate:600,
    selected_precursor_unidentified:200,
    coisolated_ms2_unresolved:100,
    no_ms2_acquired:0,
  };
  const candidates=(state.msms?.feature_evidence||[]).filter(evidence=>ids.has(String(evidence.feature_group_id||"")));
  candidates.sort((a,b)=>{
    const aPriority=Number(statusPriority[String(a.ms2_status||"")]||0), bPriority=Number(statusPriority[String(b.ms2_status||"")]||0);
    return bPriority-aPriority || Number(b.best_psm?.score||b.candidate_psm?.score||0)-Number(a.best_psm?.score||a.candidate_psm?.score||0) || Number(a.rank||1e9)-Number(b.rank||1e9);
  });
  return String(candidates[0]?.feature_group_id||item?.component_primary_feature_id||item?.feature_group_id||"");
}
function componentDisplayLabel(item){
  const base=String(item?.component_label||item?.feature_group_id||"");
  if(item?.identified_component)return base;
  const candidate=componentCandidateEvidence(item), sequence=String(candidate?.sequence||"");
  return sequence&&!base.includes(sequence)?`${base} | candidate ${sequence}`:base;
}
const residueCode={
  ALA:"A",ARG:"R",ASN:"N",ASP:"D",CYS:"C",GLN:"Q",GLU:"E",GLY:"G",HIS:"H",ILE:"I",
  LEU:"L",LYS:"K",MET:"M",PHE:"F",PRO:"P",SER:"S",THR:"T",TRP:"W",TYR:"Y",VAL:"V",
  MSE:"M",SEC:"U",PYL:"O"
};
function cleanPeptideSequence(value){ return String(value||"").toUpperCase().replace(/[^A-Z]/g,""); }
function featureItemById(featureId){
  const target=String(featureId||"");
  if(!target)return null;
  const direct=(DATA?.global_feature_groups||[]).find(item=>featureGroupIds(item).includes(target));
  if(direct)return direct;
  for(const component of currentComponentRows()){
    const member=(component.members||[]).find(item=>featureGroupIds(item).includes(target));
    if(member)return member;
  }
  return null;
}
function sequenceLocationForEvidence(f,relation=""){
  if(!f||typeof f!=="object")return null;
  const candidates=[f.best_psm,f.candidate_psm,...(f.sequence_region_candidates||[]),...(f.candidate_mass_hypotheses||[])].filter(item=>item&&typeof item==="object");
  const requested=cleanPeptideSequence(f.sequence);
  const source=candidates.find(item=>cleanPeptideSequence(item.sequence)===requested&&item.chain) || candidates.find(item=>item.sequence&&item.chain) || candidates.find(item=>item.sequence) || {};
  const sequence=requested||cleanPeptideSequence(source.sequence);
  if(!sequence)return null;
  let chain=String(source.chain||f.chain||"");
  const chains=state.msms?.chains||{};
  if(!cleanPeptideSequence(chains?.[chain])){
    const exactChains=Object.entries(chains).filter(([,chainSequence])=>cleanPeptideSequence(chainSequence).includes(sequence));
    if(exactChains.length===1)chain=String(exactChains[0][0]);
  }
  let start=Number(source.start??f.start), end=Number(source.end??f.end);
  const chainSequence=cleanPeptideSequence(chains?.[chain]);
  if((!Number.isFinite(start)||!Number.isFinite(end)||end<start)&&chainSequence){
    const index=chainSequence.indexOf(sequence);
    if(index>=0){ start=index+1; end=index+sequence.length; }
  }
  return {
    featureId:String(f.feature_group_id||""),
    featureIds:[String(f.feature_group_id||"")].filter(Boolean),
    evidenceFeatureId:String(f.feature_group_id||""),
    sequence,
    chain,
    start:Number.isFinite(start)?start:null,
    end:Number.isFinite(end)?end:null,
    modification:String(f.modification||source.modification_text||source.modification||""),
    confidence:String(f.confidence||source.sequence_inference_level||source.q_value||""),
    status:String(f.ms2_status||""),
    relation,
    evidence:f,
    source
  };
}
function selectedSequenceLocation(){
  const context=state.selectedModificationFormContext;
  if(context){
    const pair=selectedPair(),areas=context.normalized_area_by_sample||{},referenceValue=Number(areas[pair.reference]||0),testValue=Number(areas[pair.test]||0);
    return {featureId:String(state.selectedFeatureGroupId||""),featureIds:[String(state.selectedFeatureGroupId||"")].filter(Boolean),evidenceFeatureId:String(state.selectedFeatureGroupId||""),sequence:cleanPeptideSequence(context.sequence),chain:String(context.chain||""),start:Number(context.start),end:Number(context.end),modification:String(context.form_label||context.modification||""),confidence:String(context.ms2_confidence||""),status:String(context.status||""),relation:"第9节选中修饰形式",evidence:context,source:context,direction:testValue>referenceValue?"higher":(testValue<referenceValue?"lower":"equal"),fold:Math.max(referenceValue,testValue)/Math.max(Math.min(referenceValue,testValue),1e-12),logRatio:Math.log(Math.max(testValue,1e-12)/Math.max(referenceValue,1e-12))};
  }
  const matched=selectedMs2Evidence();
  if(!matched)return null;
  const location=sequenceLocationForEvidence(matched.feature,matched.relation);
  if(!location)return null;
  const feature=featureItemById(location.featureId)||matched.feature, pair=pairMetrics(feature);
  location.direction=pair.direction; location.fold=pair.fold; location.logRatio=pair.logRatio;
  return location;
}
function modificationPositions(location){
  const positions=new Set();
  if(!location||!Number.isFinite(location.start))return positions;
  for(const match of String(location.modification||"").matchAll(/@(\d+)/g)){
    const relative=Number(match[1]);
    if(Number.isFinite(relative))positions.add(location.start+relative-1);
  }
  return positions;
}
function evidenceIsTentative(location){
  const confidence=String(location?.confidence||"").toLowerCase(), status=String(location?.status||"").toLowerCase();
  return confidence.startsWith("c_")||confidence.startsWith("d_")||confidence.includes("tentative")||status.startsWith("tentative_")||status==="low_evidence_sequence_candidate";
}
function sequenceEvidenceLabel(location){
  const status=String(location?.status||"").toLowerCase(), confidence=String(location?.confidence||"").toLowerCase();
  if(["selected_precursor_unidentified","coisolated_ms2_unresolved","no_ms2_acquired"].includes(status)){
    return {kind:"mass-only",label:"仅质量预测",detail:"无可用于序列定位的碎片证据；序列、修饰和位点均未确认。"};
  }
  if(status==="identified_direct_precursor"){
    return {kind:"confirmed",label:"B级：直接 MS2 定性",detail:"前体与碎片离子支持该序列。"};
  }
  if(["identified_isotope_envelope","identified_charge_state_envelope"].includes(status)){
    return {kind:"confirmed",label:"B级：MS2 包络关联",detail:"由已鉴定 MS2 前体的同位素或电荷态包络关联。"};
  }
  if(status==="low_evidence_sequence_candidate"||confidence.startsWith("d_")){
    return {kind:"low",label:"D级：低证据候选",detail:"存在有限的 MS2 线索，但不足以支持可靠定性。"};
  }
  if(status==="tentative_feature_open_mass_identification"){
    return {kind:"tentative",label:"C级：开放质量候选",detail:"有 MS2 骨架证据，但质量偏移的修饰类型或位点仍待确认。"};
  }
  if(status==="tentative_feature_glycopeptide_identification"){
    return {kind:"tentative",label:"C级：定向糖肽候选",detail:"有糖链诊断离子、核心 Y 离子及肽骨架证据，仍需人工或正交方法确认。"};
  }
  if(status.startsWith("tentative_")||confidence.startsWith("c_")){
    return {kind:"tentative",label:"C级：MS2 暂定候选",detail:"有部分 MS2 支持，但未达到高可信定性标准。"};
  }
  return {kind:"low",label:"证据等级未定",detail:"请结合对应 Feature 的 MS2 状态和谱图复核。"};
}
function updateSequenceDirectionLegend(){
  const pair=selectedPair(), higher=$("sequenceHigherLegend"), lower=$("sequenceLowerLegend");
  const reference=sampleShort(pair.reference||"样本1"), test=sampleShort(pair.test||"样本2");
  if(higher){higher.innerHTML=`<i class="legend-swatch higher"></i>红色：${escapeHtml(test)} 高于 ${escapeHtml(reference)}`;higher.title=`${pair.test||test} > ${pair.reference||reference}`;}
  if(lower){lower.innerHTML=`<i class="legend-swatch lower"></i>蓝色：${escapeHtml(reference)} 高于 ${escapeHtml(test)}`;lower.title=`${pair.reference||reference} > ${pair.test||test}`;}
}
function globalSequenceLocations(){
  const groups=new Map();
  for(const evidence of state.msms?.feature_evidence||[]){
    const location=sequenceLocationForEvidence(evidence,"global differential evidence");
    if(!location||!location.chain||!Number.isFinite(location.start)||!Number.isFinite(location.end))continue;
    const feature=featureItemById(location.featureId)||evidence;
    const differenceType=pairDifferenceType(feature);
    if(["common_feature","low_confidence"].includes(differenceType))continue;
    const pair=pairMetrics(feature);
    if(pair.direction==="equal")continue;
    location.direction=pair.direction;
    location.fold=pair.fold;
    location.logRatio=pair.logRatio;
    location.differenceType=differenceType;
    location.tentative=evidenceIsTentative(location);
    location.abundance=Math.max(pair.referenceValue,pair.testValue);
    location.modifiedPositions=modificationPositions(location);
    const key=[location.chain,location.start,location.end,location.sequence,location.modification,location.direction].join("|");
    const existing=groups.get(key);
    if(existing){
      existing.featureIds=[...new Set([...existing.featureIds,...location.featureIds])];
      existing.fold=Math.max(existing.fold,location.fold);
      existing.abundance=Math.max(existing.abundance,location.abundance);
      existing.tentative=existing.tentative&&location.tentative;
    }else groups.set(key,location);
  }
  return [...groups.values()].sort((a,b)=>String(a.chain).localeCompare(String(b.chain))||a.start-b.start||a.end-b.end);
}
function chainDisplayName(chain){
  const upper=String(chain||"").toUpperCase();
  if(upper.includes("_HC")||upper.endsWith("HC"))return `HC（重链） · ${chain}`;
  if(upper.includes("_LC")||upper.endsWith("LC"))return `LC（轻链） · ${chain}`;
  return chain||"未命名链";
}
async function selectMappedFeature(featureId,preserveModificationContext=false){
  const item=featureItemById(featureId);
  if(!item)return;
  await selectGlobalFeature({...item,feature_group_id:String(featureId)},preserveModificationContext);
}
function residueHasDirectionOverlap(hits){
  const directions=new Set((hits||[]).map(item=>item.direction));
  return directions.has("higher")&&directions.has("lower");
}
function sequenceEvidenceSegment(locations,chains,selection){
  if(!selection)return null;
  const chainEntry=chains.find(([chain])=>chain===selection.chain);
  if(!chainEntry)return null;
  const [chain,sequence]=chainEntry, chainLocations=locations.filter(item=>item.chain===chain);
  const coverage=Array.from({length:sequence.length},()=>[]);
  chainLocations.forEach(location=>{
    for(let position=Math.max(1,location.start);position<=Math.min(sequence.length,location.end);position++)coverage[position-1].push(location);
  });
  const selectedPosition=Math.max(1,Math.min(sequence.length,Number(selection.position)||1));
  const selectedHits=coverage[selectedPosition-1]||[];
  if(!selectedHits.length)return null;
  let start=selectedPosition, end=selectedPosition;
  const hasDirectionConflict=residueHasDirectionOverlap(selectedHits);
  if(hasDirectionConflict){
    while(start>1&&residueHasDirectionOverlap(coverage[start-2]))start-=1;
    while(end<sequence.length&&residueHasDirectionOverlap(coverage[end]))end+=1;
  }else{
    while(start>1&&coverage[start-2].length&&!residueHasDirectionOverlap(coverage[start-2]))start-=1;
    while(end<sequence.length&&coverage[end].length&&!residueHasDirectionOverlap(coverage[end]))end+=1;
  }
  return {chain,sequence,chainLocations,coverage,selectedPosition,start,end,hasDirectionConflict};
}
function renderSequenceOverlapDetails(locations,chains){
  const panel=$("sequenceOverlapDetails"), selection=state.sequenceOverlapSelection;
  if(!panel)return;
  if(!selection){panel.hidden=true;panel.innerHTML="";return;}
  const segment=sequenceEvidenceSegment(locations,chains,selection);
  if(!segment){state.sequenceOverlapSelection=null;panel.hidden=true;panel.innerHTML="";return;}
  const {chain,sequence,chainLocations,start:segmentStart,end:segmentEnd,hasDirectionConflict}=segment;
  const viewStart=Math.max(1,segmentStart-5), viewEnd=Math.min(sequence.length,segmentEnd+5), pair=selectedPair();
  const relevant=chainLocations.filter(item=>item.start<=segmentEnd&&item.end>=segmentStart).sort((a,b)=>a.direction.localeCompare(b.direction)||a.start-b.start||b.fold-a.fold);
  const referenceSequence=sequence.slice(viewStart-1,viewEnd), referenceTrack=`${viewStart}  ${referenceSequence}  ${viewEnd}`;
  const rows=relevant.map(item=>{
    const clippedStart=Math.max(viewStart,item.start), clippedEnd=Math.min(viewEnd,item.end);
    const prefix=" ".repeat(clippedStart-viewStart+String(viewStart).length+2);
    const lower=item.direction==="lower", higherSample=lower?pair.reference:pair.test;
    const featureId=item.featureIds?.[0]||"", modification=String(item.modification||"");
    const evidenceLabel=sequenceEvidenceLabel(item);
    const annotations=[...modification.matchAll(/(?:^|;\s*)([^;@]+)@(\d+)/g)].map(match=>({name:String(match[1]||"").trim(),relative:Number(match[2]),absolute:item.start+Number(match[2])-1})).filter(annotation=>annotation.name&&Number.isFinite(annotation.relative));
    const peptideHtml=sequence.slice(clippedStart-1,clippedEnd).split("").map((aa,index)=>{
      const absolute=clippedStart+index, labels=annotations.filter(annotation=>annotation.absolute===absolute).map(annotation=>annotation.name);
      return labels.length?`<span class="sequence-overlap-modified-site" title="${escapeHtml(labels.join("; "))} @ ${absolute}">${escapeHtml(aa)}</span>`:escapeHtml(aa);
    }).join("");
    const namedModification=modification&&modification!=="Unmodified";
    const annotationHtml=annotations.length?annotations.map(annotation=>`<span class="sequence-overlap-modification">${escapeHtml(annotation.name)}：${annotation.absolute}</span>`).join(""):(namedModification?`<span class="sequence-overlap-modification">${escapeHtml(modification)}（位点未确定）</span>`:"");
    const title=`点击选择 ${(item.featureIds||[]).join(", ")}${namedModification?` | ${modification}`:""} | ${evidenceLabel.label}`;
    return `<div class="sequence-overlap-row"><span class="sequence-overlap-direction ${lower?"lower":"higher"}"><span class="sequence-overlap-direction-label">${escapeHtml(sampleShort(higherSample))} 较高</span><span class="sequence-overlap-fold">${nice(item.fold,2)}x</span></span><span class="sequence-overlap-evidence ${escapeHtml(evidenceLabel.kind)}" title="${escapeHtml(evidenceLabel.detail)}">${escapeHtml(evidenceLabel.label)}</span><span class="sequence-overlap-track">${prefix}<button class="sequence-overlap-peptide ${lower?"lower":"higher"}" data-overlap-feature-id="${escapeHtml(featureId)}" title="${escapeHtml(title)}">${peptideHtml}</button>${annotationHtml?`<span class="sequence-overlap-modifications">${annotationHtml}</span>`:""}</span></div>`;
  }).join("");
  panel.hidden=false;
  const title=hasDirectionConflict?"重叠区（不同 Feature 方向相反）":"选中区域（无方向冲突）";
  panel.innerHTML=`<div class="sequence-overlap-heading"><b>${title}：${escapeHtml(chainDisplayName(chain))} ${segmentStart}-${segmentEnd}</b><button id="closeSequenceOverlap" type="button">关闭</button></div><div class="sequence-overlap-note">按一级序列真实位置对齐；红色表示 ${escapeHtml(sampleShort(pair.test))} 较高，蓝色表示 ${escapeHtml(sampleShort(pair.reference))} 较高，黄色表示修饰位点。中间的证据标签说明该肽段属于 B/C/D 级，或仅质量预测。点击肽段可反选对应 Feature 和 MS/MS。</div><div class="sequence-overlap-scroll"><div class="sequence-overlap-reference"><b>原序列</b><span class="sequence-overlap-track">${escapeHtml(referenceTrack)}</span></div>${rows}</div>`;
  $("closeSequenceOverlap").onclick=()=>{state.sequenceOverlapSelection=null;renderGlobalSequenceOverview(locations);};
  panel.querySelectorAll("[data-overlap-feature-id]").forEach(element=>element.onclick=()=>selectMappedFeature(element.dataset.overlapFeatureId));
}
function renderGlobalSequenceOverview(locations){
  const status=$("globalSequenceStatus"), container=$("globalSequenceTracks");
  if(!status||!container)return;
  updateSequenceDirectionLegend();
  if(state.msmsLoading||!state.msmsLoaded){
    status.innerHTML='<span class="structure-badge">载入MS2证据</span>正在整理全部差异肽段的LC/HC位置。';
    container.innerHTML="";
    return;
  }
  if(state.msmsError){
    const skipped=String(DATA?.analysis_metadata?.ms2?.status||"").startsWith("skipped_");
    status.innerHTML=`<span class="structure-badge missing">${skipped?"MS2未执行":"MS2载入失败"}</span>${escapeHtml(state.msmsError)}`;
    container.innerHTML="";
    return;
  }
  const chains=Object.entries(state.msms?.chains||{}).map(([name,sequence])=>[String(name),cleanPeptideSequence(sequence)]).filter(([,sequence])=>sequence);
  if(!chains.length){
    status.innerHTML='<span class="structure-badge missing">缺少完整序列</span>MS2报告中没有可用的LC/HC FASTA序列，无法绘制全序列定位。';
    container.innerHTML="";
    return;
  }
  const selectedIds=selectedComponentFeatureIds(), pair=selectedPair(), selectedOverlap=sequenceEvidenceSegment(locations,chains,state.sequenceOverlapSelection);
  let coveredResidues=0;
  const cards=chains.map(([chain,sequence])=>{
    const chainLocations=locations.filter(item=>item.chain===chain);
    const coverage=Array.from({length:sequence.length},()=>[]);
    chainLocations.forEach(location=>{
      const from=Math.max(1,location.start), to=Math.min(sequence.length,location.end);
      for(let position=from;position<=to;position++)coverage[position-1].push(location);
    });
    const chainCovered=coverage.filter(items=>items.length).length;
    coveredResidues+=chainCovered;
    const blocks=[];
    for(let offset=0;offset<sequence.length;offset+=50){
      const residues=sequence.slice(offset,offset+50).split("").map((aa,index)=>{
        const position=offset+index+1, hits=coverage[position-1]||[];
        const classes=["sequence-residue"], directions=new Set(hits.map(item=>item.direction)), conflict=residueHasDirectionOverlap(hits);
        if(hits.length){
          classes.push("global-difference",conflict?"diff-conflict":(directions.has("higher")?"diff-higher":"diff-lower"));
          if(hits.every(item=>item.tentative))classes.push("diff-tentative");
          if(hits.some(item=>item.featureIds.some(id=>selectedIds.has(id))))classes.push("diff-selected");
          if(hits.some(item=>item.modifiedPositions?.has(position)))classes.push("modified");
          if(selectedOverlap?.chain===chain&&position>=selectedOverlap.start&&position<=selectedOverlap.end){
            classes.push("diff-overlap-selected");
            if(position===selectedOverlap.start)classes.push("diff-overlap-start");
            if(position===selectedOverlap.end)classes.push("diff-overlap-end");
          }
        }
        const details=hits.slice(0,5).map(item=>`${item.direction==="higher"?"↑":"↓"} ${item.direction==="higher"?`${sampleShort(pair.test)} > ${sampleShort(pair.reference)}`:`${sampleShort(pair.reference)} > ${sampleShort(pair.test)}`} ${nice(item.fold,2)}x | ${item.featureIds.join(", ")} | ${item.sequence}${item.modification&&item.modification!=="Unmodified"?` | ${item.modification}`:""} | ${item.confidence||item.status}`).join("\n");
        const title=`${chain} ${position} ${aa}${details?`\n${details}`:""}`;
        const data=hits.length?` data-overlap-chain="${escapeHtml(chain)}" data-overlap-position="${position}"`:"";
        return `<span class="${classes.join(" ")}"${data} title="${escapeHtml(title)}">${aa}</span>`;
      }).join("");
      blocks.push(`<div class="sequence-block"><span class="sequence-position">${offset+1}</span><span class="sequence-residues">${residues}</span></div>`);
    }
    return `<div class="global-chain-card"><div class="global-chain-heading"><b>${escapeHtml(chainDisplayName(chain))}</b><span class="global-chain-summary">${sequence.length} aa；${chainLocations.length}个非重复差异区域；覆盖${chainCovered}个残基（${nice(chainCovered/sequence.length*100,1)}%）</span></div><div class="global-chain-track">${blocks.join("")}</div></div>`;
  });
  container.innerHTML=cards.join("");
  container.querySelectorAll("[data-overlap-position]").forEach(element=>element.onclick=()=>{
    state.sequenceOverlapSelection={chain:element.dataset.overlapChain,position:Number(element.dataset.overlapPosition)};
    renderGlobalSequenceOverview(locations);
  });
  renderSequenceOverlapDetails(locations,chains);
  const mappedFeatureIds=new Set(locations.flatMap(item=>item.featureIds));
  const sequenceEvidenceCount=(state.msms?.feature_evidence||[]).reduce((count,item)=>count+Number(Boolean(sequenceLocationForEvidence(item))),0);
  status.innerHTML=`<span class="structure-badge mapped">全序列差异定位</span><b>${escapeHtml(sampleShort(pair.test))}</b> 相对 <b>${escapeHtml(sampleShort(pair.reference))}</b>：${locations.length}个非重复差异肽段区域，关联${mappedFeatureIds.size}个Feature，覆盖LC/HC共${coveredResidues}个残基。MS2报告中共有${sequenceEvidenceCount}个含主鉴定或候选序列的Feature证据；common/low-confidence及无法解析链位置的证据不着色。`;
}
async function renderSequenceTrack(location){
  const status=$("sequenceMappingStatus"), track=$("sequenceTrack");
  if(!status||!track)return;
  if(!state.selectedFeatureGroupId){
    status.innerHTML='<span class="structure-badge missing">未选择Feature</span>请从热图或Feature列表中选择差异组分。';
    track.textContent="尚无可映射的MS2序列。";
    return;
  }
  if(state.msmsLoading||!state.msmsLoaded){
    status.innerHTML='<span class="structure-badge">载入MS2证据</span>正在读取候选肽段及残基位置。';
    track.textContent="正在载入序列证据…";
    return;
  }
  if(!location){
    status.innerHTML='<span class="structure-badge missing">无序列候选</span>该Feature当前没有通过阈值的MS2肽段，因此不能映射到蛋白序列或三级结构。';
    track.textContent="可继续查看覆盖隔离窗的MS2谱图，但不能据此标记具体残基。";
    return;
  }
  const fullSequence=cleanPeptideSequence(state.msms?.chains?.[location.chain]);
  let displaySequence=fullSequence||location.sequence, start=location.start, end=location.end;
  if(!fullSequence){ start=1; end=location.sequence.length; }
  const modifiedPositions=new Set();
  for(const match of location.modification.matchAll(/@(\d+)/g)){
    const relative=Number(match[1]);
    if(Number.isFinite(relative)&&Number.isFinite(start))modifiedPositions.add(start+relative-1);
  }
  const blocks=[];
  for(let offset=0;offset<displaySequence.length;offset+=50){
    const residues=displaySequence.slice(offset,offset+50).split("").map((aa,index)=>{
      const position=offset+index+1;
      const classes=["sequence-residue"];
      if(Number.isFinite(start)&&Number.isFinite(end)&&position>=start&&position<=end)classes.push("peptide",location.direction==="lower"?"diff-lower":"diff-higher");
      if(modifiedPositions.has(position))classes.push("modified");
      return `<span class="${classes.join(" ")}" title="${escapeHtml(`${location.chain||"chain"} ${position} ${aa}`)}">${aa}</span>`;
    }).join("");
    blocks.push(`<div class="sequence-block"><span class="sequence-position">${offset+1}</span><span class="sequence-residues">${residues}</span></div>`);
  }
  track.innerHTML=blocks.join("");
  const positionText=Number.isFinite(start)&&Number.isFinite(end)?`${start}-${end}`:"位置未解析";
  const confidence=location.confidence?` | 证据等级 ${escapeHtml(location.confidence)}`:"";
  const modification=location.modification?` | ${escapeHtml(location.modification)}`:"";
  const pair=selectedPair(), directionText=location.direction==="higher"?`${sampleShort(pair.test)} 高于 ${sampleShort(pair.reference)}`:(location.direction==="lower"?`${sampleShort(pair.reference)} 高于 ${sampleShort(pair.test)}`:"两样本接近");
  status.innerHTML=`<span class="structure-badge mapped">已映射到序列</span><b>${escapeHtml(location.chain||"未命名链")}</b> ${escapeHtml(positionText)} | ${escapeHtml(location.sequence)} | ${escapeHtml(directionText)}${modification}${confidence}`;
  const first=track.querySelector(".sequence-residue.peptide");
  if(first)track.scrollTop=Math.max(0,first.offsetTop-track.clientHeight/2);
}
function ensureStructureViewer(){
  if(state.structureViewer)return state.structureViewer;
  if(!window.$3Dmol)throw new Error("3Dmol.js未加载，无法创建三级结构视图");
  state.structureViewer=$3Dmol.createViewer($("structureViewer"),{backgroundColor:"#ffffff",antialias:true});
  return state.structureViewer;
}
function structureChainsFromModel(model){
  const byChain=new Map();
  const atoms=model?.selectedAtoms?.({atom:"CA"})||[];
  atoms.forEach(atom=>{
    const chain=String(atom.chain||"_"), resi=atom.resi, resn=String(atom.resn||"").toUpperCase();
    const aa=residueCode[resn]; if(!aa)return;
    const key=`${chain}|${resi}|${atom.icode||""}`;
    if(!byChain.has(chain))byChain.set(chain,{id:chain,residues:[],seen:new Set()});
    const item=byChain.get(chain); if(item.seen.has(key))return;
    item.seen.add(key); item.residues.push({resi,icode:String(atom.icode||""),resn,aa});
  });
  return [...byChain.values()].map(item=>({id:item.id,residues:item.residues,sequence:item.residues.map(res=>res.aa).join("")}));
}
function updateStructureChainOptions(){
  const select=$("structureChainSelect"); if(!select)return;
  const current=state.structureSelectedChain;
  select.innerHTML='<option value="">自动匹配</option>'+state.structureChains.map(item=>`<option value="${escapeHtml(item.id)}">${escapeHtml(item.id)} (${item.residues.length} aa)</option>`).join("");
  state.structureSelectedChain=state.structureChains.some(item=>item.id===current)?current:"";
  select.value=state.structureSelectedChain;
}
function localSequenceAlignment(query,target){
  const q=cleanPeptideSequence(query), t=cleanPeptideSequence(target);
  if(!q||!t)return null;
  const exactIndex=t.indexOf(q);
  if(exactIndex>=0)return {score:q.length*3,identity:1,coverage:1,exact:true,targetIndices:Array.from({length:q.length},(_,i)=>exactIndex+i),queryToTarget:Object.fromEntries(Array.from({length:q.length},(_,i)=>[i,exactIndex+i]))};
  const rows=q.length+1, cols=t.length+1;
  const score=Array.from({length:rows},()=>new Int16Array(cols));
  const trace=Array.from({length:rows},()=>new Uint8Array(cols));
  let bestScore=0,bestI=0,bestJ=0;
  for(let i=1;i<rows;i++)for(let j=1;j<cols;j++){
    const diagonal=score[i-1][j-1]+(q[i-1]===t[j-1]?3:-2);
    const up=score[i-1][j]-3, left=score[i][j-1]-3;
    const value=Math.max(0,diagonal,up,left); score[i][j]=value;
    trace[i][j]=value===0?0:(value===diagonal?1:(value===up?2:3));
    if(value>bestScore){bestScore=value;bestI=i;bestJ=j;}
  }
  let i=bestI,j=bestJ,matches=0,pairs=0; const queryIndices=new Set(),targetIndices=[],alignedPairs=[];
  while(i>0&&j>0&&score[i][j]>0){
    const direction=trace[i][j];
    if(direction===1){ i--;j--;pairs++;queryIndices.add(i);targetIndices.push(j);alignedPairs.push([i,j]);if(q[i]===t[j])matches++; }
    else if(direction===2){ i--; }
    else if(direction===3){ j--; }
    else break;
  }
  targetIndices.reverse(); alignedPairs.reverse();
  const identity=pairs?matches/pairs:0, coverage=queryIndices.size/q.length;
  return {score:bestScore,identity,coverage,exact:false,targetIndices:[...new Set(targetIndices)],queryToTarget:Object.fromEntries(alignedPairs)};
}
function preferredStructureChains(sourceChain){
  if(state.structureSelectedChain)return state.structureChains.filter(item=>item.id===state.structureSelectedChain);
  const sourceName=String(sourceChain||""), sourceSequence=cleanPeptideSequence(state.msms?.chains?.[sourceName]);
  if(!sourceSequence)return state.structureChains;
  const cached=state.structureSourceChainMap[sourceName];
  if(Array.isArray(cached)&&cached.length)return state.structureChains.filter(item=>cached.includes(item.id));
  const ranked=state.structureChains.map(chain=>{
    const alignment=localSequenceAlignment(sourceSequence,chain.sequence);
    return {chain,alignment,score:Number(alignment?.score)||0};
  }).filter(item=>item.score>0).sort((a,b)=>b.score-a.score||Number(b.alignment?.identity||0)-Number(a.alignment?.identity||0));
  const bestScore=ranked[0]?.score||0;
  const preferred=ranked.filter(item=>item.score>=bestScore*0.9&&Number(item.alignment?.identity||0)>=0.55).map(item=>item.chain.id);
  state.structureSourceChainMap[sourceName]=preferred.length?preferred:ranked.slice(0,1).map(item=>item.chain.id);
  return state.structureChains.filter(item=>state.structureSourceChainMap[sourceName].includes(item.id));
}
function peptideStructureMatch(location){
  if(!location||!state.structureChains.length)return null;
  const allowed=preferredStructureChains(location.chain);
  const matches=allowed.map(chain=>{
    const alignment=localSequenceAlignment(location.sequence,chain.sequence);
    if(!alignment)return null;
    return {...alignment,chain:chain.id,residues:alignment.targetIndices.map(index=>chain.residues[index]).filter(Boolean)};
  }).filter(Boolean).sort((a,b)=>Number(b.exact)-Number(a.exact)||b.score-a.score||b.coverage-a.coverage||b.identity-a.identity);
  const best=matches[0];
  if(!best||best.coverage<0.55||best.identity<0.55||!best.residues.length)return null;
  return best;
}
function applyStructureMapping(locations,selectedLocation){
  const viewer=state.structureViewer,status=$("structureMappingStatus");
  if(!viewer||!state.structureModel||!status)return;
  viewer.setStyle({}, {cartoon:{color:"#cbd5e1",opacity:0.22}});
  const residueDirections=new Map(), matchedRegions=[];
  for(const location of locations||[]){
    const match=peptideStructureMatch(location);
    if(!match)continue;
    matchedRegions.push({location,match});
    match.residues.forEach(residue=>{
      const key=`${match.chain}|${residue.resi}`, entry=residueDirections.get(key)||{chain:match.chain,resi:residue.resi,directions:new Set()};
      entry.directions.add(location.direction); residueDirections.set(key,entry);
    });
  }
  const targetStructureChains=[...new Set(Object.values(state.structureSourceChainMap).flat())];
  targetStructureChains.forEach(chain=>viewer.setStyle({chain},{cartoon:{color:"#cbd5e1",opacity:0.9}}));
  const styleGroups=new Map([["higher",new Map()],["lower",new Map()],["conflict",new Map()]]);
  residueDirections.forEach(entry=>{
    const key=entry.directions.size>1?"conflict":(entry.directions.has("higher")?"higher":"lower");
    const byChain=styleGroups.get(key); if(!byChain.has(entry.chain))byChain.set(entry.chain,[]);
    byChain.get(entry.chain).push(entry.resi);
  });
  const structureColors={higher:"#dc2626",lower:"#2563eb",conflict:"#7c3aed"};
  styleGroups.forEach((byChain,key)=>byChain.forEach((resi,chain)=>{
    viewer.setStyle({chain,resi:[...new Set(resi)]},{cartoon:{color:structureColors[key]}});
  }));
  let selectedMatch=null, selectedSelection=null, localizedResidues=[];
  if(selectedLocation){
    selectedMatch=peptideStructureMatch(selectedLocation);
    if(selectedMatch){
      const resi=[...new Set(selectedMatch.residues.map(item=>item.resi))];
      selectedSelection={chain:selectedMatch.chain,resi};
      const selectedColor=selectedLocation.direction==="lower"?structureColors.lower:(selectedLocation.direction==="higher"?structureColors.higher:structureColors.conflict);
      viewer.setStyle(selectedSelection,{cartoon:{color:selectedColor},stick:{colorscheme:"orangeCarbon",radius:0.18}});
      for(const siteMatch of String(selectedLocation.modification||"").matchAll(/@(\d+)/g)){
        const queryIndex=Number(siteMatch[1])-1, targetIndex=selectedMatch.queryToTarget?.[queryIndex];
        if(Number.isFinite(targetIndex)){
          const chain=state.structureChains.find(item=>item.id===selectedMatch.chain), residue=chain?.residues?.[targetIndex];
          if(residue)localizedResidues.push(residue);
        }
      }
      if(localizedResidues.length){
        viewer.setStyle({chain:selectedMatch.chain,resi:[...new Set(localizedResidues.map(item=>item.resi))]},{cartoon:{color:"#f59e0b"},stick:{color:"#f59e0b",radius:0.28},sphere:{color:"#f59e0b",scale:0.35}});
      }
    }
  }
  if(selectedSelection)viewer.zoomTo(selectedSelection);
  else if(matchedRegions.length)viewer.zoomTo({chain:[...new Set(matchedRegions.map(item=>item.match.chain))]});
  else viewer.zoomTo();
  viewer.render();
  const meta=state.structureSourceMeta||{}, templateBadge=meta.kind==="homologous"?'<span class="structure-badge candidate">同源模板，非药物精确结构</span>':'<span class="structure-badge mapped">结构已载入</span>';
  const templateNote=meta.note?` ${escapeHtml(meta.note)}`:"";
  let selectedText="";
  if(selectedLocation&&selectedMatch){
    const first=selectedMatch.residues[0],last=selectedMatch.residues[selectedMatch.residues.length-1];
    selectedText=`；选中Feature映射到结构链 <b>${escapeHtml(selectedMatch.chain)}</b> ${escapeHtml(first.resi)}-${escapeHtml(last.resi)}，覆盖率${nice(selectedMatch.coverage*100,1)}%、一致性${nice(selectedMatch.identity*100,1)}%`;
  }else if(selectedLocation){
    selectedText="；选中Feature在该结构中未达到55%覆盖率和55%一致性";
  }
  const pair=selectedPair(), directionLegend=`红色：${sampleShort(pair.test)} 高于 ${sampleShort(pair.reference)}；蓝色：${sampleShort(pair.reference)} 高于 ${sampleShort(pair.test)}`;
  status.innerHTML=`${templateBadge}${escapeHtml(state.structureLabel||"结构文件")}：全局${matchedRegions.length}/${(locations||[]).length}个差异肽段区域可映射。${escapeHtml(directionLegend)}；紫色表示重叠方向冲突，橙色棒状表示当前选中Feature${selectedText}。${templateNote}`;
}
function drawStructureModule(){
  const locations=globalSequenceLocations(), selectedLocation=selectedSequenceLocation();
  renderGlobalSequenceOverview(locations);
  renderSequenceTrack(selectedLocation);
  if(!state.structureModel)return;
  const pair=selectedPair();
  const signature=JSON.stringify([state.comparison,pair.reference,pair.test,state.selectedFeatureGroupId,selectedLocation?.sequence||"",state.structureLabel,state.structureSelectedChain,locations.map(item=>[item.chain,item.start,item.end,item.direction,item.featureIds])]);
  if(signature===state.structureRenderSignature)return;
  state.structureRenderSignature=signature;
  applyStructureMapping(locations,selectedLocation);
}
function clearLoadedStructure(){
  state.structureRequestId+=1;
  if(state.structureViewer)state.structureViewer.clear();
  state.structureModel=null; state.structureText=null; state.structureFormat=null; state.structureLabel=null; state.structureSourceMeta=null; state.structureChains=[]; state.structureSourceChainMap={}; state.structureSelectedChain=""; state.structureRenderSignature=""; state.structureLoading=false;
  updateStructureChainOptions();
  const status=$("structureMappingStatus"); if(status)status.textContent="请上传结构文件或输入PDB ID。";
  const link=$("structureSourceLink"); if(link){link.hidden=true;link.removeAttribute("href");}
}
function loadStructureText(text,format,label,sourceUrl="",sourceMeta=null){
  const viewer=ensureStructureViewer();
  viewer.clear();
  const model=viewer.addModel(text,format);
  if(!model)throw new Error("结构文件无法解析");
  state.structureModel=model; state.structureText=text; state.structureFormat=format; state.structureLabel=label; state.structureSourceMeta=sourceMeta; state.structureChains=structureChainsFromModel(model); state.structureSourceChainMap={}; state.structureSelectedChain=""; state.structureRenderSignature="";
  if(!state.structureChains.length){viewer.clear();state.structureModel=null;throw new Error("结构中没有可用于蛋白序列映射的CA原子");}
  updateStructureChainOptions();
  const link=$("structureSourceLink");
  if(link){link.hidden=!sourceUrl;if(sourceUrl)link.href=sourceUrl;}
  drawStructureModule();
}
async function persistPdbStructure(pdbId,label){
  if(String(DATA?.analysis_metadata?.ms2?.status||"")!=="completed")throw new Error("该任务没有可用的 MS2/FASTA 结果，不能保存结构映射。");
  const response=await fetch(apiUrl("/api/task-reference/pdb"),{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({pdb_id:pdbId,name:label||`PDB ${pdbId}`}),cache:"no-store"});
  if(!response.ok)throw new Error((await response.text())||`PDB结构保存失败 HTTP ${response.status}`);
  return response.json();
}
async function loadPdbStructure(preset=null){
  const input=$("pdbIdInput"),button=$("loadPdbStructure"),status=$("structureMappingStatus");
  const sourceMeta=preset&&preset.pdbId?preset:null;
  const pdbId=String(sourceMeta?.pdbId||input?.value||"").trim().toUpperCase();
  if(!/^[A-Z0-9]{4,12}$/.test(pdbId)){status.textContent="请输入4-12位PDB ID。";return;}
  if(input)input.value=pdbId;
  const requestId=++state.structureRequestId; state.structureLoading=true; button.disabled=true; status.textContent=`正在从RCSB PDB获取 ${pdbId}…`;
  try{
    const response=await fetch(apiUrl(`/api/pdb-structure?pdb_id=${encodeURIComponent(pdbId)}`),{cache:"no-store"});
    const text=await response.text();
    if(!response.ok)throw new Error(`PDB获取失败（HTTP ${response.status}）${text?`: ${text.slice(0,160)}`:""}`);
    if(requestId!==state.structureRequestId)return;
    const label=sourceMeta?.label?`PDB ${pdbId} · ${sourceMeta.label}`:`PDB ${pdbId}`;
    loadStructureText(text,"cif",label,`https://www.rcsb.org/structure/${encodeURIComponent(pdbId)}`,sourceMeta);
    if(!preset){
      status.textContent=`PDB ${pdbId} 已载入，正在保存到任务…`;
      const saved=await persistPdbStructure(pdbId,label);
      DATA.analysis_metadata=DATA.analysis_metadata||{};
      DATA.analysis_metadata.references=saved.references||DATA.analysis_metadata.references;
      DATA.analysis_metadata.structure_mapping={status:"available",source:"pdb_id",pdb_id:pdbId,path:saved.structure?.path||""};
      state.structureSourceMeta=saved.structure?.source_meta||state.structureSourceMeta;
      status.textContent=`PDB ${pdbId} 已载入并保存，重新打开报告时会自动恢复。`;
    }
  }catch(error){
    if(requestId===state.structureRequestId)status.textContent=String(error?.message||error);
  }finally{
    if(requestId===state.structureRequestId){state.structureLoading=false;button.disabled=false;}
  }
}
async function maybeLoadProjectStructure(){
  if(state.structureModel||state.structureLoading||state.structureAutoAttemptedComparison===state.comparison)return;
  if(!DATA?.analysis_metadata){
    const preset=PROJECT_STRUCTURE_DEFAULTS[String(DATA?.project_id||state.msms?.project_id||"")];
    if(preset){ state.structureAutoAttemptedComparison=state.comparison; await loadPdbStructure(preset); }
    return;
  }
  state.structureAutoAttemptedComparison=state.comparison;
  const mapping=DATA?.analysis_metadata?.structure_mapping||{};
  if(mapping.status!=="available"){
    const status=$("structureMappingStatus");
    if(status)status.textContent=String(mapping.reason||"未提供同时满足映射条件的 FASTA 与结构资料。");
    return;
  }
  try{
    const payload=await fetchJson(apiUrl("/api/task-reference"));
    if(payload.available&&payload.text)loadStructureText(payload.text,payload.format||"pdb",payload.label||"任务结构",payload.source_url||"",payload.source_meta||null);
    else { const status=$("structureMappingStatus"); if(status)status.textContent=String(payload.reason||mapping.reason||"结构不可用。"); }
  }catch(error){ const status=$("structureMappingStatus"); if(status)status.textContent=String(error?.message||error); }
}
async function loadUploadedStructure(file,save=true){
  const status=$("structureMappingStatus");
  if(!file)return;
  if(file.size>40*1024*1024){status.textContent="结构文件超过40 MB，未载入。";return;}
  status.textContent=`正在读取 ${file.name}…`;
  try{
    const text=await file.text(), extension=String(file.name.split(".").pop()||"").toLowerCase();
    const format=["cif","mmcif"].includes(extension)?"cif":"pdb";
    loadStructureText(text,format,file.name);
    if(save){
      if(String(DATA?.analysis_metadata?.ms2?.status||"")!=="completed")throw new Error("该任务没有可用的 MS2/FASTA 结果，不能保存结构映射。");
      const body=new FormData(); body.append("structure_file",file,file.name);
      status.textContent=`正在自动保存 ${file.name} 并更新任务结构映射…`;
      const response=await fetch(apiUrl("/api/task-reference/upload"),{method:"POST",body,cache:"no-store"});
      if(!response.ok)throw new Error((await response.text())||`结构上传失败 HTTP ${response.status}`);
      const payload=await response.json();
      DATA.analysis_metadata=DATA.analysis_metadata||{};
      DATA.analysis_metadata.structure_mapping={status:"available",source:"report_upload",path:payload.structure?.path||""};
      state.structureAutoAttemptedComparison=state.comparison;
      state.structureRenderSignature="";
      status.textContent="结构已保存，重新打开报告时会自动恢复。";
      drawStructureModule();
    }
  }catch(error){status.textContent=String(error?.message||error);}
}
function wireStructureControls(){
  const mappingStatus=String(DATA?.analysis_metadata?.structure_mapping?.status||"");
  const hasFasta=String(DATA?.analysis_metadata?.ms2?.status||"")==="completed";
  const mappingAllowed=!DATA?.analysis_metadata||hasFasta;
  [$("pdbIdInput"),$("loadPdbStructure"),$("structureFileInput"),$("saveStructureFile")].forEach(element=>{
    if(element){ element.disabled=!mappingAllowed; element.title=mappingAllowed?"":"需要在任务入口提供 FASTA 序列后才能进行映射"; }
  });
  $("loadPdbStructure").onclick=()=>loadPdbStructure();
  $("pdbIdInput").onkeydown=event=>{if(event.key==="Enter")loadPdbStructure();};
  $("structureFileInput").onchange=event=>loadUploadedStructure(event.target.files?.[0],true);
  $("saveStructureFile").onclick=()=>loadUploadedStructure($("structureFileInput").files?.[0],true);
  $("structureChainSelect").onchange=event=>{state.structureSelectedChain=event.target.value;state.structureRenderSignature="";drawStructureModule();};
  $("resetStructureView").onclick=()=>{if(state.structureViewer){state.structureViewer.zoomTo();state.structureViewer.render();}};
}
function modificationQuantitationStatusLabel(status){
  return {formal_relative_quantitation:"正式相对定量",tentative_relative_quantitation:"暂定相对定量",insufficient_complementary_forms:"互补形式不足"}[status]||status||"";
}
function renderProteolyticFamilies(pair,query){
  const container=$("proteolyticFamilySummary"),families=Array.isArray(state.msms?.proteolytic_modification_families)?state.msms.proteolytic_modification_families:[];
  if(!container)return new Set();
  const visible=families.filter(family=>!query||[family.chain,...(family.members||[]).map(member=>member.sequence),...(family.site_conclusions||[]).map(site=>site.event_label)].join(" ").toLowerCase().includes(query));
  const memberIds=new Set(visible.flatMap(family=>(family.members||[]).map(member=>String(member.quantitation_id||""))));
  if(!visible.length){container.innerHTML="";return memberIds;}
  const orientedSiteComparison=site=>{
    const comparisons=site.pairwise_comparisons||[];
    const direct=comparisons.find(item=>item.reference_sample===pair.reference&&item.test_sample===pair.test);
    if(direct)return {...direct,evidence:(direct.evidence||[]).map(item=>({...item}))};
    const reverse=comparisons.find(item=>item.reference_sample===pair.test&&item.test_sample===pair.reference);
    if(!reverse)return null;
    return {...reverse,reference_sample:pair.reference,test_sample:pair.test,direction:reverse.direction==="increased"?"decreased":(reverse.direction==="decreased"?"increased":reverse.direction),min_percentage_point_difference:reverse.max_percentage_point_difference==null?null:-Number(reverse.max_percentage_point_difference),max_percentage_point_difference:reverse.min_percentage_point_difference==null?null:-Number(reverse.min_percentage_point_difference),median_percentage_point_difference:reverse.median_percentage_point_difference==null?null:-Number(reverse.median_percentage_point_difference),evidence:(reverse.evidence||[]).map(item=>({...item,reference_level:item.test_level,test_level:item.reference_level,percentage_point_difference:-Number(item.percentage_point_difference||0),direction:item.direction==="increased"?"decreased":(item.direction==="decreased"?"increased":item.direction)}))};
  };
  container.innerHTML=visible.map(family=>{
    const sites=(family.site_conclusions||[]).map(site=>({site,comparison:orientedSiteComparison(site)})).filter(item=>item.comparison).sort((a,b)=>Math.abs(Number(b.comparison.median_percentage_point_difference||0))-Math.abs(Number(a.comparison.median_percentage_point_difference||0)));
    const lead=sites[0],comparison=lead?.comparison,site=lead?.site;
    let conclusion='<span class="note">嵌套肽段已归组，但当前样本对尚无共同定位修饰结论。</span>';
    if(comparison){
      if(comparison.consistent){
        const delta=Number(comparison.median_percentage_point_difference||0),higher=delta>=0?pair.test:pair.reference;
        const magnitudes=[Math.abs(Number(comparison.min_percentage_point_difference||0)),Math.abs(Number(comparison.max_percentage_point_difference||0))].sort((a,b)=>a-b);
        conclusion=`<b>${escapeHtml(site.event_label)}</b>：${escapeHtml(sampleShort(higher))} 较高；${comparison.support_count}/${comparison.support_count} 种酶切形态方向一致，变化范围 <b>${(magnitudes[0]*100).toFixed(2)}–${(magnitudes[1]*100).toFixed(2)} 个百分点</b>。`;
      }else conclusion=`<b>${escapeHtml(site.event_label)}</b>：不同酶切形态的修饰比例方向不一致，暂不形成位点级高价值结论。`;
    }
    const totals=(family.members||[]).map(member=>{const area=member.total_area_by_sample||{},referenceArea=Number(area[pair.reference]||0),testArea=Number(area[pair.test]||0),fold=referenceArea>0&&testArea>0?testArea/referenceArea:null;return {member,referenceArea,testArea,fold};});
    const folds=totals.map(item=>item.fold).filter(Boolean),redistribution=folds.some(fold=>fold>=1.2)&&folds.some(fold=>fold<=1/1.2);
    const combinedReference=totals.reduce((sum,item)=>sum+item.referenceArea,0),combinedTest=totals.reduce((sum,item)=>sum+item.testArea,0),combinedFold=combinedReference>0&&combinedTest>0?combinedTest/combinedReference:null;
    const digestionText=redistribution?`<span class="digestion-warning">检测到酶切分配差异：</span>长短肽总响应呈相反方向，不能把任一肽段的总响应差异直接解释为蛋白或修饰丰度差异。`:'各酶切形态总响应未出现明显相反方向。';
    let combinedText="";
    if(combinedFold!=null){const higher=combinedFold>=1?pair.test:pair.reference,ratio=combinedFold>=1?combinedFold:1/combinedFold;combinedText=`不同序列未校正合计仅作辅助：${escapeHtml(sampleShort(higher))} ${ratio.toFixed(2)}x。`;}
    const memberRows=totals.map(({member,referenceArea,testArea,fold})=>{
      const memberEvidence=comparison?.evidence?.find(item=>item.quantitation_id===member.quantitation_id),delta=memberEvidence?Number(memberEvidence.percentage_point_difference||0):null;
      let totalDirection="NA";if(fold!=null){const higher=fold>=1?pair.test:pair.reference,ratio=fold>=1?fold:1/fold;totalDirection=`${escapeHtml(sampleShort(higher))} ${ratio.toFixed(2)}x`;}
      return `<tr><td><b>${escapeHtml(member.sequence||"")}</b><br><span class="note">${member.start||""}-${member.end||""}</span></td><td>${memberEvidence?`${(Number(memberEvidence.reference_level||0)*100).toFixed(2)}% → ${(Number(memberEvidence.test_level||0)*100).toFixed(2)}%`:"NA"}</td><td class="${delta!=null&&delta>0.005?"mod-delta-up":(delta!=null&&delta<-0.005?"mod-delta-down":"mod-delta-stable")}">${delta==null?"NA":`${delta>0?"+":""}${(delta*100).toFixed(2)} pp`}</td><td>${totalDirection}</td><td>${memberEvidence?.evidence_grade||"—"}级</td></tr>`;
    }).join("");
    const high=Boolean(site?.high_value_consensus&&comparison?.consistent&&comparison?.all_b_grade);
    return `<div class="proteolytic-family-card ${high?"high-value":""}"><div class="proteolytic-family-title"><span>位点结论优先</span>${high?'<span class="mod-quant-badge high-value">跨酶切形态高价值</span>':''}${redistribution?'<span class="mod-quant-badge high-value">疑似酶切分配差异</span>':''}<span class="note">${escapeHtml(family.chain||"")} ${family.start||""}-${family.end||""}</span></div><div class="proteolytic-family-conclusion">${conclusion}</div><div class="proteolytic-family-grid"><div class="proteolytic-family-note">${digestionText}<br><span class="note">${combinedText} 不同长度肽段没有响应因子时不计算绝对合并修饰率。</span></div><div class="scroll"><table class="proteolytic-member-table"><tr><th>酶切形态</th><th>${escapeHtml(site?.event_label||"共同修饰")} 构成</th><th>变化</th><th>该肽段总响应</th><th>证据</th></tr>${memberRows}</table></div></div></div>`;
  }).join("");
  return memberIds;
}
function renderModificationQuantitationLegacy(){
  const table=$("modQuantTable"),summary=$("modQuantSummary"),filter=$("modQuantFilter"),statusSelect=$("modQuantStatus");
  if(!table||!summary)return;
  if(filter&&!filter.dataset.bound){filter.dataset.bound="1";filter.oninput=renderModificationQuantitation;}
  if(statusSelect&&!statusSelect.dataset.bound){statusSelect.dataset.bound="1";statusSelect.onchange=renderModificationQuantitation;}
  if(state.msmsLoading||!state.msmsLoaded){
    summary.innerHTML='<span class="mod-quant-badge">正在载入修饰形式定量…</span>';
    if($("proteolyticFamilySummary"))$("proteolyticFamilySummary").innerHTML="";
    table.innerHTML="";
    return;
  }
  if(state.msmsError){
    summary.innerHTML=`<span class="mod-quant-badge">不可用</span><span class="note">${escapeHtml(state.msmsError)}</span>`;
    if($("proteolyticFamilySummary"))$("proteolyticFamilySummary").innerHTML="";
    table.innerHTML="";
    return;
  }
  const all=Array.isArray(state.msms?.modification_level_quantitation)?state.msms.modification_level_quantitation:[];
  const pair=selectedPair(), query=String(filter?.value||"").trim().toLowerCase(), mode=statusSelect?.value||"quantifiable";
  const nestedFamilyIds=renderProteolyticFamilies(pair,query);
  const evaluated=all.map(item=>{
    const changes=(item.forms||[]).filter(form=>form.included_in_denominator).map(form=>{
      const levels=form.relative_level_by_sample||{},areas=form.normalized_area_by_sample||{},referenceLevel=levels[pair.reference],testLevel=levels[pair.test];
      if(referenceLevel==null||testLevel==null)return null;
      const referenceArea=Number(areas[pair.reference]||0),testArea=Number(areas[pair.test]||0),delta=Number(testLevel)-Number(referenceLevel);
      const responseFold=referenceArea>0&&testArea>0?testArea/referenceArea:null;
      return {form,referenceLevel:Number(referenceLevel),testLevel:Number(testLevel),referenceArea,testArea,responseFold,delta};
    }).filter(Boolean).sort((a,b)=>Math.abs(b.delta)-Math.abs(a.delta));
    const totalAreas=item.total_area_by_sample||{},familyReferenceArea=Number(totalAreas[pair.reference]||0),familyTestArea=Number(totalAreas[pair.test]||0);
    const familyResponseFold=familyReferenceArea>0&&familyTestArea>0?familyTestArea/familyReferenceArea:null;
    const familyResponseImbalance=familyResponseFold!=null&&(familyResponseFold>2||familyResponseFold<0.5);
    return {...item,currentChanges:changes,currentStrongest:changes[0]||null,familyReferenceArea,familyTestArea,familyResponseFold,familyResponseImbalance};
  }).filter(item=>{
    if(mode==="formal"&&item.quantitation_status!=="formal_relative_quantitation")return false;
    if(mode==="quantifiable"&&!item.quantifiable)return false;
    if(!query)return true;
    return [item.event_name,item.event_type,item.sequence,item.chain,...(item.linked_feature_ids||[]),...(item.forms||[]).map(form=>form.label)].join(" ").toLowerCase().includes(query);
  }).sort((a,b)=>Math.abs(Number(b.currentStrongest?.delta||0))-Math.abs(Number(a.currentStrongest?.delta||0))||Number(b.high_value_candidate)-Number(a.high_value_candidate)||Number(a.rank)-Number(b.rank));
  const formal=all.filter(item=>item.quantitation_status==="formal_relative_quantitation").length,standaloneHigh=all.filter(item=>item.high_value_candidate&&!nestedFamilyIds.has(String(item.quantitation_id||""))).length,siteHigh=(state.msms?.proteolytic_modification_families||[]).filter(item=>item.high_value_consensus).length;
  summary.innerHTML=`<span class="mod-quant-badge">${all.length} 个肽段形式族</span><span class="mod-quant-badge formal">${formal} 个正式相对定量</span><span class="mod-quant-badge high-value">${siteHigh} 个跨酶切形态高价值结论</span>${standaloneHigh?`<span class="mod-quant-badge high-value">${standaloneHigh} 个独立高价值候选</span>`:""}<span class="note">当前显示 ${evaluated.length} 个；${escapeHtml(sampleShort(pair.test))} 对 ${escapeHtml(sampleShort(pair.reference))}</span>`;
  if(!evaluated.length){table.innerHTML='<tr><td class="note">当前筛选条件下没有可显示的修饰形式族。</td></tr>';return;}
  const formPalette=["#2563eb","#f59e0b","#10b981","#8b5cf6","#ef4444","#06b6d4","#84cc16","#f97316","#64748b"];
  table.innerHTML='<tr><th></th><th>rank</th><th>事件 / 位置</th><th>肽段</th><th>主要变化</th><th>定量等级</th><th>形式数</th></tr>'+evaluated.map(item=>{
    const changes=item.currentChanges||[], strongest=item.currentStrongest;
    let conclusion='<span class="note">暂不能形成组间比例结论</span>';
    if(strongest){
      const higher=strongest.delta>=0?pair.test:pair.reference;
      const higherLevel=strongest.delta>=0?strongest.testLevel:strongest.referenceLevel, lowerLevel=strongest.delta>=0?strongest.referenceLevel:strongest.testLevel;
      let responseText="";
      if(strongest.responseFold!=null){
        const responseHigher=strongest.responseFold>=1?pair.test:pair.reference,responseRatio=strongest.responseFold>=1?strongest.responseFold:1/strongest.responseFold;
        responseText=`；${escapeHtml(sampleShort(responseHigher))} 的该形式归一化响应为 ${responseRatio.toFixed(2)}x`;
      }
      conclusion=`<b>${escapeHtml(sampleShort(higher))}</b> 的“${escapeHtml(strongest.form.label)}”构成较高 <b>${(Math.abs(strongest.delta)*100).toFixed(2)} 个百分点</b><br><span class="note">${(lowerLevel*100).toFixed(2)}% → ${(higherLevel*100).toFixed(2)}%${responseText}</span>`;
    }
    const location=`${escapeHtml(item.chain||"")} ${item.start||""}-${item.end||""}`;
    const familyKey=String(item.quantitation_id||item.family_id||item.rank||""), expanded=state.expandedModificationFamilies.has(familyKey);
    const confidenceGrade=changes.some(change=>String(change.form.ms2_confidence||"").startsWith("C_"))?"含C级定位":"B级形式";
    const isNestedMember=nestedFamilyIds.has(String(item.quantitation_id||"")),standaloneHighValue=item.high_value_candidate&&!isNestedMember;
    const imbalanceBadge=item.familyResponseImbalance?' <span class="mod-quant-badge high-value">该肽段总响应不平衡</span>':"";
    const contextBadge=isNestedMember?' <span class="mod-quant-badge">酶切族成员</span>':"";
    const summaryRow=`<tr class="${standaloneHighValue?"high-value":""}"><td><button type="button" class="mod-quant-toggle" data-mod-quant-toggle="${escapeHtml(familyKey)}" aria-expanded="${expanded}">${expanded?"−":"+"}</button></td><td>${item.rank||""}</td><td class="wrap"><div class="mod-quant-event">${escapeHtml(item.event_name||"")}${standaloneHighValue?' <span class="mod-quant-badge high-value">独立高价值候选</span>':""}${contextBadge}${imbalanceBadge}</div><div class="mod-quant-site">${location}</div></td><td class="wrap"><b>${escapeHtml(item.sequence||"")}</b></td><td class="mod-quant-conclusion ${standaloneHighValue?"high-value":""}">${conclusion}${isNestedMember?'<br><span class="note">该肽段属于上方酶切形态族，最终判断以位点级一致性结论为准。</span>':""}${item.familyResponseImbalance?'<br><span class="note">该肽段自身总响应相差超过2倍，不单独解释为修饰或蛋白丰度变化。</span>':""}</td><td class="wrap"><b>${escapeHtml(modificationQuantitationStatusLabel(item.quantitation_status))}</b><br><span class="note">${escapeHtml(confidenceGrade)}</span></td><td>${changes.length}</td></tr>`;
    if(!expanded)return summaryRow;
    const barForSample=sample=>`<div class="mod-composition-row"><div class="mod-composition-sample" title="${escapeHtml(sample)}">${escapeHtml(sampleShort(sample))}</div><div class="mod-composition-bar">${changes.map((change,index)=>{const level=sample===pair.reference?change.referenceLevel:change.testLevel;return `<span class="mod-composition-segment" style="width:${Math.max(0,level*100)}%;background:${formPalette[index%formPalette.length]}" title="${escapeHtml(change.form.label)}：${(level*100).toFixed(2)}%"></span>`;}).join("")}</div></div>`;
    const matrixRows=changes.map((change,index)=>{
      const sign=change.delta>0?"+":"", deltaClass=change.delta>0.005?"mod-delta-up":(change.delta<-0.005?"mod-delta-down":"mod-delta-stable");
      const evidence=String(change.form.ms2_confidence||"").startsWith("B_")?"B级":(String(change.form.ms2_confidence||"").startsWith("C_")?"C级":"D级");
      let responseFoldText="NA";
      if(change.responseFold!=null){const responseHigher=change.responseFold>=1?pair.test:pair.reference,responseRatio=change.responseFold>=1?change.responseFold:1/change.responseFold;responseFoldText=`${escapeHtml(sampleShort(responseHigher))} ${responseRatio.toFixed(2)}x`;}
      else if(change.referenceArea<=0&&change.testArea>0)responseFoldText=`仅 ${escapeHtml(sampleShort(pair.test))} 检出`;
      else if(change.testArea<=0&&change.referenceArea>0)responseFoldText=`仅 ${escapeHtml(sampleShort(pair.reference))} 检出`;
      return `<tr><td><div class="mod-form-label"><span class="mod-form-swatch" style="background:${formPalette[index%formPalette.length]}"></span><span>${escapeHtml(change.form.label)}</span></div></td><td>${(change.referenceLevel*100).toFixed(2)}%</td><td>${(change.testLevel*100).toFixed(2)}%</td><td class="${deltaClass}">${sign}${(change.delta*100).toFixed(2)}</td><td>${responseFoldText}</td><td>${evidence}</td></tr>`;
    }).join("")||'<tr><td colspan="6" class="note">互补形式不足，暂不能计算构成比例</td></tr>';
    let familyResponseText="";
    if(item.familyResponseFold!=null){const familyHigher=item.familyResponseFold>=1?pair.test:pair.reference,familyRatio=item.familyResponseFold>=1?item.familyResponseFold:1/item.familyResponseFold;familyResponseText=`同肽段形式族总归一化响应：${escapeHtml(sampleShort(familyHigher))} 为 ${familyRatio.toFixed(2)}x。`;}
    const featureButtons=(item.linked_feature_ids||[]).map(id=>`<button type="button" class="mod-quant-feature" data-mod-quant-feature="${escapeHtml(id)}">${escapeHtml(id)}</button>`).join("");
    const detailRow=`<tr class="mod-quant-detail-row"><td colspan="7"><div class="mod-quant-detail"><div class="mod-quant-detail-grid"><div class="mod-composition-panel"><div class="mod-detail-title">两个样本的形式构成（每条合计100%）</div>${barForSample(pair.reference)}${barForSample(pair.test)}<div class="note">色块对应右侧各形式；极低比例形式可能仅在数值表中可见。${familyResponseText}</div></div><div class="mod-form-panel"><div class="mod-detail-title">逐形式定量：构成比例与该形式自身响应分开显示</div><table class="mod-form-matrix"><tr><th>形式</th><th>${escapeHtml(sampleShort(pair.reference))} 构成</th><th>${escapeHtml(sampleShort(pair.test))} 构成</th><th>构成变化/百分点</th><th>形式响应倍数</th><th>证据</th></tr>${matrixRows}</table></div></div><details class="mod-feature-details"><summary>关联 Feature（${(item.linked_feature_ids||[]).length}）</summary><div>${featureButtons||'<span class="note">无</span>'}</div></details></div></td></tr>`;
    return summaryRow+detailRow;
  }).join("");
  table.querySelectorAll("[data-mod-quant-toggle]").forEach(button=>button.onclick=()=>{const key=String(button.dataset.modQuantToggle||"");if(state.expandedModificationFamilies.has(key))state.expandedModificationFamilies.delete(key);else state.expandedModificationFamilies.add(key);renderModificationQuantitation();});
  table.querySelectorAll("[data-mod-quant-feature]").forEach(button=>button.onclick=()=>selectMappedFeature(button.dataset.modQuantFeature));
}
function modOrientedSiteComparison(site,pair){
  const comparisons=site?.pairwise_comparisons||[];
  const direct=comparisons.find(item=>item.reference_sample===pair.reference&&item.test_sample===pair.test);
  if(direct)return {...direct,evidence:(direct.evidence||[]).map(item=>({...item}))};
  const reverse=comparisons.find(item=>item.reference_sample===pair.test&&item.test_sample===pair.reference);
  if(!reverse)return null;
  const flipDirection=direction=>direction==="increased"?"decreased":(direction==="decreased"?"increased":direction);
  return {...reverse,reference_sample:pair.reference,test_sample:pair.test,direction:flipDirection(reverse.direction),min_percentage_point_difference:reverse.max_percentage_point_difference==null?null:-Number(reverse.max_percentage_point_difference),max_percentage_point_difference:reverse.min_percentage_point_difference==null?null:-Number(reverse.min_percentage_point_difference),median_percentage_point_difference:reverse.median_percentage_point_difference==null?null:-Number(reverse.median_percentage_point_difference),evidence:(reverse.evidence||[]).map(item=>({...item,reference_level:item.test_level,test_level:item.reference_level,percentage_point_difference:-Number(item.percentage_point_difference||0),direction:flipDirection(item.direction)}))};
}
function modFormatArea(value){
  const number=Number(value||0);
  if(!Number.isFinite(number))return "NA";
  return number.toLocaleString("zh-CN",{maximumFractionDigits:0});
}
function modEvaluateItem(item,pair){
  const changes=(item.forms||[]).filter(form=>form.included_in_denominator).map(form=>{
    const levels=form.relative_level_by_sample||{},areas=form.normalized_area_by_sample||{},referenceLevel=levels[pair.reference],testLevel=levels[pair.test];
    if(referenceLevel==null||testLevel==null)return null;
    const referenceArea=Number(areas[pair.reference]||0),testArea=Number(areas[pair.test]||0),delta=Number(testLevel)-Number(referenceLevel);
    return {form,referenceLevel:Number(referenceLevel),testLevel:Number(testLevel),referenceArea,testArea,responseFold:referenceArea>0&&testArea>0?testArea/referenceArea:null,delta};
  }).filter(Boolean).sort((a,b)=>Math.abs(b.delta)-Math.abs(a.delta));
  const totals=item.total_area_by_sample||{},familyReferenceArea=Number(totals[pair.reference]||0),familyTestArea=Number(totals[pair.test]||0),familyResponseFold=familyReferenceArea>0&&familyTestArea>0?familyTestArea/familyReferenceArea:null;
  return {...item,currentChanges:changes,currentStrongest:changes[0]||null,familyReferenceArea,familyTestArea,familyResponseFold,familyResponseImbalance:familyResponseFold!=null&&(familyResponseFold>2||familyResponseFold<0.5)};
}
function modItemDetailRow(item,pair,formPalette){
  const changes=item.currentChanges||[],mode=state.modQuantViewMode||"composition",maxTotal=Math.max(item.familyReferenceArea,item.familyTestArea,1);
  const barForSample=sample=>{
    const isReference=sample===pair.reference,total=isReference?item.familyReferenceArea:item.familyTestArea;
    const scaleWidth=mode==="composition"?100:Math.max(0,total/maxTotal*100);
    const segments=changes.map((change,index)=>{
      const value=isReference?change.referenceArea:change.testArea,level=isReference?change.referenceLevel:change.testLevel;
      const segmentWidth=mode==="composition"?Math.max(0,level*100):(total>0?Math.max(0,value/total*100):0);
      const title=mode==="composition"?`${change.form.label}：${(level*100).toFixed(2)}%`:`${change.form.label}：${modFormatArea(value)}`;
      return `<span class="mod-composition-segment" style="width:${segmentWidth}%;background:${formPalette[index%formPalette.length]}" title="${escapeHtml(title)}"></span>`;
    }).join("");
    return `<div class="mod-composition-row"><div class="mod-composition-sample" title="${escapeHtml(sample)}">${escapeHtml(sampleShort(sample))}</div><div class="mod-composition-bar"><div class="mod-composition-scale" data-target-width="${scaleWidth.toFixed(4)}">${segments}</div></div></div>`;
  };
  const matrixRows=changes.map((change,index)=>{
    const evidence=String(change.form.ms2_confidence||"").startsWith("B_")?"B级":(String(change.form.ms2_confidence||"").startsWith("C_")?"C级":"D级");
    let responseFoldText="NA";
    if(change.responseFold!=null){const responseHigher=change.responseFold>=1?pair.test:pair.reference,responseRatio=change.responseFold>=1?change.responseFold:1/change.responseFold;responseFoldText=`${escapeHtml(sampleShort(responseHigher))} ${responseRatio.toFixed(2)}x`;}
    else if(change.referenceArea<=0&&change.testArea>0)responseFoldText=`仅 ${escapeHtml(sampleShort(pair.test))} 检出`;
    else if(change.testArea<=0&&change.referenceArea>0)responseFoldText=`仅 ${escapeHtml(sampleShort(pair.reference))} 检出`;
    const linkedIds=(change.form.linked_feature_ids||[]).map(String).filter(Boolean),formKey=`${item.quantitation_id||item.family_id||""}::${change.form.form_id||change.form.label||index}`,selected=state.selectedModificationFormKey===formKey,targetable=linkedIds.length||Number.isFinite(Number(change.form.neutral_mass))&&Boolean((change.form.consensus_rts||[]).length);
    const formText=targetable?`<button type="button" class="mod-form-link" data-mod-form-key="${escapeHtml(formKey)}" data-mod-form-links="${escapeHtml(linkedIds.join("|"))}" title="${linkedIds.length?"点击后联动 Feature、热图、XIC、MS/MS、序列及结构模块":"该形式不在差异Feature表中；点击后联动目标XIC、MS/MS、序列及结构模块"}">${escapeHtml(change.form.label)}</button><span class="mod-form-link-count">↗ ${linkedIds.length||"XIC"}</span>`:`<span title="该形式暂无可关联信号">${escapeHtml(change.form.label)}</span>`;
    const label=`<div class="mod-form-label"><span class="mod-form-swatch" style="background:${formPalette[index%formPalette.length]}"></span>${formText}</div>`;
    if(mode==="xic")return `<tr class="${selected?"mod-form-selected":""}"><td>${label}</td><td>${modFormatArea(change.referenceArea)}</td><td>${modFormatArea(change.testArea)}</td><td>${responseFoldText}</td><td>${evidence}</td></tr>`;
    const sign=change.delta>0?"+":"",deltaClass=change.delta>0.005?"mod-delta-up":(change.delta<-0.005?"mod-delta-down":"mod-delta-stable");
    return `<tr class="${selected?"mod-form-selected":""}"><td>${label}</td><td>${(change.referenceLevel*100).toFixed(2)}%</td><td>${(change.testLevel*100).toFixed(2)}%</td><td class="${deltaClass}">${sign}${(change.delta*100).toFixed(2)}</td><td>${responseFoldText}</td><td>${evidence}</td></tr>`;
  }).join("")||`<tr><td colspan="${mode==="xic"?5:6}" class="note">互补形式不足，暂不能计算</td></tr>`;
  let familyResponseText="";
  if(item.familyResponseFold!=null){const higher=item.familyResponseFold>=1?pair.test:pair.reference,ratio=item.familyResponseFold>=1?item.familyResponseFold:1/item.familyResponseFold;familyResponseText=`同肽段形式族总归一化 XIC 响应：${escapeHtml(sampleShort(higher))} 为 ${ratio.toFixed(2)}x。`;}
  const panelTitle=mode==="xic"?"两个样本的归一化 XIC 面积（条形总长度使用同一尺度）":"两个样本的形式构成（每条合计100%）";
  const panelNote=mode==="xic"?`色块宽度表示各形式的 XIC 面积；较长样本条代表该肽段总响应更高。数值已按总 TIC 归一化以便跨样本比较，不是未经校正的原始面积。${familyResponseText}`:`色块对应右侧各形式；极低比例形式可能仅在数值表中可见。${familyResponseText}`;
  const tableTitle=mode==="xic"?"逐形式归一化 XIC 面积":"逐形式定量：构成比例与形式自身响应分开显示";
  const tableHead=mode==="xic"?`<tr><th>形式</th><th>${escapeHtml(sampleShort(pair.reference))} XIC面积</th><th>${escapeHtml(sampleShort(pair.test))} XIC面积</th><th>XIC面积比</th><th>证据</th></tr>`:`<tr><th>形式</th><th>${escapeHtml(sampleShort(pair.reference))} 构成</th><th>${escapeHtml(sampleShort(pair.test))} 构成</th><th>构成变化/百分点</th><th>形式响应倍数</th><th>证据</th></tr>`;
  const featureButtons=(item.linked_feature_ids||[]).map(id=>`<button type="button" class="mod-quant-feature" data-mod-quant-feature="${escapeHtml(id)}">${escapeHtml(id)}</button>`).join("");
  return `<tr class="mod-quant-detail-row"><td colspan="7"><div class="mod-quant-detail"><div class="mod-quant-detail-grid"><div class="mod-composition-panel"><div class="mod-detail-title">${panelTitle}</div>${barForSample(pair.reference)}${barForSample(pair.test)}<div class="note">${panelNote}</div></div><div class="mod-form-panel"><div class="mod-detail-title">${tableTitle}</div><table class="mod-form-matrix">${tableHead}${matrixRows}</table></div></div><details class="mod-feature-details"><summary>关联 Feature（${(item.linked_feature_ids||[]).length}）</summary><div>${featureButtons||'<span class="note">无</span>'}</div></details></div></td></tr>`;
}
function modItemRows(item,pair,formPalette,nestedMember=false){
  const changes=item.currentChanges||[],strongest=item.currentStrongest;
  let conclusion='<span class="note">暂不能形成组间比例结论</span>';
  if(strongest){
    const higher=strongest.delta>=0?pair.test:pair.reference,higherLevel=strongest.delta>=0?strongest.testLevel:strongest.referenceLevel,lowerLevel=strongest.delta>=0?strongest.referenceLevel:strongest.testLevel;
    let responseText="";
    if(strongest.responseFold!=null){const responseHigher=strongest.responseFold>=1?pair.test:pair.reference,responseRatio=strongest.responseFold>=1?strongest.responseFold:1/strongest.responseFold;responseText=`；${escapeHtml(sampleShort(responseHigher))} 的该形式归一化响应为 ${responseRatio.toFixed(2)}x`;}
    conclusion=`<b>${escapeHtml(sampleShort(higher))}</b> 的“${escapeHtml(strongest.form.label)}”构成较高 <b>${(Math.abs(strongest.delta)*100).toFixed(2)} 个百分点</b><br><span class="note">${(lowerLevel*100).toFixed(2)}% → ${(higherLevel*100).toFixed(2)}%${responseText}</span>`;
  }
  const familyKey=String(item.quantitation_id||item.family_id||item.rank||""),expanded=state.expandedModificationFamilies.has(familyKey),standaloneHigh=item.high_value_candidate&&!nestedMember;
  const confidenceGrade=changes.some(change=>String(change.form.ms2_confidence||"").startsWith("C_"))?"含C级定位":"B级形式";
  const badges=`${standaloneHigh?'<span class="mod-quant-badge high-value">独立高价值候选</span>':""}${nestedMember?'<span class="mod-quant-badge">分组内序列</span>':""}${item.familyResponseImbalance?'<span class="mod-quant-badge high-value">该肽段总响应不平衡</span>':""}`;
  const rowClass=`${nestedMember?"mod-proteolytic-member-row ":""}${standaloneHigh?"high-value":""}`;
  const summaryRow=`<tr class="${rowClass}"><td><button type="button" class="mod-quant-toggle" data-mod-quant-toggle="${escapeHtml(familyKey)}" aria-expanded="${expanded}">${expanded?"−":"+"}</button></td><td>${item.rank||""}</td><td class="wrap"><div class="mod-quant-event">${escapeHtml(item.event_name||"")} ${badges}</div><div class="mod-quant-site">${escapeHtml(item.chain||"")} ${item.start||""}-${item.end||""}</div></td><td class="wrap"><b>${escapeHtml(item.sequence||"")}</b></td><td class="mod-quant-conclusion ${standaloneHigh?"high-value":""}">${conclusion}${item.familyResponseImbalance?'<br><span class="note">该肽段自身总响应不单独解释为修饰或蛋白丰度变化。</span>':""}</td><td class="wrap"><b>${escapeHtml(modificationQuantitationStatusLabel(item.quantitation_status))}</b><br><span class="note">${escapeHtml(confidenceGrade)}</span></td><td>${changes.length}</td></tr>`;
  return summaryRow+(expanded?modItemDetailRow(item,pair,formPalette):"");
}
function modProteolyticGroupRows(family,members,pair,formPalette){
  const key=String(family.family_id||family.rank||""),expanded=state.expandedProteolyticFamilies.has(key);
  const sites=(family.site_conclusions||[]).map(site=>({site,comparison:modOrientedSiteComparison(site,pair)})).filter(item=>item.comparison).sort((a,b)=>Math.abs(Number(b.comparison.median_percentage_point_difference||0))-Math.abs(Number(a.comparison.median_percentage_point_difference||0)));
  const lead=sites[0],comparison=lead?.comparison,site=lead?.site;
  const totals=members.map(member=>{const fold=member.familyReferenceArea>0&&member.familyTestArea>0?member.familyTestArea/member.familyReferenceArea:null;return {member,fold};}),folds=totals.map(item=>item.fold).filter(Boolean),redistribution=folds.some(fold=>fold>=1.2)&&folds.some(fold=>fold<=1/1.2);
  let conclusion='<span class="note">展开查看各条重叠肽段的定量结果。</span>';
  if(comparison?.consistent){
    const delta=Number(comparison.median_percentage_point_difference||0),higher=delta>=0?pair.test:pair.reference,magnitudes=[Math.abs(Number(comparison.min_percentage_point_difference||0)),Math.abs(Number(comparison.max_percentage_point_difference||0))].sort((a,b)=>a-b);
    conclusion=`<b>${escapeHtml(site.event_label)}</b>：${escapeHtml(sampleShort(higher))} 较高；${comparison.support_count}/${comparison.support_count} 条重叠肽段方向一致，变化范围 <b>${(magnitudes[0]*100).toFixed(2)}–${(magnitudes[1]*100).toFixed(2)} 个百分点</b>。`;
  }else if(comparison)conclusion=`<b>${escapeHtml(site.event_label)}</b>：重叠肽段间方向不一致，暂不形成位点级结论。`;
  const high=Boolean(site?.high_value_consensus&&comparison?.consistent&&comparison?.all_b_grade),rank=Math.min(...members.map(member=>Number(member.rank||9999))),sequences=members.map(member=>member.sequence).join(" / ");
  const groupRow=`<tr class="mod-proteolytic-group-row ${high?"high-value":""}"><td><button type="button" class="mod-quant-toggle" data-proteolytic-toggle="${escapeHtml(key)}" aria-expanded="${expanded}">${expanded?"−":"+"}</button></td><td>${rank}</td><td class="wrap"><div class="mod-quant-event">重叠肽段修饰组 ${high?'<span class="mod-quant-badge high-value">跨序列高价值</span>':""}${redistribution?'<span class="mod-quant-badge high-value">疑似酶切分配差异</span>':""}</div><div class="mod-quant-site">${escapeHtml(family.chain||"")} ${family.start||""}-${family.end||""}</div></td><td class="wrap"><b>${members.length} 条互相包含序列</b><div class="mod-group-sequences">${escapeHtml(sequences)}</div></td><td class="mod-quant-conclusion ${high?"high-value":""}">${conclusion}${redistribution?'<br><span class="note">长短肽总响应呈相反方向，展开后分别复核；不同序列面积不直接合并为绝对修饰率。</span>':""}</td><td class="wrap"><b>位点一致性判断</b><br><span class="note">各序列独立定量</span></td><td>${members.length} 序列</td></tr>`;
  if(!expanded)return groupRow;
  return groupRow+members.map(member=>modItemRows(member,pair,formPalette,true)).join("");
}
function modBestLinkedFeatureId(featureIds){
  return (featureIds||[]).map(id=>({id:String(id),item:featureItemById(id)})).filter(entry=>entry.item).sort((first,second)=>{
    const significant=item=>["presence_absence","area_changed","moderate_difference"].includes(String(item?.difference_type||""))?1:0;
    return significant(second.item)-significant(first.item)||Number(second.item.ranking_score||0)-Number(first.item.ranking_score||0)||Math.max(...Object.values(second.item.normalized_area_by_sample||second.item.area_by_sample||{}).map(Number),0)-Math.max(...Object.values(first.item.normalized_area_by_sample||first.item.area_by_sample||{}).map(Number),0);
  })[0]?.id||null;
}
function modFormContextByKey(formKey){
  for(const item of state.msms?.modification_level_quantitation||[]){
    for(const form of item.forms||[]){
      const key=`${item.quantitation_id||item.family_id||""}::${form.form_id||form.label||""}`;
      if(key===formKey)return {...form,quantitation_id:item.quantitation_id,family_id:item.family_id,sequence:item.sequence,chain:item.chain,start:item.start,end:item.end,reference_sample:item.reference_sample,form_label:form.label};
    }
  }
  return null;
}
async function selectTargetedModificationForm(context){
  if(!context)return;
  const bestPsm=context.best_psm||{},charges=[Number(bestPsm.precursor_charge||0),...Object.values(context.charge_states_by_sample||{}).flat().map(Number)].filter(charge=>Number.isFinite(charge)&&charge>0),charge=charges[0]||1,neutralMass=Number(context.neutral_mass),rts=(context.consensus_rts||[]).map(Number).filter(Number.isFinite),rt=rts.length?rts.sort((a,b)=>a-b)[Math.floor(rts.length/2)]:Number(bestPsm.rt),mz=neutralMass/charge+1.007276466621;
  if(!Number.isFinite(mz)||!Number.isFinite(rt))return;
  const peaks=DATA.peak_results||[],containing=peaks.find(peak=>Number(peak.rt_start)<=rt&&rt<=Number(peak.rt_end)),nearest=containing||[...peaks].sort((first,second)=>Math.abs(Number(first.rt_apex)-rt)-Math.abs(Number(second.rt_apex)-rt))[0];
  state.selectedPeakId=nearest?.tic_peak_id||null;state.selectedFeatureGroupId=null;state.scrollGlobalSelectionIntoView=false;state.detailZoom=null;state.xic=null;state.xicZoom=null;
  await selectSpectrumMz(mz,rt,null,false,true);
}
function renderModificationQuantitation(){
  const table=$("modQuantTable"),summary=$("modQuantSummary"),filter=$("modQuantFilter"),statusSelect=$("modQuantStatus"),compositionButton=$("modViewComposition"),xicButton=$("modViewXic");
  if(!table||!summary)return;
  if(filter&&!filter.dataset.bound){filter.dataset.bound="1";filter.oninput=renderModificationQuantitation;}
  if(statusSelect&&!statusSelect.dataset.bound){statusSelect.dataset.bound="1";statusSelect.onchange=renderModificationQuantitation;}
  if(compositionButton&&!compositionButton.dataset.bound){compositionButton.dataset.bound="1";compositionButton.onclick=()=>{state.modQuantViewMode="composition";renderModificationQuantitation();};}
  if(xicButton&&!xicButton.dataset.bound){xicButton.dataset.bound="1";xicButton.onclick=()=>{state.modQuantViewMode="xic";renderModificationQuantitation();};}
  if(compositionButton)compositionButton.classList.toggle("active",state.modQuantViewMode==="composition");
  if(xicButton)xicButton.classList.toggle("active",state.modQuantViewMode==="xic");
  if(state.msmsLoading||!state.msmsLoaded){summary.innerHTML='<span class="mod-quant-badge">正在载入修饰形式定量…</span>';table.innerHTML="";return;}
  if(state.msmsError){summary.innerHTML=`<span class="mod-quant-badge">不可用</span><span class="note">${escapeHtml(state.msmsError)}</span>`;table.innerHTML="";return;}
  const all=Array.isArray(state.msms?.modification_level_quantitation)?state.msms.modification_level_quantitation:[],pair=selectedPair(),query=String(filter?.value||"").trim().toLowerCase(),filterMode=statusSelect?.value||"quantifiable";
  const evaluated=all.map(item=>modEvaluateItem(item,pair)).filter(item=>{
    if(filterMode==="formal"&&item.quantitation_status!=="formal_relative_quantitation")return false;
    if(filterMode==="quantifiable"&&!item.quantifiable)return false;
    if(!query)return true;
    return [item.event_name,item.event_type,item.sequence,item.chain,...(item.linked_feature_ids||[]),...(item.forms||[]).map(form=>form.label)].join(" ").toLowerCase().includes(query);
  }).sort((a,b)=>Math.abs(Number(b.currentStrongest?.delta||0))-Math.abs(Number(a.currentStrongest?.delta||0))||Number(b.high_value_candidate)-Number(a.high_value_candidate)||Number(a.rank)-Number(b.rank));
  const byId=new Map(evaluated.map(item=>[String(item.quantitation_id||""),item])),consumed=new Set(),entries=[];
  for(const family of state.msms?.proteolytic_modification_families||[]){
    const members=(family.members||[]).map(member=>byId.get(String(member.quantitation_id||""))).filter(Boolean);
    if(members.length<2)continue;
    members.forEach(member=>consumed.add(String(member.quantitation_id||"")));
    entries.push({kind:"group",family,members,sortRank:Math.min(...members.map(member=>Number(member.rank||9999))),sortScore:Math.max(...members.map(member=>Math.abs(Number(member.currentStrongest?.delta||0))))});
  }
  evaluated.filter(item=>!consumed.has(String(item.quantitation_id||""))).forEach(item=>entries.push({kind:"item",item,sortRank:Number(item.rank||9999),sortScore:Math.abs(Number(item.currentStrongest?.delta||0))}));
  entries.sort((a,b)=>b.sortScore-a.sortScore||a.sortRank-b.sortRank);
  const formal=all.filter(item=>item.quantitation_status==="formal_relative_quantitation").length,siteHigh=(state.msms?.proteolytic_modification_families||[]).filter(item=>item.high_value_consensus).length,standaloneHigh=all.filter(item=>item.high_value_candidate&&!Array.from(consumed).includes(String(item.quantitation_id||""))).length;
  summary.innerHTML=`<span class="mod-quant-badge">${entries.length} 个表格分组（${all.length} 个序列形式族）</span><span class="mod-quant-badge formal">${formal} 个正式相对定量</span><span class="mod-quant-badge high-value">${siteHigh} 个跨序列高价值结论</span>${standaloneHigh?`<span class="mod-quant-badge high-value">${standaloneHigh} 个独立高价值候选</span>`:""}<span class="note">当前图表：${state.modQuantViewMode==="xic"?"归一化 XIC 面积":"组成比例"}；${escapeHtml(sampleShort(pair.test))} 对 ${escapeHtml(sampleShort(pair.reference))}</span>`;
  if(!entries.length){table.innerHTML='<tr><td class="note">当前筛选条件下没有可显示的修饰形式族。</td></tr>';return;}
  const formPalette=["#2563eb","#f59e0b","#10b981","#8b5cf6","#ef4444","#06b6d4","#84cc16","#f97316","#64748b"];
  table.innerHTML='<tr><th></th><th>rank</th><th>事件 / 位置</th><th>肽段 / 分组</th><th>主要变化</th><th>定量等级</th><th>形式数</th></tr>'+entries.map(entry=>entry.kind==="group"?modProteolyticGroupRows(entry.family,entry.members,pair,formPalette):modItemRows(entry.item,pair,formPalette,false)).join("");
  table.querySelectorAll("[data-proteolytic-toggle]").forEach(button=>button.onclick=()=>{const key=String(button.dataset.proteolyticToggle||"");if(state.expandedProteolyticFamilies.has(key))state.expandedProteolyticFamilies.delete(key);else state.expandedProteolyticFamilies.add(key);renderModificationQuantitation();});
  table.querySelectorAll("[data-mod-quant-toggle]").forEach(button=>button.onclick=()=>{const key=String(button.dataset.modQuantToggle||"");if(state.expandedModificationFamilies.has(key))state.expandedModificationFamilies.delete(key);else state.expandedModificationFamilies.add(key);renderModificationQuantitation();});
  table.querySelectorAll("[data-mod-form-key]").forEach(button=>button.onclick=async event=>{event.stopPropagation();const formKey=String(button.dataset.modFormKey||""),context=modFormContextByKey(formKey),ids=String(button.dataset.modFormLinks||"").split("|").filter(Boolean),featureId=modBestLinkedFeatureId(ids);state.selectedModificationFormKey=formKey;state.selectedModificationFormContext=context;renderModificationQuantitation();if(featureId)await selectMappedFeature(featureId,true);else await selectTargetedModificationForm(context);});
  table.querySelectorAll("[data-mod-quant-feature]").forEach(button=>button.onclick=()=>selectMappedFeature(button.dataset.modQuantFeature));
  requestAnimationFrame(()=>requestAnimationFrame(()=>table.querySelectorAll(".mod-composition-scale[data-target-width]").forEach(scale=>{scale.style.width=`${Math.max(0,Math.min(100,Number(scale.dataset.targetWidth||0)))}%`;scale.classList.add("animated");})));
}
async function ensureMsmsData(){
  if(state.msmsLoaded||state.msmsLoading)return state.msmsPromise;
  const configuredStatus=String(DATA?.analysis_metadata?.ms2?.status||"");
  if(DATA?.analysis_metadata&&configuredStatus!=="completed"){
    state.msmsError=String(DATA?.analysis_metadata?.ms2?.reason||"未提供 FASTA，未执行 MS2 计算。");
    state.msmsLoaded=true;
    drawFeatureMs2(); renderModificationQuantitation(); drawStructureModule(); maybeLoadProjectStructure();
    return Promise.resolve();
  }
  const requestComparison=state.comparison;
  state.msmsLoading=true; state.msmsError=null; drawFeatureMs2(); renderModificationQuantitation(); drawStructureModule();
  state.msmsPromise=fetchJson(apiUrl("/api/msms")).then(payload=>{
    if(state.comparison!==requestComparison)return;
    state.msms=payload; state.msmsLoaded=true;
  }).catch(err=>{
    if(state.comparison!==requestComparison)return;
    state.msmsError=String(err?.message||err); state.msmsLoaded=true;
  }).finally(()=>{
    if(state.comparison!==requestComparison)return;
    state.msmsLoading=false; state.msmsPromise=null; drawFeatureMs2(); renderModificationQuantitation(); drawStructureModule(); maybeLoadProjectStructure();
    renderGlobalFeatureTable(); drawFeatureMap();
  });
  return state.msmsPromise;
}
function renderAgentAnnotations(feature=null){
  const container=$("featureMs2Annotations");
  if(!container)return;
  const featureIds=new Set([String(state.selectedFeatureGroupId||""),String(feature?.feature_group_id||"")].filter(Boolean));
  const annotations=(state.agentAnnotations||[]).filter(item=>featureIds.has(String(item?.feature_group_id||"")));
  if(!annotations.length){ container.hidden=true; container.innerHTML=""; return; }
  const rows=annotations.map(item=>{
    const agent=escapeHtml(item.agent_name||"MCP agent"), created=escapeHtml(item.created_at||"");
    const suggestion=escapeHtml(item.suggestion||""), rationale=escapeHtml(item.rationale||"");
    return `<div class="agent-annotation"><b>待人工确认</b> · ${agent}${created?` · ${created}`:""}<br>建议：${suggestion}${rationale?`<br>依据：${rationale}`:""}</div>`;
  });
  container.hidden=false;
  container.innerHTML=`<b>Agent 待人工确认注释（${annotations.length}）</b>${rows.join("")}`;
}
function drawFeatureMs2(){
  const canvas=$("featureMs2Canvas"),detail=$("featureMs2Detail"),info=$("featureMs2Info");
  if(!canvas||!detail||!info)return;
  renderAgentAnnotations();
  const ctx=canvas.getContext("2d"); ctx.clearRect(0,0,canvas.width,canvas.height);
  canvas._featureMs2Peaks=[]; canvas._featureMs2Domain=null;
  const featureId=String(state.selectedFeatureGroupId||"");
  const formContext=state.selectedModificationFormContext;
  if(!featureId&&!formContext){
    detail.textContent="Select a Feature from the heatmap or Feature list.";
    info.textContent="The corresponding best identified spectrum, or the best covering unresolved MS2 scan, will appear here.";
    return;
  }
  if(state.msmsLoading){
    detail.textContent=`${featureId} | loading corresponding MS/MS evidence...`;
    info.textContent="";
    return;
  }
  if(state.msmsError){
    const skipped=String(DATA?.analysis_metadata?.ms2?.status||"").startsWith("skipped_");
    detail.textContent=`${featureId} | ${skipped?"MS/MS skipped":"MS/MS report unavailable"}`;
    info.textContent=state.msmsError;
    return;
  }
  if(!state.msmsLoaded){
    detail.textContent=`${featureId} | preparing corresponding MS/MS evidence...`;
    info.textContent="";
    return;
  }
  const matched=featureId?selectedMs2Evidence():(formContext?.best_psm?{feature:{feature_group_id:"",best_psm:formContext.best_psm,sequence:formContext.sequence,chain:formContext.chain,start:formContext.start,end:formContext.end,modification:formContext.form_label,confidence:formContext.ms2_confidence,ms2_status:formContext.status},relation:"第9节目标形式"}:null);
  if(!matched){
    if(formContext){detail.innerHTML=`<span class="feature-ms2-status unresolved">目标 XIC 形式</span><b>${escapeHtml(formContext.form_label||"")}</b><br>sequence <b>${escapeHtml(formContext.sequence||"")}</b>`;info.textContent="该形式已联动目标 XIC 和序列/结构定位，但没有可显示的代表 MS/MS 谱图，也未进入差异 Feature 表。";}
    else {detail.innerHTML=`<span class="feature-ms2-status missing">MS2 evidence index unavailable</span>${escapeHtml(featureId)}`;info.textContent="该 Feature 未出现在当前报告的 MS2 关联索引中，通常表示这是尚未重新计算的旧报告，而不是已经确认未采集 MS2。重新运行 MS2 后会区分“未采集 MS2”和“已采集但未定性”。";}
    return;
  }
  const f=matched.feature, evidenceP=f.best_psm||f.candidate_psm||f.coverage_scan;
  const spectrumP=[f.best_psm,f.candidate_psm,f.coverage_scan].find(item=>Array.isArray(item?.spectrum_peaks)&&item.spectrum_peaks.length)||evidenceP;
  const confidence=String(f.confidence||""), identified=Boolean(f.best_psm)&&confidence.startsWith("B_"), candidate=Boolean(f.candidate_psm)||confidence.startsWith("C_")||confidence.startsWith("D_");
  // Older task databases labelled some D-level sequence matches as identified.
  // Derive the display status from confidence so reopening an old report cannot
  // overstate a low-evidence candidate.
  const effectiveStatus=confidence.startsWith("D_")?"low_evidence_sequence_candidate":f.ms2_status;
  const status=ms2StatusText[effectiveStatus]||effectiveStatus||"unknown MS2 status";
  const statusClass=identified?"identified":((candidate||evidenceP)?"unresolved":"missing");
  const sequence=f.sequence||evidenceP?.sequence||"", modification=f.modification||evidenceP?.modification_text||"";
  const displayId=featureId||(formContext?.form_label||"目标 XIC 形式"),relation=!featureId||String(f.feature_group_id||"")===featureId?"":` | evidence Feature ${f.feature_group_id} (${matched.relation})`;
  detail.innerHTML=`<span class="feature-ms2-status ${statusClass}">${escapeHtml(status)}</span><b>${escapeHtml(displayId)}</b>${escapeHtml(relation)}<br>${sequence?`sequence <b>${escapeHtml(sequence)}</b>${modification?` | ${escapeHtml(modification)}`:""}`:"No peptide sequence assigned"}`;
  const reason=f.unidentified_reason||f.hypothesis||"";
  const scanText=spectrumP?`${spectrumP.sample_id||""} | scan ${spectrumP.scan_id||""} | RT ${nice(spectrumP.rt,3)} min | precursor ${nice(spectrumP.precursor_mz,5)} z${spectrumP.precursor_charge||"?"}`:"No covering MS2 scan";
  const glycanText=evidenceP?.glycan_name?` | glycan diagnostic ions ${evidenceP.glycan_diagnostic_ion_count||0} | core Y ions ${evidenceP.glycan_core_y_ion_count||0} | HexNAc-retaining fragments ${evidenceP.glycan_hexnac_fragment_count||0}`:"";
  const scoreText=(identified||candidate)&&evidenceP?` | ${candidate&&!identified?"candidate ":""}score ${nice(evidenceP.score,1)} | q ${nice(evidenceP.q_value,4)} | matched ions ${evidenceP.matched_ion_count||0} | fragment coverage ${nice(evidenceP.fragment_coverage,2)}${evidenceP.sequence_tag_length!=null?` | sequence tag ${evidenceP.sequence_tag_length} aa`:""}${evidenceP.mass_delta!=null?` | ΔMass ${nice(evidenceP.mass_delta,4)} Da`:""}${glycanText}`:"";
  const fallbackText=spectrumP&&evidenceP&&spectrumP!==evidenceP?" | displaying the best stored covering scan":"";
  info.textContent=`${scanText}${scoreText}${fallbackText}${reason?` | ${reason}`:""}`;
  if(featureId)renderAgentAnnotations(f);
  const peaks=(spectrumP?.spectrum_peaks||[]).map(v=>({mz:Number(v.mz),intensity:Number(v.intensity),label:String(v.label||"")})).filter(v=>Number.isFinite(v.mz)&&Number.isFinite(v.intensity)&&v.intensity>0);
  canvas._featureMs2Peaks=peaks;
  if(!peaks.length){
    ctx.fillStyle="#667085"; ctx.font="14px Arial"; ctx.textAlign="center";
    ctx.fillText(evidenceP?"No displayable fragment peaks are stored for this evidence.":"No MS2 scan was acquired for this Feature.",canvas.width/2,canvas.height/2);
    ctx.textAlign="left";
    return;
  }
  const pad={l:62,r:18,t:36,b:48}, minMz=0, maxMz=Math.max(...peaks.map(v=>v.mz))*1.03, maxI=Math.max(...peaks.map(v=>v.intensity));
  const px=mz=>pad.l+(mz-minMz)/Math.max(maxMz-minMz,1e-9)*(canvas.width-pad.l-pad.r);
  const py=intensity=>canvas.height-pad.b-intensity/maxI*(canvas.height-pad.t-pad.b);
  canvas._featureMs2Domain={pad,minMz,maxMz,maxI,px,py};
  ctx.strokeStyle="#98a2b3"; ctx.lineWidth=1; ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,canvas.height-pad.b); ctx.lineTo(canvas.width-pad.r,canvas.height-pad.b); ctx.stroke();
  ctx.fillStyle="#667085"; ctx.font="11px Arial";
  for(let i=0;i<=5;i++){
    const mz=maxMz*i/5,x=px(mz); ctx.textAlign=i===0?"left":(i===5?"right":"center"); ctx.fillText(axisValue(mz),x,canvas.height-24);
    ctx.strokeStyle="#e4e7ec"; ctx.beginPath(); ctx.moveTo(x,pad.t); ctx.lineTo(x,canvas.height-pad.b); ctx.stroke();
  }
  [0,25,50,75,100].forEach(percent=>{
    const y=canvas.height-pad.b-percent/100*(canvas.height-pad.t-pad.b);
    ctx.fillStyle="#667085"; ctx.textAlign="right"; ctx.fillText(`${percent}%`,pad.l-10,y+4);
  });
  peaks.forEach(v=>{
    const x=px(v.mz), y=py(v.intensity), isB=/^b/i.test(v.label), isY=/^y/i.test(v.label);
    ctx.strokeStyle=isB?"#2563eb":(isY?"#dc2626":(v.label?"#7c3aed":"#667085"));
    ctx.lineWidth=v.label?1.8:1;
    ctx.beginPath(); ctx.moveTo(x,canvas.height-pad.b); ctx.lineTo(x,y); ctx.stroke();
  });
  const labelled=peaks.filter(v=>v.label).sort((a,b)=>b.intensity-a.intensity);
  labelled.slice(0,24).forEach((v,index)=>{
    const x=px(v.mz), y=py(v.intensity), isB=/^b/i.test(v.label), isY=/^y/i.test(v.label);
    ctx.fillStyle=isB?"#1d4ed8":(isY?"#b42318":"#6d28d9"); ctx.font="11px Arial";
    ctx.fillText(v.label,Math.min(canvas.width-pad.r-24,x+2),Math.max(pad.t+10,y-3-(index%2)*10));
  });
  ctx.fillStyle="#475467"; ctx.font="12px Arial"; ctx.textAlign="center"; ctx.fillText("fragment m/z",(pad.l+canvas.width-pad.r)/2,canvas.height-7);
  ctx.textAlign="left"; ctx.fillText("relative intensity",pad.l+4,16);
  ctx.fillStyle="#2563eb"; ctx.fillText("b ions",canvas.width-168,18);
  ctx.fillStyle="#dc2626"; ctx.fillText("y ions",canvas.width-112,18);
  ctx.fillStyle="#667085"; ctx.fillText("unmatched",canvas.width-62,18);
}
function featureMs2HitAt(canvas,clientX,clientY){
  const domain=canvas._featureMs2Domain, peaks=canvas._featureMs2Peaks||[]; if(!domain||!peaks.length)return null;
  const rect=canvas.getBoundingClientRect(), x=(clientX-rect.left)*canvas.width/rect.width, y=(clientY-rect.top)*canvas.height/rect.height;
  if(x<domain.pad.l||x>canvas.width-domain.pad.r||y<domain.pad.t||y>canvas.height-domain.pad.b)return null;
  let best=null;
  peaks.forEach(peak=>{ const px=domain.px(peak.mz), py=domain.py(peak.intensity), dx=Math.abs(px-x); if(dx>10||y<py-10||y>canvas.height-domain.pad.b+8)return; const distance=Math.hypot(dx,Math.max(0,py-y)); if(!best||distance<best.distance)best={...peak,px,py,distance}; });
  return best;
}
function showFeatureMs2Tooltip(event,hit){
  const tip=$("detailTooltip"); if(!tip||!hit){hideDetailTooltip();return;}
  const label=hit.label?`<span class="ion-badge">${escapeHtml(hit.label)}</span>`:`<span class="tooltip-muted">未匹配到推断离子</span>`;
  tip.innerHTML=`<div class="tooltip-title">MS/MS fragment</div>${label}<div>m/z <b>${nice(hit.mz,5)}</b></div><div>信号强度 <b>${axisValue(hit.intensity)}</b></div>`;
  positionTooltip(tip,event); tip.style.display="block";
}
function pairMetrics(item){
  const {reference,test}=selectedPair(), values=featureValueMap(item);
  const referenceValue=Math.max(0,Number(values?.[reference])||0), testValue=Math.max(0,Number(values?.[test])||0), maxValue=Math.max(referenceValue,testValue);
  if(!reference||!test||maxValue<=0)return {reference,test,referenceValue,testValue,logRatio:0,fold:1,direction:"equal"};
  const zeroFloor=Math.max(maxValue*1e-6,1e-12);
  const rawLogRatio=Math.log(Math.max(testValue,zeroFloor)/Math.max(referenceValue,zeroFloor));
  const logRatio=Math.max(-MAX_FEATURE_LN_FOLD,Math.min(MAX_FEATURE_LN_FOLD,rawLogRatio));
  return {reference,test,referenceValue,testValue,rawLogRatio,logRatio,fold:Math.exp(Math.abs(logRatio)),direction:rawLogRatio>1e-9?"higher":(rawLogRatio<-1e-9?"lower":"equal")};
}
function pairDifferenceType(item){
  const pair=pairMetrics(item), presenceFraction=Number(DATA?.params?.feature_presence_relative_area_fraction)||0.02;
  const commonFold=Number(DATA?.params?.feature_common_fold_change_threshold)||2, strongFold=Number(DATA?.params?.feature_strong_fold_change_threshold)||4;
  if(Math.max(pair.referenceValue,pair.testValue)<=0)return "low_confidence";
  if(Math.min(pair.referenceValue,pair.testValue)/Math.max(pair.referenceValue,pair.testValue)<=presenceFraction)return "presence_absence";
  if(pair.fold>=strongFold)return "area_changed";
  if(pair.fold>=commonFold)return "moderate_difference";
  return "common_feature";
}
function abundanceOrderHtml(item){
  const pair=pairMetrics(item), entries=[{sample:pair.reference,value:pair.referenceValue},{sample:pair.test,value:pair.testValue}].filter(entry=>entry.sample).sort((a,b)=>b.value-a.value);
  const title=entries.map(entry=>`${entry.sample}: ${nice(entry.value,1)}`).join(" > ");
  return `<span class="sample-order" title="${escapeHtml(title)}">${entries.map(entry=>`<span class="sample-dot" style="background:${sampleColor(entry.sample)}" title="${escapeHtml(entry.sample)}: ${nice(entry.value,1)}"></span>`).join("")}</span>`;
}
function pairDirectionHtml(item){
  const pair=pairMetrics(item), fold=nice(pair.fold,2);
  if(pair.direction==="higher")return `<span class="direction-up" title="Test / reference">&#8593; ${escapeHtml(sampleShort(pair.test))} ${fold}x</span>`;
  if(pair.direction==="lower")return `<span class="direction-down" title="Test / reference">&#8595; ${escapeHtml(sampleShort(pair.test))} ${fold}x</span>`;
  return `<span class="direction-equal">&#8776; 1.00x</span>`;
}
function sortValue(item,key){
  if(key==="pair_log_ratio")return pairMetrics(item).logRatio;
  if(key==="pair_fold")return pairMetrics(item).fold;
  if(key==="pair_type")return pairDifferenceType(item);
  return item?.[key];
}
function sortedTableRows(rows,sort){
  const direction=sort?.dir==="asc"?1:-1;
  return [...(rows||[])].sort((a,b)=>{
    const av=sortValue(a,sort?.key), bv=sortValue(b,sort?.key);
    if(typeof av==="number"||typeof bv==="number")return ((Number(av)||0)-(Number(bv)||0))*direction;
    return String(av??"").localeCompare(String(bv??""))*direction;
  });
}
function sortableTh(label,key,scope){
  const sort=scope==="global"?state.globalSort:state.mzSort, arrow=sort.key===key?(sort.dir==="asc"?" &#9650;":" &#9660;"):"";
  return `<th class="sortable" data-sort-scope="${scope}" data-sort-key="${key}">${label}${arrow}</th>`;
}
function bindSortableHeaders(scope,render){
  document.querySelectorAll(`th[data-sort-scope="${scope}"]`).forEach(header=>header.onclick=event=>{
    event.stopPropagation(); const sort=scope==="global"?state.globalSort:state.mzSort, key=header.dataset.sortKey;
    if(sort.key===key)sort.dir=sort.dir==="asc"?"desc":"asc"; else {sort.key=key;sort.dir=["feature_group_id","component_label","difference_type","pair_type","quantitation_confidence"].includes(key)?"asc":"desc";}
    render();
  });
}
function updatePairLegend(){
  const pair=selectedPair(), legend=$("featureMapPairLegend"); if(!legend)return;
  legend.innerHTML=`<span><span class="sample-dot" style="background:${sampleColor(pair.reference)}"></span>Reference: ${escapeHtml(pair.reference||"")}</span><span><span class="sample-dot" style="background:${sampleColor(pair.test)}"></span>Test: ${escapeHtml(pair.test||"")}</span>`;
}
function redrawPairViews(){
  hideDetailTooltip(); updatePairLegend(); renderMzTable(); renderGlobalFeatureTable(); drawFeatureMap(); renderModificationQuantitation(); drawStructureModule();
}
function renderPairControls(reset=false){
  const samples=DATA?.sample_ids||[], referenceSelect=$("featureMapReference"), testSelect=$("featureMapTest"), swap=$("swapFeatureMapSamples"), mode=$("featureMapMode");
  if(!samples.length)return;
  if(reset||!samples.includes(state.pairReference))state.pairReference=samples.includes(DATA.reference_sample)?DATA.reference_sample:samples[0];
  if(reset||!samples.includes(state.pairTest)||state.pairTest===state.pairReference)state.pairTest=samples.find(sample=>sample!==state.pairReference)||state.pairReference;
  const options=samples.map(sample=>`<option value="${escapeHtml(sample)}">${escapeHtml(sample)}</option>`).join("");
  referenceSelect.innerHTML=options; testSelect.innerHTML=options;
  referenceSelect.value=state.pairReference; testSelect.value=state.pairTest;
  referenceSelect.disabled=samples.length<2; testSelect.disabled=samples.length<2; swap.disabled=samples.length<2;
  mode.value=state.featureMapMode;
  referenceSelect.onchange=()=>{
    const oldReference=state.pairReference; state.pairReference=referenceSelect.value;
    if(state.pairReference===state.pairTest)state.pairTest=samples.includes(oldReference)&&oldReference!==state.pairReference?oldReference:(samples.find(sample=>sample!==state.pairReference)||state.pairReference);
    testSelect.value=state.pairTest; redrawPairViews();
  };
  testSelect.onchange=()=>{
    state.pairTest=testSelect.value;
    if(state.pairTest===state.pairReference)state.pairReference=samples.find(sample=>sample!==state.pairTest)||state.pairTest;
    referenceSelect.value=state.pairReference; redrawPairViews();
  };
  swap.onclick=()=>{
    const reference=state.pairReference; state.pairReference=state.pairTest; state.pairTest=reference;
    referenceSelect.value=state.pairReference; testSelect.value=state.pairTest; redrawPairViews();
  };
  mode.onchange=()=>{state.featureMapMode=mode.value==="magnitude"?"magnitude":"direction"; hideDetailTooltip(); drawFeatureMap();};
  updatePairLegend();
}
function plot(canvas, domain){ const p={left:66,right:canvas.width-24,top:18,bottom:canvas.height-48,...domain}; p.toX=x=>p.left+(x-p.xmin)/Math.max(p.xmax-p.xmin,1e-12)*(p.right-p.left); p.toY=y=>p.bottom-(y-p.ymin)/Math.max(p.ymax-p.ymin,1e-12)*(p.bottom-p.top); p.fromX=x=>p.xmin+(x-p.left)/Math.max(p.right-p.left,1e-12)*(p.xmax-p.xmin); p.fromY=y=>p.ymin+(p.bottom-y)/Math.max(p.bottom-p.top,1e-12)*(p.ymax-p.ymin); return p; }
function axisValue(value){
  const n=Number(value); if(!Number.isFinite(n))return "";
  const abs=Math.abs(n), units=[[1e9,"B"],[1e6,"M"],[1e3,"k"]];
  for(const [scale,suffix] of units)if(abs>=scale)return `${(n/scale).toFixed(abs/scale>=100?0:2)}${suffix}`;
  if(abs>=100)return n.toFixed(0);
  if(abs>=1)return n.toFixed(2);
  return n.toFixed(3);
}
function axes(ctx,p,xl,yl){
  ctx.save();
  ctx.strokeStyle="#d9e0ea"; ctx.lineWidth=1; ctx.beginPath(); ctx.moveTo(p.left,p.top); ctx.lineTo(p.left,p.bottom); ctx.lineTo(p.right,p.bottom); ctx.stroke();
  ctx.fillStyle="#667085"; ctx.font="11px Inter,Arial,sans-serif";
  ctx.save(); ctx.translate(15,(p.top+p.bottom)/2); ctx.rotate(-Math.PI/2); ctx.textAlign="center"; ctx.fillText(yl,0,0); ctx.restore();
  ctx.textBaseline="top";
  for(let i=0;i<=5;i++){
    const x=p.left+(p.right-p.left)*i/5, y=p.bottom-(p.bottom-p.top)*i/5;
    ctx.textAlign=i===0?"left":(i===5?"right":"center"); ctx.fillText(axisValue(p.xmin+(p.xmax-p.xmin)*i/5),x,p.bottom+9);
    ctx.textAlign="right"; ctx.fillText(axisValue(p.ymin+(p.ymax-p.ymin)*i/5),p.left-10,y-5);
  }
  ctx.textAlign="center"; ctx.textBaseline="alphabetic"; ctx.font="12px Inter,Arial,sans-serif"; ctx.fillText(xl,(p.left+p.right)/2,p.bottom+34);
  ctx.restore();
}
function xDomainWidth(z){ return Math.max((z?.xmax||0)-(z?.xmin||0),1e-9); }
function setPan(id, zoom, full, onPan){
  const el=$(id); if(!el||!full){return;}
  const width=xDomainWidth(zoom); const fullWidth=Math.max(full.xmax-full.xmin,1e-9);
  if(!zoom||width>=fullWidth*.999){ el.disabled=true; el.min=full.xmin; el.max=full.xmin; el.value=full.xmin; return; }
  el.disabled=false; el.min=full.xmin; el.max=full.xmax-width; el.step=Math.max(fullWidth/1000,1e-6); el.value=Math.max(full.xmin,Math.min(zoom.xmin,full.xmax-width));
  el.oninput=()=>onPan(Number(el.value),width);
}
function drawDragBox(canvas,start,end){
  const ctx=canvas.getContext("2d"); const x=Math.min(start.x,end.x), y=Math.min(start.y,end.y), w=Math.abs(end.x-start.x), h=Math.abs(end.y-start.y);
  ctx.save(); ctx.fillStyle="rgba(37,99,235,.16)"; ctx.strokeStyle="rgba(37,99,235,.75)"; ctx.setLineDash([6,4]); ctx.fillRect(x,y,w,h); ctx.strokeRect(x,y,w,h); ctx.restore();
}
function positionTooltip(tip,event){
  const width=Math.min(360,Math.max(220,window.innerWidth-24)), left=Math.max(8,Math.min(window.innerWidth-width-8,event.clientX+14)), top=Math.max(8,Math.min(window.innerHeight-180,event.clientY+14));
  tip.style.left=`${left}px`; tip.style.top=`${top}px`;
}
function hideDetailTooltip(){ const tip=$("detailTooltip"); if(tip)tip.style.display="none"; }
function chromSeries(){
  const aligned=$("showAligned").checked;
  const pk=selectedPeak();
  const localMode=localPeakAlignmentEnabled();
  return DATA.sample_ids.map((sample,i)=>{
    const pts=DATA.chromatograms[sample][aligned?"aligned":"raw"];
    const localShift=aligned&&localMode&&pk?Number(((pk.local_alignment||{}).local_peak_shift_by_sample||{})[sample]||0):0;
    return {sample,color:color(i),localShift,x:pts.map(p=>p[0]+localShift),y:pts.map(p=>$("chromMode").value==="bpc"?p[2]:p[1])};
  });
}
function localAlignmentGatePassed(pk){ return (pk?.local_alignment||{}).local_alignment_quality_gate_passed!==false; }
function localPeakAlignmentEnabled(){ const pk=selectedPeak(); return !!($("showAligned")?.checked && $("showLocalPeakAligned")?.checked && pk && localAlignmentGatePassed(pk)); }
function localAlignmentStatus(pk){ const local=pk?.local_alignment||{}; return local.local_alignment_status||(local.local_alignment_quality_gate_passed===false?"skipped_no_clear_peak":"not_needed"); }
function samplePeakBand(pk, sample){
  const aligned=$("showAligned").checked;
  const localMode=localPeakAlignmentEnabled();
  const bounds=(pk.sample_cut_bounds||{})[sample]||{};
  const globalShift=Number(((DATA.alignment||{}).rt_shift_by_sample||{})[sample]||0);
  if(aligned){
    if(localMode) return {start:Number(bounds.aligned_rt_start ?? pk.rt_start), end:Number(bounds.aligned_rt_end ?? pk.rt_end)};
    return {start:Number(bounds.sample_corrected_rt_start ?? pk.rt_start), end:Number(bounds.sample_corrected_rt_end ?? pk.rt_end)};
  }
  return {
    start:Number(bounds.sample_corrected_rt_start ?? pk.rt_start)-globalShift,
    end:Number(bounds.sample_corrected_rt_end ?? pk.rt_end)-globalShift,
  };
}
function drawSelectedPeakBand(ctx,p,pk,series){
  if(!pk)return;
  const bands=series.map(s=>{
    const band=samplePeakBand(pk,s.sample);
    return {start:Number(band.start),end:Number(band.end)};
  }).filter(band=>Number.isFinite(band.start)&&Number.isFinite(band.end));
  if(!bands.length)return;
  const x0=p.toX(Math.min(...bands.map(band=>Math.min(band.start,band.end))));
  const x1=p.toX(Math.max(...bands.map(band=>Math.max(band.start,band.end))));
  if(x1>=p.left&&x0<=p.right){
    const left=Math.max(p.left,Math.min(x0,x1)), right=Math.min(p.right,Math.max(x0,x1));
    ctx.save();
    ctx.fillStyle="rgba(59,130,246,.11)";
    ctx.fillRect(left,p.top,Math.max(1,right-left),p.bottom-p.top);
    ctx.strokeStyle="rgba(37,99,235,.48)";
    ctx.strokeRect(left,p.top,Math.max(1,right-left),p.bottom-p.top);
    ctx.restore();
  }
}
function selectedFeatureRt(){
  if(state.selectedFeatureRt!==null&&state.selectedFeatureRt!==undefined&&state.selectedFeatureRt!==""){
    const explicit=Number(state.selectedFeatureRt);
    if(Number.isFinite(explicit))return explicit;
  }
  const group=selectedFeatureGroup(), inferred=Number(group?.representative_rt);
  return Number.isFinite(inferred)?inferred:null;
}
function selectedFeatureDisplayRt(){
  const rt=selectedFeatureRt();
  if(!Number.isFinite(rt))return null;
  if($("showAligned")?.checked)return rt;
  const reference=DATA?.reference_sample||DATA?.sample_ids?.[0];
  const shift=Number(((DATA?.alignment||{}).rt_shift_by_sample||{})[reference]||0);
  return rt-shift;
}
function drawSelectedFeatureRt(ctx,p){
  const featureRt=selectedFeatureDisplayRt();
  if(!Number.isFinite(featureRt))return;
  const x=p.toX(featureRt);
  if(x<p.left||x>p.right)return;
  const label=`Feature RT ${nice(selectedFeatureRt(),4)} min`;
  ctx.save();
  ctx.strokeStyle="#7c3aed";
  ctx.lineWidth=2;
  ctx.setLineDash([8,5]);
  ctx.beginPath();
  ctx.moveTo(x,p.top);
  ctx.lineTo(x,p.bottom);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.font="bold 12px Arial";
  const textWidth=ctx.measureText(label).width;
  const labelX=Math.max(p.left+3,Math.min(x+7,p.right-textWidth-9));
  ctx.fillStyle="rgba(255,255,255,.90)";
  ctx.fillRect(labelX-3,p.top+4,textWidth+6,18);
  ctx.fillStyle="#6d28d9";
  ctx.fillText(label,labelX,p.top+17);
  ctx.restore();
}
function chromFullDomain(){
  const series=chromSeries(); const xs=series.flatMap(s=>s.x), ys=series.flatMap(s=>s.y);
  return {series,xmin:Math.min(...xs),xmax:Math.max(...xs),ymin:0,ymax:Math.max(...ys,1)};
}
function updatePlotLegend(containerId,series,visible=true){
  const container=$(containerId); if(!container)return;
  container.hidden=!visible;
  container.innerHTML=(series||[]).map(item=>`<span class="chrom-legend-item"><i class="chrom-legend-swatch" style="background:${item.color}"></i>${escapeHtml(item.sample)}</span>`).join("");
}
function updateChromLegend(series){ updatePlotLegend("chromLegend",series,reportSettings.legends.chrom); }
function updateSpectrumLegend(){ updatePlotLegend("spectrumLegend",(DATA?.sample_ids||[]).map((sample,index)=>({sample,color:color(index)})),reportSettings.legends.spectrum); }
function updateXicLegend(){ updatePlotLegend("xicLegend",(DATA?.sample_ids||[]).map((sample,index)=>({sample,color:color(index)})),reportSettings.legends.xic); }
function drawChrom(){
  const canvas=$("chromCanvas"),ctx=canvas.getContext("2d"); ctx.clearRect(0,0,canvas.width,canvas.height);
  const full=chromFullDomain(); const series=full.series; updateChromLegend(series);
  const domain=state.chromZoom||full;
  const offset=(domain.ymax-domain.ymin)*Number($("offset").value||0); const p=plot(canvas,{...domain,ymax:domain.ymax+offset*(series.length-1)}); axes(ctx,p,"RT (min)",$("chromMode").value.toUpperCase());
  const selected=selectedPeak();
  $("showLocalPeakAligned").disabled=!$("showAligned").checked||!selected||!localAlignmentGatePassed(selected);
  drawSelectedPeakBand(ctx,p,selected,series);
  series.forEach((s,i)=>{ const yoff=(series.length-1-i)*offset; ctx.strokeStyle=s.color; ctx.lineWidth=1.8; ctx.beginPath(); s.x.forEach((x,j)=>{ const px=p.toX(x),py=p.toY(s.y[j]+yoff); if(j===0)ctx.moveTo(px,py); else ctx.lineTo(px,py); }); ctx.stroke(); });
  DATA.peak_results.forEach(pk=>{ const x=p.toX(pk.rt_apex); if(x<p.left||x>p.right)return; ctx.strokeStyle=pk.status==="high_consistency"?"rgba(22,163,74,.35)":"rgba(220,38,38,.45)"; ctx.beginPath(); ctx.moveTo(x,p.bottom); ctx.lineTo(x,p.bottom-9); ctx.stroke(); });
  drawSelectedFeatureRt(ctx,p);
  const localText=selected&&!localAlignmentGatePassed(selected)?"local peak RT correction SKIPPED because no clear two-sided peak shape was found; global RT is retained.":(localPeakAlignmentEnabled()?"local peak RT correction ON: curves are shifted by the selected peak local RT shift.":"global RT alignment is shown.");
  const featureText=Number.isFinite(selectedFeatureRt())?` Purple dashed line marks selected Feature RT ${nice(selectedFeatureRt(),4)} min.`:"";
  canvas._plot=p; $("chromInfo").textContent=`${DATA.peak_results.length} compared TIC peaks; ${localText} The single blue band is the selected TIC integration window.${featureText} Small ticks mark other compared peak apices.`;
  setPan("chromPan",state.chromZoom,full,(xmin,width)=>{state.chromZoom={...state.chromZoom,xmin,xmax:xmin+width};drawAll();});
}
function nearestPointIndex(values,target){
  if(!values?.length)return -1;
  let lo=0,hi=values.length-1;
  while(lo<hi){ const mid=Math.floor((lo+hi)/2); if(Number(values[mid])<target)lo=mid+1; else hi=mid; }
  if(lo<=0)return 0;
  return Math.abs(Number(values[lo])-target)<Math.abs(Number(values[lo-1])-target)?lo:lo-1;
}
function chromHitAt(canvas,clientX,clientY){
  if(!canvas._plot||!DATA)return null;
  const rect=canvas.getBoundingClientRect(), x=(clientX-rect.left)*canvas.width/rect.width, y=(clientY-rect.top)*canvas.height/rect.height, p=canvas._plot;
  if(x<p.left||x>p.right||y<p.top||y>p.bottom)return null;
  const series=chromSeries(), offset=(p.ymax-p.ymin)*Number($("offset").value||0), targetRt=p.fromX(x); let best=null;
  series.forEach((item,index)=>{
    const pointIndex=nearestPointIndex(item.x,targetRt); if(pointIndex<0)return;
    const yOffset=(series.length-1-index)*offset, signal=Number(item.y[pointIndex]||0), pointY=p.toY(signal+yOffset), distance=Math.abs(pointY-y);
    if(distance<=18&&(!best||distance<best.distance))best={sample:item.sample,rt:Number(item.x[pointIndex]),signal,displaySignal:signal+yOffset,color:item.color,distance,mode:$("chromMode").value.toUpperCase(),aligned:$("showAligned").checked};
  });
  return best;
}
function showChromTooltip(event,hit){
  const tip=$("detailTooltip"); if(!tip||!hit){hideDetailTooltip();return;}
  tip.innerHTML=`<div class="tooltip-title"><i class="chrom-legend-swatch" style="background:${hit.color}"></i>${escapeHtml(hit.sample)}</div><div>RT <b>${nice(hit.rt,4)}</b> min</div><div>${escapeHtml(hit.mode)} signal <b>${axisValue(hit.signal)}</b></div>${hit.aligned?"<div class=tooltip-muted>aligned RT</div>":""}`;
  positionTooltip(tip,event); tip.style.display="block";
}
function renderPeakTable(){
  const rows=[...DATA.peak_results].sort((a,b)=>Number(a.rt_apex)-Number(b.rt_apex)).map(pk=>{const local=pk.local_alignment||{},reason=local.local_alignment_rejection_reason||"",localStatus=localAlignmentStatus(pk);return `<tr data-peak="${pk.tic_peak_id}" class="${pk.tic_peak_id===state.selectedPeakId?'selected':''}"><td>${pk.tic_peak_id}</td><td class="status-${pk.status}">${pk.status}</td><td class="local-align-${localStatus}" title="${escapeHtml(reason)}">${localStatus}</td><td>${nice(pk.rt_start)}</td><td>${nice(pk.rt_apex)}</td><td>${nice(pk.rt_end)}</td><td>${nice(pk.peak_consistency_score,3)}</td><td>${nice(pk.spectrum_score,3)}</td><td>${nice(pk.chromatogram_score,3)}</td><td>${nice(pk.presence_absence_score,3)}</td><td>${nice(pk.area,1)}</td><td>${nice(pk.snr,1)}</td></tr>`;});
  $("peakTable").innerHTML=`<tr><th>peak_id</th><th>status</th><th>local RT</th><th>rt_start</th><th>rt_apex</th><th>rt_end</th><th>consistency</th><th>spectrum</th><th>chrom</th><th>presence</th><th>area</th><th>S/N</th></tr>${rows.join("")}`;
  applyTableColumnVisibility("peakTable","peakTable");
  document.querySelectorAll("tr[data-peak]").forEach(row=>row.onclick=()=>{state.selectedPeakId=row.dataset.peak; state.selectedMz=null; state.selectedFeatureRt=null; state.selectedFeatureGroupId=null; state.xic=null; state.detailZoom=null; state.xicZoom=null; drawAll();});
}
function selectedPeak(){ return DATA.peak_results.find(p=>p.tic_peak_id===state.selectedPeakId)||null; }
function detailFullDomain(pk){
  const bins=pk.spectrum_mz_bins||[], matrix=pk.spectrum_raw_matrix||{}; const ymax=Math.max(1,...Object.values(matrix).flat());
  return {xmin:Math.min(...bins,DATA.params.mz_tolerance_da),xmax:Math.max(...bins,2000),ymin:0,ymax};
}
function inferredIonForMz(mz){
  const matched=selectedMs2Evidence(); if(!matched)return null;
  const feature=matched.feature||{}, p=feature.best_psm||feature.candidate_psm||feature.coverage_scan, peaks=p?.spectrum_peaks||[];
  let best=null;
  peaks.forEach(peak=>{
    const label=String(peak.label||""); if(!label)return;
    const observed=Number(peak.mz), delta=Math.abs(observed-Number(mz)); if(!Number.isFinite(observed)||!Number.isFinite(delta))return;
    if(!best||delta<best.delta)best={label,observed,intensity:Number(peak.intensity)||0,delta};
  });
  const tolerance=Math.max(0.03,Math.abs(Number(mz))*25/1e6);
  return best&&best.delta<=tolerance?best:null;
}
function spectrumHitAt(canvas, clientX, clientY){
  const pk=selectedPeak(); if(!pk||!canvas._plot)return null;
  const r=canvas.getBoundingClientRect(); const x=(clientX-r.left)*canvas.width/r.width, y=(clientY-r.top)*canvas.height/r.height;
  const p=canvas._plot; if(x<p.left||x>p.right||y<p.top||y>p.bottom)return null;
  const bins=pk.spectrum_mz_bins||[], matrix=pk.spectrum_raw_matrix||{};
  let best=null;
  bins.forEach((mz,j)=>{
    const px=p.toX(mz); const dx=Math.abs(px-x);
    if(dx>8)return;
    const values=DATA.sample_ids.map(sample=>Number((matrix[sample]||[])[j]||0));
    const maxI=Math.max(...values);
    if(maxI<=0)return;
    const py=p.toY(maxI);
    if(y<py-10||y>p.bottom+8)return;
    if(!best||dx<best.dx){ const ion=inferredIonForMz(Number(mz)); best={mz, index:j, values, maxIntensity:maxI, dx, px, ion}; }
  });
  return best;
}
function spectrumInfoRows(hit){
  if(!hit)return "";
  const rows=DATA.sample_ids.map((sample,i)=>`<tr><td>${sample}</td><td>${nice(hit.values[i],1)}</td></tr>`).join("");
  const positive=hit.values.filter(v=>v>0);
  const max=Math.max(...hit.values,0), min=positive.length?Math.min(...positive):0;
  const fold=min>0?max/min:null;
  const ion=hit.ion?`<span class="ion-badge">推断离子 ${escapeHtml(hit.ion.label)} · Δ ${nice(hit.ion.delta,4)} Da</span>`:"";
  return `<div class="tooltip-title">Selected m/z <b>${nice(hit.mz,5)}</b></div>${ion}<div>max signal <b>${axisValue(max)}</b>${fold?` · fold ${nice(fold,2)}`:""}</div><div class="mini-table"><table><tr><th>sample</th><th>summed intensity</th></tr>${rows}</table></div>`;
}
function selectedSpectrumHit(){
  const pk=selectedPeak(); if(!pk||!state.selectedMz)return null;
  const bins=pk.spectrum_mz_bins||[], matrix=pk.spectrum_raw_matrix||{};
  if(!bins.length)return null;
  let index=0, best=Infinity;
  bins.forEach((mz,j)=>{ const d=Math.abs(Number(mz)-Number(state.selectedMz)); if(d<best){best=d; index=j;} });
  const values=DATA.sample_ids.map(sample=>Number((matrix[sample]||[])[index]||0));
  return {mz:Number(bins[index]), index, values, maxIntensity:Math.max(...values,0)};
}
function selectedFeatureGroup(){
  const pk=selectedPeak(); if(!pk||!state.selectedMz)return null;
  const groups=pk.feature_groups||[];
  let best=null, bestDelta=Infinity;
  groups.forEach(group=>{ const d=Math.abs(Number(group.representative_mz)-Number(state.selectedMz)); if(d<bestDelta){bestDelta=d; best=group;} });
  return best;
}
function showDetailTooltip(e, hit){
  const tip=$("detailTooltip"); if(!tip||!hit){hideDetailTooltip(); return;}
  tip.innerHTML=spectrumInfoRows(hit);
  positionTooltip(tip,e); tip.style.display="block";
}
function featureMapRank(item){
  const rows=featureMapRows();
  const idx=rows.findIndex(row=>row.feature_group_id===item.feature_group_id);
  return idx>=0 ? idx+1 : "";
}
function featureMapTooltipRows(item){
  if(!item)return "";
  const merged=(item.merged_parent_tic_peak_ids||[]).join(", ");
  const pair=pairMetrics(item);
  return `<div><b>${item.feature_group_id||""}</b></div>
    <table>
      <tr><th>rank</th><td>${featureMapRank(item)}</td></tr>
      <tr><th>RT</th><td>${nice(item.representative_rt,4)} min</td></tr>
      <tr><th>true peak m/z</th><td>${nice(truePeakMz(item),5)}</td></tr>
      <tr><th>envelope m/z</th><td>${nice(envelopeRepresentativeMz(item),5)}</td></tr>
      <tr><th>neutral mass</th><td>${nice(item.component_neutral_mass,4)}</td></tr>
      <tr><th>TIC peak</th><td>${item.parent_tic_peak_id||""}</td></tr>
      <tr><th>merged TIC</th><td>${merged||item.parent_tic_peak_id||""}</td></tr>
      <tr><th>reference</th><td>${escapeHtml(sampleShort(pair.reference))}: ${nice(pair.referenceValue,1)}</td></tr>
      <tr><th>test</th><td>${escapeHtml(sampleShort(pair.test))}: ${nice(pair.testValue,1)}</td></tr>
      <tr><th>direction</th><td>${pairDirectionHtml(item)}</td></tr>
      <tr><th>abundance order</th><td>${abundanceOrderHtml(item)}</td></tr>
      <tr><th>pair fold</th><td>${nice(pair.fold,2)}</td></tr>
      <tr><th>ln(test/reference)</th><td>${nice(pair.logRatio,3)}</td></tr>
      <tr><th>pair type</th><td>${pairDifferenceType(item)}</td></tr>
      <tr><th>quantitation</th><td>${item.quantitation_confidence||""}</td></tr>
      <tr><th>ranking</th><td>${nice(item.ranking_score,3)}</td></tr>
      <tr><th>square-size abundance</th><td>${nice(featureMapPairAbundance(item),1)}</td></tr>
      <tr><th>max raw area</th><td>${nice(item.max_area,1)}</td></tr>
    </table>`;
}
function showFeatureMapTooltip(e, item){
  const tip=$("detailTooltip"); if(!tip||!item){hideDetailTooltip(); return;}
  tip.innerHTML=featureMapTooltipRows(item);
  positionTooltip(tip,e); tip.style.display="block";
}
async function selectSpectrumMz(mz, representativeRt=null, featureGroupId=null, allowFeatureAutoMatch=true, preserveModificationContext=false){
  const requestId=++state.xicRequestId;
  if(!preserveModificationContext){state.selectedModificationFormKey=null;state.selectedModificationFormContext=null;}
  state.selectedMz=Number(mz);
  state.selectedFeatureGroupId=featureGroupId?String(featureGroupId):null;
  if(!state.selectedFeatureGroupId&&allowFeatureAutoMatch){
    const matchedGroup=selectedFeatureGroup();
    state.selectedFeatureGroupId=matchedGroup?.feature_group_id?String(matchedGroup.feature_group_id):null;
  }
  const suppliedRt=representativeRt===null||representativeRt===undefined||representativeRt===""?Number.NaN:Number(representativeRt);
  state.selectedFeatureRt=Number.isFinite(suppliedRt)?suppliedRt:null;
  if(!Number.isFinite(state.selectedFeatureRt)){
    const groupRt=Number(selectedFeatureGroup()?.representative_rt);
    state.selectedFeatureRt=Number.isFinite(groupRt)?groupRt:null;
  }
  if(state.selectedFeatureGroupId)ensureMsmsData();
  state.xicZoom=null; state.xic=null; state.xicLoading=true; drawAll();
  try {
    await loadXic(requestId);
    if(requestId===state.xicRequestId){ state.xicLoading=false; drawAll(); }
  } catch(err) {
    if(requestId===state.xicRequestId){ state.xicLoading=false; $("xicInfo").textContent=`XIC request failed: ${err.message||err}`; }
  }
}
function drawDetail(){
  const canvas=$("detailCanvas"),ctx=canvas.getContext("2d"); ctx.clearRect(0,0,canvas.width,canvas.height);
  updateSpectrumLegend();
  const pk=selectedPeak(); if(!pk){ $("detailInfo").textContent="Select one TIC peak from the table or TIC plot to compare its summed spectrum."; canvas._plot=null; $("detailPan").disabled=true; return; }
  state.selectedPeakId=pk.tic_peak_id;
  const bins=pk.spectrum_mz_bins||[], matrix=pk.spectrum_raw_matrix||{}, full=detailFullDomain(pk);
  const p=plot(canvas,state.detailZoom||full); axes(ctx,p,"m/z","summed spectrum");
  DATA.sample_ids.forEach((sample,i)=>{ const v=matrix[sample]||[]; ctx.strokeStyle=color(i); v.forEach((inten,j)=>{ if(inten<=0)return; const x=p.toX(bins[j]); if(x<p.left||x>p.right)return; ctx.beginPath(); ctx.moveTo(x,p.bottom); ctx.lineTo(x,p.toY(inten)); ctx.stroke(); }); });
  if(state.selectedMz){ const x=p.toX(state.selectedMz); if(x>=p.left&&x<=p.right){ ctx.save(); ctx.strokeStyle="#111827"; ctx.setLineDash([6,4]); ctx.beginPath(); ctx.moveTo(x,p.top); ctx.lineTo(x,p.bottom); ctx.stroke(); ctx.fillStyle="#111827"; ctx.fillText(`m/z ${nice(state.selectedMz,4)}`,Math.min(x+6,p.right-95),p.top+14); ctx.restore(); }}
  canvas._plot=p;
  setPan("detailPan",state.detailZoom,full,(xmin,width)=>{state.detailZoom={...state.detailZoom,xmin,xmax:xmin+width};drawAll();});
  const scoreText=`shape ${nice(pk.shape_score)}; area ${nice(pk.area_score)}; alignment ${nice(pk.apex_score)}; width ${nice(pk.width_score)}; cosine ${nice(pk.cosine_score)}; top m/z ${nice(pk.top_mz_overlap_score)}; rel abundance ${nice(pk.relative_abundance_score)}; presence ${nice(pk.presence_absence_score)}; apex spectrum ${nice(pk.apex_spectrum_score)}`;
  const local=pk.local_alignment||{}, shifts=local.local_peak_shift_by_sample||{}, profile=local.profile_score_by_sample||{};
  const localText=DATA.sample_ids.map(sample=>`${sample}: shift ${nice(shifts[sample]||0,4)} min, profile ${nice(profile[sample]||0,3)}`).join("; ");
  const gateText=local.local_alignment_quality_gate_passed===false?`local alignment skipped (${local.local_alignment_rejection_reason||"no clear peak shape"}); global RT retained`:`local ${local.alignment_method||"apex"}`;
  $("detailInfo").textContent=`${pk.tic_peak_id} ${pk.status}; RT ${nice(pk.rt_start)}-${nice(pk.rt_end)}; consistency ${nice(pk.peak_consistency_score)}; ${scoreText}; ${gateText}; ${localText}`;
}
function renderMzTable(){
  const pk=selectedPeak(); const changes=sortedTableRows(pk?.top_changed_mz||[],state.mzSort);
  const hit=selectedSpectrumHit();
  const group=selectedFeatureGroup();
  const groupRows=group?DATA.sample_ids.map(sample=>{ const f=(group.features_by_sample||{})[sample]||{}; return `<tr><td>${sample}</td><td>${nice((group.area_by_sample||{})[sample],1)}</td><td>${nice((group.normalized_area_by_sample||{})[sample],1)}</td><td>${(group.presence_by_sample||{})[sample]?"present":"absent/weak"}</td><td>${f.quantitation_status||""}</td><td>${nice(f.signal_to_noise,1)}</td><td>${nice((group.match_score_by_sample||{})[sample],3)}</td><td>${nice((group.rt_correction_by_sample||{})[sample],4)}</td></tr>`; }).join(""):"";
  const rows=changes.map(item=>`<tr data-mz="${item.mz}" data-rt="${item.representative_rt??""}" data-feature-id="${escapeHtml(item.feature_group_id||"")}" class="${Number(state.selectedMz)===Number(item.mz)?'selected':''}"><td>${nice(item.mz,4)}</td><td>${pairDifferenceType(item)}</td><td>${pairDirectionHtml(item)}</td><td>${abundanceOrderHtml(item)}</td><td>${item.quantitation_confidence||""}</td><td>${nice(item.ranking_score,3)}</td><td>${nice(item.cv,3)}</td><td>${nice(pairMetrics(item).fold,2)}</td></tr>`);
  $("mzTable").innerHTML=`<tr>${sortableTh("true peak m/z","mz","mz")}${sortableTh("pair type","pair_type","mz")}${sortableTh("test direction","pair_log_ratio","mz")}${sortableTh("abundance order","pair_log_ratio","mz")}${sortableTh("confidence","quantitation_confidence","mz")}${sortableTh("ranking","ranking_score","mz")}${sortableTh("CV","cv","mz")}${sortableTh("pair fold","pair_fold","mz")}</tr>${rows.join("")}`;
  applyTableColumnVisibility("mzTable","mzTable");
  bindSortableHeaders("mz",renderMzTable);
  document.querySelectorAll("tr[data-mz]").forEach(row=>row.onclick=async()=>{await selectSpectrumMz(Number(row.dataset.mz),row.dataset.rt,row.dataset.featureId);});
}
async function loadXic(requestId){
  const pk=selectedPeak(); if(!pk||!state.selectedMz)return;
  const full=state.xicFull ? "&full=1" : "";
  const payload=await fetchJson(apiUrl(`/api/xic?peak_id=${encodeURIComponent(pk.tic_peak_id)}&mz=${encodeURIComponent(state.selectedMz)}${full}`));
  if(requestId&&requestId!==state.xicRequestId)return;
  state.xic=payload;
  state.xicZoom=state.xicFull?null:{xmin:Number(payload.rt_start),xmax:Number(payload.rt_end),ymin:0,ymax:null};
}
function xicFullDomain(){
  if(!state.xic)return null;
  const series=state.xic.xic_by_sample||{}; const xs=Object.values(series).flatMap(s=>s.rt), ys=Object.values(series).flatMap(s=>s.intensity);
  return {xmin:Math.min(...xs),xmax:Math.max(...xs),ymin:0,ymax:Math.max(...ys,1),xs,ys};
}
function drawXic(){
  const canvas=$("xicCanvas"),ctx=canvas.getContext("2d"); ctx.clearRect(0,0,canvas.width,canvas.height);
  updateXicLegend();
  if(state.xicLoading){ $("xicInfo").textContent=`Loading XIC for m/z ${nice(state.selectedMz,5)}...`; canvas._plot=null; $("xicPan").disabled=true; return; }
  if(!state.xic){ $("xicInfo").textContent="Select a Top changed m/z to extract XIC."; canvas._plot=null; $("xicPan").disabled=true; return; }
  const series=state.xic.xic_by_sample||{}, full=xicFullDomain();
  const z=state.xicZoom||{xmin:Number(state.xic.rt_start),xmax:Number(state.xic.rt_end),ymin:0,ymax:null};
  const visibleYs=Object.values(series).flatMap(s=>s.rt.map((rt,j)=>rt>=z.xmin&&rt<=z.xmax?s.intensity[j]:null).filter(v=>v!==null));
  const domain={xmin:z.xmin,xmax:z.xmax,ymin:0,ymax:z.ymax||Math.max(1,...visibleYs)};
  const p=plot(canvas,domain); axes(ctx,p,"aligned RT (min)","XIC intensity");
  const x0=p.toX(state.xic.integration_rt_start), x1=p.toX(state.xic.integration_rt_end);
  if(reportSettings.xicIntegrationBand&&x1>=p.left&&x0<=p.right){ ctx.save(); ctx.fillStyle="rgba(34,197,94,.12)"; ctx.fillRect(Math.max(p.left,x0),p.top,Math.min(p.right,x1)-Math.max(p.left,x0),p.bottom-p.top); ctx.restore(); }
  DATA.sample_ids.forEach((sample,i)=>{ const s=series[sample]; if(!s)return; ctx.strokeStyle=color(i); ctx.beginPath(); let started=false; s.rt.forEach((x,j)=>{ if(x<domain.xmin||x>domain.xmax){started=false; return;} const px=p.toX(x),py=p.toY(s.intensity[j]); if(!started){ctx.moveTo(px,py); started=true;} else ctx.lineTo(px,py); }); ctx.stroke(); });
  canvas._plot=p;
  setPan("xicPan",state.xicZoom,full,(xmin,width)=>{state.xicZoom={...state.xicZoom,xmin,xmax:xmin+width,ymax:null};drawAll();});
  $("xicInfo").textContent=`Single-centroid XIC m/z ${nice(state.xic.target_mz,5)} +/- ${nice(state.xic.mz_tolerance,5)} Da; parent peak ${state.xic.peak_id}; ${state.xic.full_run?"full XIC":"local XIC"}; green band is feature/TIC integration RT ${nice(state.xic.integration_rt_start)}-${nice(state.xic.integration_rt_end)}`;
}
function xicHitAt(canvas,clientX,clientY){
  if(!canvas._plot||!state.xic)return null;
  const rect=canvas.getBoundingClientRect(), x=(clientX-rect.left)*canvas.width/rect.width, y=(clientY-rect.top)*canvas.height/rect.height, p=canvas._plot;
  if(x<p.left||x>p.right||y<p.top||y>p.bottom)return null;
  const targetRt=p.fromX(x), series=state.xic.xic_by_sample||{}; let best=null;
  DATA.sample_ids.forEach((sample,index)=>{
    const item=series[sample]; if(!item?.rt?.length)return;
    const pointIndex=nearestPointIndex(item.rt,targetRt), rt=Number(item.rt[pointIndex]), intensity=Number(item.intensity?.[pointIndex]||0), pointY=p.toY(intensity), distance=Math.abs(pointY-y);
    if(distance<=18&&(!best||distance<best.distance))best={sample,rt,intensity,area:Number((state.xic.integration_by_sample?.[sample]||{}).area||0),color:color(index),distance};
  });
  return best;
}
function showXicTooltip(event,hit){
  const tip=$("detailTooltip"); if(!tip||!hit){hideDetailTooltip();return;}
  tip.innerHTML=`<div class="tooltip-title"><i class="chrom-legend-swatch" style="background:${hit.color}"></i>${escapeHtml(hit.sample)}</div><div>RT <b>${nice(hit.rt,4)}</b> min</div><div>XIC signal <b>${axisValue(hit.intensity)}</b></div><div>integrated area <b>${axisValue(hit.area)}</b></div><div class="tooltip-muted">m/z ${nice(state.xic.target_mz,5)} · ${state.xic.full_run?"overall":"local"} XIC</div>`;
  positionTooltip(tip,event); tip.style.display="block";
}
function featureMapRows(){ return (DATA.global_feature_groups||[]).filter(item=>Number(item.representative_rt)>0&&Number(item.representative_mz)>0&&Number(item.max_area)>0); }
function featureMapFullDomain(){
  const rows=featureMapRows(); if(!rows.length)return {xmin:0,xmax:1,ymin:0,ymax:1};
  const xs=rows.map(item=>Number(item.representative_rt)), ys=rows.map(item=>Number(item.representative_mz));
  const xpad=Math.max((Math.max(...xs)-Math.min(...xs))*0.04,0.05), ypad=Math.max((Math.max(...ys)-Math.min(...ys))*0.04,1);
  return {xmin:Math.min(...xs)-xpad,xmax:Math.max(...xs)+xpad,ymin:Math.max(0,Math.min(...ys)-ypad),ymax:Math.max(...ys)+ypad};
}
function lerp(a,b,t){ return a+(b-a)*t; }
function heatmapMaxFold(rows){
  const finite=(rows||[]).map(item=>Math.min(MAX_FEATURE_FOLD,Number(item.max_fold_change))).filter(value=>Number.isFinite(value)&&value>=1);
  return Math.min(MAX_FEATURE_FOLD,Math.max(4,...finite));
}
function heatmapFoldNorm(item,maxFold){
  const fold=Math.min(MAX_FEATURE_FOLD,Number(item?.max_fold_change));
  if(!Number.isFinite(fold)||fold<1)return item?.difference_type==="presence_absence"?1:0;
  return Math.max(0,Math.min(1,Math.log(fold)/Math.max(Math.log(Math.max(4,maxFold)),1e-9)));
}
function featureFoldColor(item, rankNorm=0, scaleOnly=false, maxFold=4){
  const fold=Number(item?.max_fold_change);
  let strength=heatmapFoldNorm(item,maxFold);
  if(item?.difference_type==="presence_absence"&&!Number.isFinite(fold))strength=1;
  const palette=[[22,163,74],[234,179,8],[249,115,22],[220,38,38],[88,28,135]], slate=[100,116,139];
  const palettePosition=strength*(palette.length-1), leftIndex=Math.min(palette.length-2,Math.floor(palettePosition));
  const t=palettePosition-leftIndex, left=palette[leftIndex], right=palette[leftIndex+1];
  let r=Math.round(lerp(left[0],right[0],t)), g=Math.round(lerp(left[1],right[1],t)), b=Math.round(lerp(left[2],right[2],t));
  if(!scaleOnly){
    const confidence=item?.quantitation_confidence||"high";
    const mute=confidence==="partial_detection"?0.55:(confidence==="low"?0.85:0);
    r=Math.round(lerp(r,slate[0],mute)); g=Math.round(lerp(g,slate[1],mute)); b=Math.round(lerp(b,slate[2],mute));
  }
  const alpha=Math.max(0.28,Math.min(0.95,0.35+0.40*strength+0.20*Math.sqrt(Math.max(0,Math.min(1,rankNorm)))));
  return `rgba(${r},${g},${b},${alpha})`;
}
function pairMaxAbsLn(rows){ return Math.min(MAX_FEATURE_LN_FOLD,Math.max(Math.log(2),...(rows||[]).map(item=>Math.abs(pairMetrics(item).logRatio)).filter(Number.isFinite))); }
function featureMapPairAbundance(item){
  const pair=selectedPair(), values=featureValueMap(item)||{};
  const hasReference=Object.prototype.hasOwnProperty.call(values,pair.reference);
  const hasTest=Object.prototype.hasOwnProperty.call(values,pair.test);
  if(hasReference||hasTest){
    const reference=Number(values[pair.reference]), test=Number(values[pair.test]);
    return Math.max(0,Number.isFinite(reference)?reference:0,Number.isFinite(test)?test:0);
  }
  return Math.max(0,Number(item?.max_area)||0);
}
function featureMapAbundanceScale(rows){
  const values=(rows||[]).map(item=>Math.log1p(featureMapPairAbundance(item))).filter(Number.isFinite).sort((a,b)=>a-b);
  if(!values.length)return {low:0,high:1};
  const percentile=fraction=>{
    const position=(values.length-1)*fraction, lower=Math.floor(position), upper=Math.ceil(position);
    return lower===upper?values[lower]:lerp(values[lower],values[upper],position-lower);
  };
  const low=percentile(0.05), high=percentile(0.95);
  return {low,high:Math.max(high,low+1e-9)};
}
function featureMapAbundanceNorm(item,scale){
  const value=Math.log1p(featureMapPairAbundance(item));
  return Math.max(0,Math.min(1,(value-scale.low)/Math.max(scale.high-scale.low,1e-9)));
}
function featureMapHalfSide(item,scale){
  const minimum=3, maximum=15, normalized=featureMapAbundanceNorm(item,scale);
  return Math.sqrt(minimum*minimum+(maximum*maximum-minimum*minimum)*normalized);
}
function logFoldTicks(maxFold,includeOne=true){
  const ticks=[1,10,100,1000].filter(value=>value<=maxFold+1e-9);
  if(maxFold>1&&maxFold<10&&!ticks.some(value=>Math.abs(value-maxFold)<1e-6))ticks.push(maxFold);
  return ticks.filter(value=>includeOne||value>1).sort((a,b)=>a-b);
}
function pairDirectionColor(item,rankNorm=0,scaleOnly=false,maxAbsLn=Math.log(2)){
  const signed=Math.max(-1,Math.min(1,pairMetrics(item).logRatio/Math.max(maxAbsLn,1e-9)));
  const palette=[[30,64,175],[96,165,250],[226,232,240],[248,113,113],[153,27,27]], position=(signed+1)*0.5*(palette.length-1);
  const leftIndex=Math.min(palette.length-2,Math.floor(position)), t=position-leftIndex, left=palette[leftIndex], right=palette[leftIndex+1], slate=[100,116,139];
  let r=Math.round(lerp(left[0],right[0],t)),g=Math.round(lerp(left[1],right[1],t)),b=Math.round(lerp(left[2],right[2],t));
  if(!scaleOnly){
    const confidence=item?.quantitation_confidence||"high", mute=confidence==="partial_detection"?0.50:(confidence==="low"?0.82:0);
    r=Math.round(lerp(r,slate[0],mute));g=Math.round(lerp(g,slate[1],mute));b=Math.round(lerp(b,slate[2],mute));
  }
  const strength=Math.abs(signed),alpha=Math.max(0.28,Math.min(0.96,0.34+0.43*strength+0.18*Math.sqrt(Math.max(0,Math.min(1,rankNorm)))));
  return `rgba(${r},${g},${b},${alpha})`;
}
function featureMapPriority(item, maxRanking=1, maxScale=4){
  const quality=Math.max(0,Math.min(1,Number(item.ranking_confidence_weight??1)));
  const foldNorm=(state.featureMapMode==="direction"?Math.abs(pairMetrics(item).logRatio)/Math.max(maxScale,1e-9):heatmapFoldNorm(item,maxScale))*quality;
  const rankNorm=Math.max(0,Math.min(1,(Number(item.ranking_score)||0)/Math.max(maxRanking,1e-9)));
  return foldNorm*0.78 + Math.sqrt(rankNorm)*0.22;
}
function updateFeatureMapLegend(){
  const container=$("featureMapLegend"); if(!container)return;
  const rows=featureMapRows(), directionMode=state.featureMapMode==="direction";
  const maxFold=directionMode?Math.exp(pairMaxAbsLn(rows)):heatmapMaxFold(rows);
  const pair=selectedPair();
  const labels=directionMode
    ? [`${escapeHtml(sampleShort(pair.test))} lower`,"equal",`${escapeHtml(sampleShort(pair.test))} higher`]
    : ["1x","10x",`${nice(maxFold,maxFold>=10?0:1)}x`];
  container.innerHTML=`<span class="legend-gradient ${directionMode?"direction":"magnitude"}"></span><span class="legend-label">${labels[0]}</span><span class="legend-label">${labels[1]}</span><span class="legend-label">${labels[2]}</span><span class="legend-caption">${directionMode?"red = test higher · blue = test lower":"color = ln(fold)"}</span>`;
}
function featureMapHitAt(canvas, clientX, clientY){
  if(!canvas._plot)return null;
  const r=canvas.getBoundingClientRect(); const x=(clientX-r.left)*canvas.width/r.width, y=(clientY-r.top)*canvas.height/r.height;
  const p=canvas._plot; if(x<p.left||x>p.right||y<p.top||y>p.bottom)return null;
  const rows=featureMapRows();
  const maxRanking=Math.max(...rows.map(item=>Number(item.ranking_score)||0),1e-9);
  const maxScale=state.featureMapMode==="direction"?pairMaxAbsLn(rows):heatmapMaxFold(rows);
  const abundanceScale=featureMapAbundanceScale(rows);
  let best=null;
  rows.forEach(item=>{
    const px=p.toX(Number(item.representative_rt)), py=p.toY(Number(item.representative_mz));
    const d=Math.hypot(px-x,py-y);
    const halfSide=featureMapHalfSide(item,abundanceScale);
    const priority=featureMapPriority(item,maxRanking,maxScale);
    const inside=Math.abs(px-x)<=halfSide+3&&Math.abs(py-y)<=halfSide+3;
    if(inside&&(!best||priority>best.priority+0.03||(Math.abs(priority-best.priority)<=0.03&&d<best.d)))best={item,d,priority};
  });
  return best?.item||null;
}
async function selectGlobalFeature(item,preserveModificationContext=false){
  if(!item)return;
  if(!preserveModificationContext){state.selectedModificationFormKey=null;state.selectedModificationFormContext=null;}
  state.selectedPeakId=item.parent_tic_peak_id;
  state.selectedMz=null;
  state.selectedFeatureRt=null;
  state.selectedFeatureGroupId=item.feature_group_id?String(item.feature_group_id):null;
  state.scrollGlobalSelectionIntoView=true;
  state.detailZoom=null;
  state.xic=null;
  state.xicZoom=null;
  await selectSpectrumMz(Number(item.representative_mz),item.representative_rt,state.selectedFeatureGroupId,true,preserveModificationContext);
}
function drawFeatureMap(){
  const canvas=$("featureMapCanvas"),ctx=canvas.getContext("2d"); ctx.clearRect(0,0,canvas.width,canvas.height);
  const rows=featureMapRows(), full=featureMapFullDomain();
  if(!rows.length){ $("featureMapInfo").textContent="No global feature groups."; canvas._plot=null; $("featureMapPan").disabled=true; updateFeatureMapLegend(); return; }
  const domain=state.featureMapZoom||full; const p=plot(canvas,domain); axes(ctx,p,"RT (min)","m/z");
  const maxRanking=Math.max(...rows.map(item=>Number(item.ranking_score)||0),1e-9);
  const directionMode=state.featureMapMode==="direction", maxFold=heatmapMaxFold(rows), maxAbsLn=pairMaxAbsLn(rows), maxScale=directionMode?maxAbsLn:maxFold;
  const abundanceScale=featureMapAbundanceScale(rows);
  const drawRows=[...rows].sort((a,b)=>featureMapPriority(a,maxRanking,maxScale)-featureMapPriority(b,maxRanking,maxScale));
  drawRows.forEach(item=>{
    const rt=Number(item.representative_rt), mz=Number(item.representative_mz);
    if(rt<domain.xmin||rt>domain.xmax||mz<domain.ymin||mz>domain.ymax)return;
    const x=p.toX(rt), y=p.toY(mz);
    const radius=featureMapHalfSide(item,abundanceScale);
    ctx.fillStyle=directionMode?pairDirectionColor(item,0,false,maxAbsLn):featureFoldColor(item,0,false,maxFold);
    ctx.fillRect(x-radius,y-radius,radius*2,radius*2);
    const selectedByFeatureId=state.selectedFeatureGroupId&&featureGroupIds(item).includes(String(state.selectedFeatureGroupId));
    if(selectedByFeatureId||(item.parent_tic_peak_id===state.selectedPeakId&&Math.abs(Number(item.representative_mz)-Number(state.selectedMz||0))<1e-6)){
      ctx.strokeStyle="#111827"; ctx.lineWidth=2; ctx.strokeRect(x-radius-2,y-radius-2,radius*2+4,radius*2+4); ctx.lineWidth=1;
    }
  });
  canvas._plot=p;
  updateFeatureMapLegend();
  setPan("featureMapPan",state.featureMapZoom,full,(xmin,width)=>{state.featureMapZoom={...state.featureMapZoom,xmin,xmax:xmin+width};drawAll();});
  const pair=selectedPair();
  $("featureMapInfo").textContent=directionMode?`${rows.length} features; ${sampleShort(pair.test)} versus ${sampleShort(pair.reference)}; red = test higher, blue = test lower; color uses capped |ln(test/reference)| (maximum ${MAX_FEATURE_FOLD}x); square area uses the higher normalized abundance in the selected pair on a robust ln scale.`:`${rows.length} features; color uses ln(fold), scaled from 1x to ${nice(maxFold,1)}x and capped at ${MAX_FEATURE_FOLD}x; square area uses the higher normalized abundance in the selected pair on a robust ln scale.`;
}
function renderGlobalFeatureTable(){
  const componentRows=currentComponentRows();
  const source=componentRows.length?componentRows:(DATA.global_feature_groups||[]).map(item=>({...item,component_group_id:`MS1SINGLE_${item.feature_group_id}`,component_label:item.feature_group_id,component_member_count:1,component_primary_feature_id:item.feature_group_id,members:[item]}));
  const sorted=sortedTableRows(source,state.globalSort), rows=[];
  sorted.forEach((item,i)=>{
    const itemTrueMz=truePeakMz(item), itemEnvelopeMz=envelopeRepresentativeMz(item);
    const pair=pairMetrics(item), componentId=String(item.component_group_id||item.feature_group_id||""), members=item.members||[], expandable=members.length>1, expanded=state.expandedComponents.has(componentId);
    const componentSelected=componentContainsFeature(item,state.selectedFeatureGroupId);
    const legacySelected=item.parent_tic_peak_id===state.selectedPeakId && Math.abs(Number(itemTrueMz)-Number(state.selectedMz||0))<1e-6;
    const parentSelected=state.selectedFeatureGroupId?componentSelected&&!expanded:legacySelected;
    const parentFeatureId=String(item.component_primary_feature_id||item.feature_group_id||"");
    const toggle=expandable?`<button class="component-toggle" data-component-toggle="${escapeHtml(componentId)}" title="${expanded?'Collapse':'Expand'} ${members.length} MS1 Features">${expanded?'&#9662;':'&#9656;'}</button>`:`<span class="component-toggle-placeholder"></span>`;
    const candidateEvidence=componentCandidateEvidence(item), candidateBadge=!item.identified_component&&candidateEvidence?`<span class="candidate-badge">MS2 D candidate</span>`:"";
    const identification=(item.identified_component?`<span class="identification-badge">MS2 ${escapeHtml(String(item.component_confidence||"").replace(/_.*/,""))}</span>`:(item.inferred_component?`<span class="inference-badge">${escapeHtml(String(item.component_confidence||"MS1 inferred").replace("_"," "))}</span>`:""))+candidateBadge;
    const chargeText=(item.component_charge_states||[]).length?`z=${item.component_charge_states.join('/')}`:"";
    const isotopeFit=Number(item.component_isotope_fit_score), isotopeFitText=Number.isFinite(isotopeFit)?`; isotope fit ${nice(isotopeFit,2)}${item.component_monoisotope_corrected?"; mono corrected":""}`:"";
    const componentMeta=(item.identified_component||item.inferred_component)?`${members.length} ions${chargeText?`; ${chargeText}`:""}${isotopeFitText}`:`unresolved singleton`;
    rows.push(`<tr data-global-feature="1" data-component-parent="1" data-component-id="${escapeHtml(componentId)}" data-feature-id="${escapeHtml(parentFeatureId)}" data-selection-target="${parentSelected?'1':'0'}" data-peak="${escapeHtml(item.parent_tic_peak_id)}" data-rt="${item.representative_rt}" data-mz="${itemTrueMz}" class="component-row ${parentSelected?'selected':''}"><td>${i+1}</td><td class="component-cell">${toggle}<span class="component-label">${escapeHtml(componentDisplayLabel(item))}</span>${identification}<span class="component-meta">${escapeHtml(componentMeta)}</span></td><td>${escapeHtml(item.parent_tic_peak_id)}</td><td>${nice(item.representative_rt,3)}</td><td>${nice(itemTrueMz,5)}</td><td>${nice(itemEnvelopeMz,5)}</td><td>${nice(item.component_neutral_mass,4)}</td><td>${pairDifferenceType(item)}</td><td>${pairDirectionHtml(item)}</td><td>${abundanceOrderHtml(item)}</td><td>${nice(pair.fold,2)}</td><td>${nice(pair.logRatio,3)}</td><td>${nice(item.ranking_score,3)}</td><td>${nice(item.max_area,1)}</td><td>${escapeHtml(item.quantitation_confidence||"")}</td><td>${members.length||1}</td><td>${item.merged_feature_count||1}</td><td>${escapeHtml(item.difference_type||"")}</td></tr>`);
    if(expandable&&expanded)members.forEach(member=>{
      const memberTrueMz=truePeakMz(member), memberEnvelopeMz=envelopeRepresentativeMz(member), childFeatureSelected=state.selectedFeatureGroupId&&featureGroupIds(member).includes(String(state.selectedFeatureGroupId)), childSelected=state.selectedFeatureGroupId?childFeatureSelected:(member.parent_tic_peak_id===state.selectedPeakId&&Math.abs(Number(memberTrueMz)-Number(state.selectedMz||0))<1e-6), childPair=pairMetrics(member), offset=Number(member.component_isotope_offset||0), isotope=offset===0?"M":(offset>0?`M+${offset}`:`M${offset}`), charge=member.component_charge?`z=${member.component_charge}`:"z=?", relation=String(member.component_link_type||"").replace("_envelope","").replaceAll("_"," ");
      const memberFit=Number(member.component_isotope_fit_score), fitText=Number.isFinite(memberFit)?`; isotope fit ${nice(memberFit,2)}${Number(member.component_monoisotopic_offset||0)>0?`; M+${member.component_monoisotopic_offset} corrected`:""}`:"";
      const resolved=member.component_observed_first_isotope_mz?`; first isotope ${nice(member.component_observed_first_isotope_mz,4)}; ${member.component_isotope_peak_count||0} peaks${fitText}`:"";
      rows.push(`<tr data-global-feature="1" data-feature-id="${escapeHtml(member.feature_group_id||"")}" data-selection-target="${childSelected?'1':'0'}" data-peak="${escapeHtml(member.parent_tic_peak_id)}" data-rt="${member.representative_rt}" data-mz="${memberTrueMz}" class="component-child ${childSelected?'selected':''}"><td></td><td class="component-cell"><span class="component-tree">&#8627;</span>${escapeHtml(member.feature_group_id)}<span class="component-meta">${escapeHtml(`${charge}; ${isotope}; ${relation}${resolved}`)}</span></td><td>${escapeHtml(member.parent_tic_peak_id)}</td><td>${nice(member.representative_rt,3)}</td><td>${nice(memberTrueMz,5)}</td><td>${nice(memberEnvelopeMz,5)}</td><td>${nice(member.component_neutral_mass,4)}</td><td>${pairDifferenceType(member)}</td><td>${pairDirectionHtml(member)}</td><td>${abundanceOrderHtml(member)}</td><td>${nice(childPair.fold,2)}</td><td>${nice(childPair.logRatio,3)}</td><td>${nice(member.ranking_score,3)}</td><td>${nice(member.max_area,1)}</td><td>${escapeHtml(member.quantitation_confidence||"")}</td><td>&mdash;</td><td>${member.merged_feature_count||1}</td><td>${escapeHtml(member.difference_type||"")}</td></tr>`);
    });
  });
  $("globalFeatureTable").innerHTML=`<tr>${sortableTh("rank","ranking_score","global")}${sortableTh("component / Feature","component_label","global")}${sortableTh("TIC peak","parent_tic_peak_id","global")}${sortableTh("RT","representative_rt","global")}${sortableTh("true peak m/z","true_peak_mz","global")}${sortableTh("envelope m/z","envelope_representative_mz","global")}${sortableTh("neutral mass","component_neutral_mass","global")}${sortableTh("pair type","pair_type","global")}${sortableTh("test direction","pair_log_ratio","global")}${sortableTh("abundance order","pair_log_ratio","global")}${sortableTh("pair fold","pair_fold","global")}${sortableTh("ln(test/ref)","pair_log_ratio","global")}${sortableTh("ranking","ranking_score","global")}${sortableTh("max raw area","max_area","global")}${sortableTh("confidence","quantitation_confidence","global")}${sortableTh("ions","component_member_count","global")}<th>merged</th>${sortableTh("cohort type","difference_type","global")}</tr>${rows.join("")}`;
  applyTableColumnVisibility("globalFeatureTable","globalFeatureTable");
  bindSortableHeaders("global",renderGlobalFeatureTable);
  document.querySelectorAll("button[data-component-toggle]").forEach(button=>button.onclick=event=>{event.stopPropagation();const id=button.dataset.componentToggle;if(state.expandedComponents.has(id))state.expandedComponents.delete(id);else state.expandedComponents.add(id);renderGlobalFeatureTable();});
  document.querySelectorAll("tr[data-global-feature]").forEach(row=>row.onclick=async()=>{
    let featureId=String(row.dataset.featureId||"");
    if(row.dataset.componentParent==="1"){
      await ensureMsmsData();
      const componentId=String(row.dataset.componentId||"");
      const component=currentComponentRows().find(item=>String(item.component_group_id||item.feature_group_id||"")===componentId);
      if(component)featureId=preferredComponentEvidenceFeatureId(component)||featureId;
    }
    await selectGlobalFeature({feature_group_id:featureId,parent_tic_peak_id:row.dataset.peak, representative_rt:row.dataset.rt, representative_mz:Number(row.dataset.mz)});
  });
  if(state.scrollGlobalSelectionIntoView){
    const target=$("globalFeatureTable").querySelector('tr[data-selection-target="1"]'), container=$("globalFeatureTable").parentElement;
    if(target&&container){
      const top=target.offsetTop, bottom=top+target.offsetHeight;
      if(top<container.scrollTop)container.scrollTop=Math.max(0,top-4);
      else if(bottom>container.scrollTop+container.clientHeight)container.scrollTop=Math.max(0,bottom-container.clientHeight+4);
    }
    state.scrollGlobalSelectionIntoView=false;
  }
}
function exportGlobalFeatureTable(){
  const componentRows=currentComponentRows();
  const source=componentRows.length?componentRows:(DATA.global_feature_groups||[]);
  const headers=["rank","component_or_feature","tic_peak","rt_min","true_peak_mz","envelope_mz","neutral_mass","isotope_fit","monoisotope_corrected","pair_type","test_direction","abundance_order","pair_fold","ln_test_reference","ranking","max_raw_area","confidence","ions","merged","cohort_type"];
  const lines=[headers];
  source.forEach((item,index)=>{
    const pair=pairMetrics(item), trueMz=truePeakMz(item), envelopeMz=envelopeRepresentativeMz(item), members=item.members||[];
    const order=pair.direction==="higher"?`${sampleShort(pair.test)} > ${sampleShort(pair.reference)}`:pair.direction==="lower"?`${sampleShort(pair.reference)} > ${sampleShort(pair.test)}`:"equal";
    lines.push([index+1,componentDisplayLabel(item),item.parent_tic_peak_id||"",item.representative_rt,trueMz,envelopeMz,item.component_neutral_mass||"",item.component_isotope_fit_score??"",Boolean(item.component_monoisotope_corrected),pairDifferenceType(item),pair.direction,order,pair.fold,pair.logRatio,item.ranking_score,item.max_area,item.quantitation_confidence||"",members.length||1,item.merged_feature_count||1,item.difference_type||""]);
    members.forEach(member=>{
      const childPair=pairMetrics(member);
      lines.push(["",member.feature_group_id||"",member.parent_tic_peak_id||"",member.representative_rt,truePeakMz(member),envelopeRepresentativeMz(member),member.component_neutral_mass||"",member.component_isotope_fit_score??"",Number(member.component_monoisotopic_offset||0)>0,pairDifferenceType(member),childPair.direction,childPair.direction==="higher"?`${sampleShort(childPair.test)} > ${sampleShort(childPair.reference)}`:childPair.direction==="lower"?`${sampleShort(childPair.reference)} > ${sampleShort(childPair.test)}`:"equal",childPair.fold,childPair.logRatio,member.ranking_score,member.max_area,member.quantitation_confidence||"",1,member.merged_feature_count||1,member.difference_type||""]);
    });
  });
  const csv=lines.map(row=>row.map(value=>`"${String(value??"").replaceAll('"','""')}"`).join(",")).join("\r\n");
  const blob=new Blob(["\ufeff"+csv],{type:"text/csv;charset=utf-8"}), url=URL.createObjectURL(blob), link=document.createElement("a");
  link.href=url; link.download=`LCMS_overall_low_similarity_${DATA.project_id||"report"}.csv`; document.body.appendChild(link); link.click(); link.remove(); URL.revokeObjectURL(url);
}
async function loadFeatures(){
  if(!$('featureTable')){ state.features=[]; return; }
  const p=await fetchJson(apiUrl("/api/features")); state.features=p.features||[]; renderFeatures();
}
async function loadAgentAnnotations(){
  try{
    const payload=await fetchJson(apiUrl("/api/agent-annotations"));
    state.agentAnnotations=Array.isArray(payload.annotations)?payload.annotations:[];
  }catch(_error){
    state.agentAnnotations=[];
  }
}
function renderFeatures(){
  if(!$('featureTable'))return;
  const rows=state.features.map((f,i)=>{
    const rt=Number(f.representative_rt??f.rt_apex);
    const peakMz=truePeakMz(f), envelopeMz=envelopeRepresentativeMz(f);
    return `<tr data-saved-feature="1" data-peak="${escapeHtml(f.parent_tic_peak_id||"")}" data-rt="${Number.isFinite(rt)?rt:""}" data-mz="${peakMz}" class="${f.parent_tic_peak_id===state.selectedPeakId&&Math.abs(Number(peakMz)-Number(state.selectedMz||0))<1e-6?'selected':''}"><td>${i+1}</td><td>${f.parent_tic_peak_id||""}</td><td>${nice(rt,4)}</td><td>${nice(peakMz,5)}</td><td>${nice(envelopeMz,5)}</td><td>${f.source||""}</td></tr>`;
  });
  $("featureTable").innerHTML=`<tr><th>#</th><th>parent peak</th><th>RT</th><th>true peak m/z</th><th>envelope m/z</th><th>source</th></tr>${rows.join("")}`;
  document.querySelectorAll("tr[data-saved-feature]").forEach(row=>row.onclick=async()=>{await selectGlobalFeature({parent_tic_peak_id:row.dataset.peak,representative_rt:row.dataset.rt,representative_mz:Number(row.dataset.mz)});});
}
async function saveFeature(){
  const pk=selectedPeak(); if(!pk||!state.selectedMz||!state.xic)return;
  const featureRt=selectedFeatureRt()??Number(pk.rt_apex);
  const group=selectedFeatureGroup();
  const feature={feature_id:`PF_${Date.now()}`, parent_tic_peak_id:pk.tic_peak_id, project_id:DATA.project_id, sample_ids:DATA.sample_ids, representative_rt:featureRt, representative_mz:state.selectedMz, true_peak_mz:state.selectedMz, envelope_representative_mz:envelopeRepresentativeMz(group)||state.selectedMz, quantitation_mz:group?.quantitation_mz??state.selectedMz, mz_tolerance:state.xic.mz_tolerance, rt_start:pk.rt_start, rt_apex:featureRt, rt_end:pk.rt_end, raw_area_by_sample:Object.fromEntries(Object.entries(state.xic.integration_by_sample).map(([s,v])=>[s,v.area])), difference_score:1-pk.peak_consistency_score, peak_consistency_score:pk.peak_consistency_score, spectrum_score:pk.spectrum_score, source:"peak_first_local_analysis", annotation_status:"unannotated", created_at:new Date().toISOString()};
  const p=await fetchJson(apiUrl("/api/features"),{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({features:[...state.features,feature]})}); state.features=p.features||[]; $("featureInfo").textContent=`Saved ${p.saved_count} features`; renderFeatures();
}
function drawAll(){ drawChrom(); renderPeakTable(); drawDetail(); renderMzTable(); drawXic(); renderGlobalFeatureTable(); drawFeatureMap(); drawFeatureMs2(); renderModificationQuantitation(); drawStructureModule(); }
function attachHover(canvas, kind){
  if(!canvas)return;
  canvas.onmousemove=e=>{
    if(kind==="featureMs2"){
      const hit=featureMs2HitAt(canvas,e.clientX,e.clientY);
      canvas.style.cursor=hit?"crosshair":"default";
      showFeatureMs2Tooltip(e,hit);
    }
  };
  canvas.onmouseleave=()=>{ hideDetailTooltip(); canvas.style.cursor="default"; };
}
function attachZoom(canvas, kind){
  let start=null;
  canvas.onmousedown=e=>{ const r=canvas.getBoundingClientRect(); start={x:(e.clientX-r.left)*canvas.width/r.width,y:(e.clientY-r.top)*canvas.height/r.height}; };
  canvas.onmousemove=e=>{
    if(!start){
      if(kind==="chrom"){ const hit=chromHitAt(canvas,e.clientX,e.clientY); canvas.style.cursor=hit?"crosshair":"default"; showChromTooltip(e,hit); }
      if(kind==="detail"){ const hit=spectrumHitAt(canvas,e.clientX,e.clientY); canvas.style.cursor=hit?"pointer":"crosshair"; showDetailTooltip(e,hit); }
      if(kind==="xic"){ const hit=xicHitAt(canvas,e.clientX,e.clientY); canvas.style.cursor=hit?"crosshair":"default"; showXicTooltip(e,hit); }
      if(kind==="featureMap"){ const hit=featureMapHitAt(canvas,e.clientX,e.clientY); canvas.style.cursor=hit?"pointer":"crosshair"; showFeatureMapTooltip(e,hit); }
      return;
    }
    if(!canvas._plot)return;
    hideDetailTooltip();
    const r=canvas.getBoundingClientRect(); const cur={x:(e.clientX-r.left)*canvas.width/r.width,y:(e.clientY-r.top)*canvas.height/r.height};
    drawAll(); drawDragBox(canvas,start,cur);
  };
  canvas.onmouseleave=()=>{ hideDetailTooltip(); canvas.style.cursor="default"; if(start){ start=null; drawAll(); } };
  canvas.onmouseup=e=>{
    if(!start||!canvas._plot)return;
    const r=canvas.getBoundingClientRect(); const end={x:(e.clientX-r.left)*canvas.width/r.width,y:(e.clientY-r.top)*canvas.height/r.height};
    const p=canvas._plot; const dx=Math.abs(end.x-start.x), dy=Math.abs(end.y-start.y);
    if(dx<5&&dy<5){
      if(kind==="chrom"){ const rt=p.fromX(end.x); const hit=DATA.peak_results.find(pk=>pk.rt_start<=rt&&rt<=pk.rt_end); if(hit){state.selectedPeakId=hit.tic_peak_id; state.selectedMz=null; state.selectedFeatureRt=null; state.selectedFeatureGroupId=null; state.xic=null; state.detailZoom=null; state.xicZoom=null;} }
      if(kind==="detail"){ const hit=spectrumHitAt(canvas,e.clientX,e.clientY); if(hit){ selectSpectrumMz(hit.mz); } }
      if(kind==="featureMap"){ const hit=featureMapHitAt(canvas,e.clientX,e.clientY); if(hit){ selectGlobalFeature(hit); } }
    } else {
      const xmin=Math.min(p.fromX(start.x),p.fromX(end.x)), xmax=Math.max(p.fromX(start.x),p.fromX(end.x));
      const ymin=Math.max(0,Math.min(p.fromY(start.y),p.fromY(end.y))), ymax=Math.max(p.fromY(start.y),p.fromY(end.y));
      if(kind==="chrom") state.chromZoom={xmin,xmax,ymin,ymax};
      if(kind==="detail") state.detailZoom={xmin,xmax,ymin,ymax};
      if(kind==="xic") state.xicZoom={xmin,xmax,ymin,ymax};
      if(kind==="featureMap") state.featureMapZoom={xmin,xmax,ymin,ymax};
    }
    start=null; drawAll();
  };
}
async function boot(){
  const comparisonPayload=await fetchJson("/api/comparisons");
  COMPARISONS=comparisonPayload.comparisons||[];
  const queryParams=new URLSearchParams(location.search);
  const requestedComparison=queryParams.get("comparison"), requestedFeatureGroupId=queryParams.get("feature_group_id");
  state.comparison=COMPARISONS.some(item=>item.id===requestedComparison)?requestedComparison:(comparisonPayload.default_comparison || COMPARISONS[0]?.id || "");
  reportSettings=loadReportSettings(); state.xicFull=reportSettings.xicFull;
  const currentComparison=COMPARISONS.find(item=>item.id===state.comparison);
  $("reportTitle").textContent=currentComparison?.label || state.comparison || "分析报告";
  await loadComparison();
  if(requestedFeatureGroupId){
    const item=(DATA.global_feature_groups||[]).find(feature=>featureGroupIds(feature).includes(String(requestedFeatureGroupId)));
    if(item)await selectGlobalFeature({...item,feature_group_id:String(requestedFeatureGroupId)});
  }
}
async function loadComparison(){
  DATA=await fetchJson(apiUrl("/api/bootstrap")); state.selectedPeakId=null; state.selectedMz=null; state.selectedFeatureRt=null; state.selectedFeatureGroupId=null; state.sequenceOverlapSelection=null; state.featureMapZoom=null;
  const typeCounts=(DATA.global_feature_groups||[]).reduce((acc,item)=>{const key=item.difference_type||"unknown";acc[key]=(acc[key]||0)+1;return acc;},{});
  setupReportSettings();
  renderPairControls(true);
  ["chromMode","showAligned","showLocalPeakAligned","offset"].forEach(id=>$(id).oninput=drawAll); $("resetChrom").onclick=()=>{state.chromZoom=null;drawAll();}; $("resetDetail").onclick=()=>{state.detailZoom=null;drawAll();}; $("resetXic").onclick=()=>{ if(state.xic){state.xicZoom=state.xicFull?null:{xmin:Number(state.xic.rt_start),xmax:Number(state.xic.rt_end),ymin:0,ymax:null};} drawAll(); }; $("resetFeatureMap").onclick=()=>{state.featureMapZoom=null;drawAll();};
  $("exportGlobalFeature").onclick=exportGlobalFeatureTable;
  $("toggleXicFull").onclick=async()=>{ state.xicFull=!state.xicFull; reportSettings.xicFull=state.xicFull; saveReportSettings(); updateXicToggleLabel(); if(state.selectedMz){ await selectSpectrumMz(state.selectedMz,state.selectedFeatureRt); } else { drawAll(); } };
  updateXicToggleLabel();
  wireStructureControls();
  attachZoom($("chromCanvas"),"chrom"); attachZoom($("detailCanvas"),"detail"); attachZoom($("xicCanvas"),"xic"); attachZoom($("featureMapCanvas"),"featureMap");
  attachHover($("featureMs2Canvas"),"featureMs2");
  await loadFeatures(); await loadAgentAnnotations(); drawAll(); ensureMsmsData();
}
boot().catch(err=>{document.body.innerHTML=`<pre>${err.stack||err}</pre>`;});
</script>
</body>
</html>
"""


def write_peak_first_sqlite(path: Path, payload: dict[str, object], spectra: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    bootstrap = dict(payload)
    with closing(sqlite3.connect(path)) as connection:
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS peak_first_artifacts (
                artifact_key TEXT PRIMARY KEY,
                payload_json TEXT NOT NULL,
                updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
            )
            """
        )
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS saved_lcms_features (
                feature_id TEXT PRIMARY KEY,
                feature_json TEXT NOT NULL,
                created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
            )
            """
        )
        connection.execute(
            """
            INSERT INTO peak_first_artifacts (artifact_key, payload_json, updated_at)
            VALUES (?, ?, CURRENT_TIMESTAMP)
            ON CONFLICT(artifact_key) DO UPDATE SET
                payload_json = excluded.payload_json,
                updated_at = excluded.updated_at
            """,
            ("bootstrap", json.dumps(bootstrap, ensure_ascii=False, separators=(",", ":"))),
        )
        for sample_id, sample_spectra in spectra.items():
            connection.execute(
                """
                INSERT INTO peak_first_artifacts (artifact_key, payload_json, updated_at)
                VALUES (?, ?, CURRENT_TIMESTAMP)
                ON CONFLICT(artifact_key) DO UPDATE SET
                    payload_json = excluded.payload_json,
                    updated_at = excluded.updated_at
                """,
                (f"spectra:{sample_id}", json.dumps(sample_spectra, ensure_ascii=False, separators=(",", ":"))),
            )
        connection.commit()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate LC-MS Peak-first Compare V2.")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-id", default="lcms_peak_first")
    parser.add_argument("--reference-sample", default="")
    parser.add_argument("--sample-names-json", default="", help="JSON mapping from input file stem to user-facing sample name.")
    parser.add_argument("--sample-contains", action="append", default=[], help="Keep samples whose sample_id contains this text. Can be repeated.")
    parser.add_argument("--include-sample", action="append", default=[], help="Keep an exact sample_id. Can be repeated.")
    parser.add_argument("--sqlite-name", default="lcms_peak_first_compare.sqlite")
    parser.add_argument("--max-peaks-per-scan", type=int, default=500, help="Maximum stored centroid points per scan; 0 keeps all points.")
    parser.add_argument("--spectrum-min-intensity", type=float, default=0.0)
    parser.add_argument("--min-snr", type=float, default=5.0)
    parser.add_argument("--min-area-ratio", type=float, default=0.00001)
    parser.add_argument("--top-n-peaks", type=int, default=0, help="Optional safety limit; 0 analyzes every TIC peak passing the quality thresholds.")
    parser.add_argument("--top-n-mz", type=int, default=40, help="Candidate m/z limit; 0 keeps all centroid bins.")
    parser.add_argument("--top-n-changed-mz", type=int, default=15, help="Changed m/z limit; 0 keeps all candidates for ranking/display.")
    parser.add_argument("--dynamic-background-candidates", action="store_true", help="Use adaptive summed-spectrum and XIC-quality screening for the background candidate pool.")
    parser.add_argument("--candidate-spectral-noise-multiplier", type=float, default=3.0)
    parser.add_argument("--candidate-min-local-tic-ppm", type=float, default=1500.0)
    parser.add_argument("--candidate-max-per-tic", type=int, default=200)
    parser.add_argument("--candidate-min-consecutive-scans", type=int, default=3)
    parser.add_argument("--candidate-min-xic-tic-area-fraction", type=float, default=0.0003)
    parser.add_argument("--max-spectrum-points-per-scan", type=int, default=500, help="Maximum points used for summed spectra; 0 keeps all points.")
    parser.add_argument("--mz-tolerance-da", type=float, default=0.16, help="Half-window used to locate and extract one centroid peak per scan.")
    parser.add_argument("--mz-tolerance-ppm", type=float, default=10.0, help="Centroid clustering and cross-sample peak matching tolerance.")
    parser.add_argument("--mz-tolerance-mode", choices=["da", "ppm"], default="da")
    parser.add_argument("--min-changed-mz-gap-da", type=float, default=0.02)
    parser.add_argument("--feature-presence-relative-area-fraction", type=float, default=0.02)
    parser.add_argument("--feature-common-fold-change-threshold", type=float, default=2.0)
    parser.add_argument("--feature-strong-fold-change-threshold", type=float, default=4.0)
    parser.add_argument("--feature-partial-detection-rank-weight", type=float, default=0.35)
    return parser.parse_args()


def apply_sample_names(
    raw_files: list[object],
    scans_by_sample: dict[str, list[object]],
    mapping_path: str,
) -> tuple[list[object], dict[str, list[object]]]:
    """Replace filename-derived sample IDs with the labels entered on upload."""
    if not mapping_path:
        return raw_files, scans_by_sample
    path = Path(mapping_path)
    if not path.exists():
        raise FileNotFoundError(f"sample names mapping was not found: {path}")
    mapping = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(mapping, dict):
        raise ValueError("sample names mapping must be a JSON object")
    renamed: dict[str, list[object]] = {}
    for raw_file in raw_files:
        source_sample_id = str(raw_file.sample_id)
        file_stem = Path(str(raw_file.file_name)).stem
        label = str(mapping.get(file_stem) or mapping.get(str(raw_file.raw_file_id)) or source_sample_id).strip()
        if not label:
            label = source_sample_id
        if label in renamed:
            raise ValueError(f"sample name is duplicated: {label}")
        raw_file.sample_id = label
        scans = scans_by_sample.get(source_sample_id, [])
        for scan in scans:
            scan.sample_id = label
        renamed[label] = scans
    return raw_files, renamed


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_files, scans_by_sample = load_lcms_directory(input_dir, args.project_id)
    raw_files, scans_by_sample = apply_sample_names(raw_files, scans_by_sample, args.sample_names_json)
    if args.include_sample or args.sample_contains:
        exact = set(args.include_sample or [])
        contains = [str(item) for item in args.sample_contains or [] if str(item)]
        keep = [
            sample_id for sample_id in scans_by_sample
            if sample_id in exact or any(token in sample_id for token in contains)
        ]
        scans_by_sample = {sample_id: scans_by_sample[sample_id] for sample_id in keep}
        raw_files = [raw_file for raw_file in raw_files if raw_file.sample_id in scans_by_sample]
    if len(scans_by_sample) < 2:
        raise ValueError(
            "Peak-first Compare needs at least two selected samples. "
            f"Selected: {list(scans_by_sample)}"
        )
    sample_ids = list(scans_by_sample)
    reference = args.reference_sample if args.reference_sample in scans_by_sample else sample_ids[0]
    params = PeakFirstParams(
        min_snr=args.min_snr,
        min_area_ratio=args.min_area_ratio,
        top_n_peaks=args.top_n_peaks,
        top_n_mz=args.top_n_mz,
        top_n_changed_mz=args.top_n_changed_mz,
        dynamic_background_candidates=args.dynamic_background_candidates,
        candidate_spectral_noise_multiplier=args.candidate_spectral_noise_multiplier,
        candidate_min_local_tic_ppm=args.candidate_min_local_tic_ppm,
        candidate_max_per_tic=args.candidate_max_per_tic,
        candidate_min_consecutive_scans=args.candidate_min_consecutive_scans,
        candidate_min_xic_tic_area_fraction=args.candidate_min_xic_tic_area_fraction,
        max_spectrum_points_per_scan=args.max_spectrum_points_per_scan,
        mz_tolerance_da=args.mz_tolerance_da,
        mz_tolerance_ppm=args.mz_tolerance_ppm,
        mz_tolerance_mode=args.mz_tolerance_mode,
        min_changed_mz_gap_da=args.min_changed_mz_gap_da,
        feature_presence_relative_area_fraction=args.feature_presence_relative_area_fraction,
        feature_common_fold_change_threshold=args.feature_common_fold_change_threshold,
        feature_strong_fold_change_threshold=args.feature_strong_fold_change_threshold,
        feature_partial_detection_rank_weight=args.feature_partial_detection_rank_weight,
    )
    payload = prepare_peak_first_payload(raw_files, scans_by_sample, args.project_id, reference, params)
    spectra = spectrum_payload(
        scans_by_sample,
        shifts=dict(payload["alignment"]["rt_shift_by_sample"]),
        mz_min=0,
        mz_max=10_000,
        max_peaks_per_scan=args.max_peaks_per_scan,
        min_intensity=args.spectrum_min_intensity,
    )
    payload["ms1_component_groups"] = build_ms1_component_groups(
        payload,
        [],
        spectra_loader=lambda sample_id: list(spectra.get(sample_id, [])),
    )
    html_path = output_dir / "lcms_peak_first_compare.html"
    sqlite_path = output_dir / args.sqlite_name
    write_peak_first_sqlite(sqlite_path, payload, spectra)
    html_path.write_text(PEAK_FIRST_TEMPLATE, encoding="utf-8")
    print(f"Samples: {len(sample_ids)}")
    print(f"Reference: {reference}")
    print(f"Confirmed TIC peaks: {len(payload['peak_results'])}")
    print(f"SQLite: {sqlite_path.resolve()}")
    print(f"Peak-first Compare: {html_path.resolve()}")


if __name__ == "__main__":
    main()
