#!/usr/bin/env python3
"""Identify peptide-map MS/MS spectra and link them to Peak-first differential features."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sqlite3
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from core.lcms_msms import (
    annotate_feature_groups,
    assign_q_values,
    build_feature_ms2_evidence,
    build_ms1_component_groups,
    build_modification_level_quantitation,
    build_proteolytic_modification_families,
    build_modified_peptide_findings,
    build_peptide_form_comparisons,
    differential_annotations,
    generate_candidates,
    generate_sequence_inference_candidates,
    read_fasta,
    read_peak_first_payload,
    search_component_consensus_scans,
    search_component_sequence_tag_scans,
    search_feature_glycopeptide_scans,
    search_feature_guided_scans,
    search_feature_open_mass_scans,
    search_scans,
    search_selected_feature_consensus_scans,
)
from core.lcms_parser import read_mzml
from core.lcms_enzymes import DEFAULT_ENZYME, ENZYMES
from core.lcms_unimod import search_unimod_rescue
from core.lcms_sample_prep import ALL_PREP_FIELDS, PREP_CHOICES, normalize_sample_prep, sample_prep_model, filter_linear_candidates


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mzml", action="append", required=True, help="Centroid mzML; repeat for each sample")
    parser.add_argument("--fasta", required=True, help="User-verified heavy/light-chain FASTA")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--feature-sqlite", default="", help="Optional Peak-first SQLite for feature linking")
    parser.add_argument("--project-id", default="vedolizumab_msms")
    parser.add_argument("--enzyme", choices=tuple(ENZYMES), default=DEFAULT_ENZYME)
    for key, label in ALL_PREP_FIELDS.items():
        parser.add_argument(f"--prep-{key}", choices=tuple(PREP_CHOICES[key]), default="unknown",
                            help=label)
    parser.add_argument("--disable-unimod-rescue", action="store_true",
                        help="Disable the separate exploratory second pass (for controlled comparisons)")
    parser.add_argument("--max-missed-cleavages", type=int, default=2)
    parser.add_argument("--max-variable-modifications", type=int, default=2)
    parser.add_argument("--semitryptic-max-trim", type=int, default=2)
    parser.add_argument("--min-peptide-length", type=int, default=6)
    parser.add_argument("--max-peptide-length", type=int, default=60)
    parser.add_argument("--no-carbamidomethyl-cys", action="store_true", help="Disable fixed Cys carbamidomethylation")
    parser.add_argument("--precursor-ppm", type=float, default=10.0)
    parser.add_argument("--fragment-ppm", type=float, default=20.0)
    parser.add_argument("--fdr", type=float, default=0.01)
    parser.add_argument("--min-score", type=float, default=20.0)
    parser.add_argument("--max-report-psms", type=int, default=2000)
    parser.add_argument("--sequence-inference-max-missed-cleavages", type=int, default=4)
    parser.add_argument("--sequence-inference-max-terminal-trim", type=int, default=20)
    parser.add_argument("--sequence-inference-min-delta-da", type=float, default=-250.0)
    parser.add_argument("--sequence-inference-max-delta-da", type=float, default=2500.0)
    parser.add_argument(
        "--ms2-workers",
        type=int,
        default=0,
        help="Parallel sample-level MS2 workers; 0 chooses one worker per sample.",
    )
    args = parser.parse_args()
    if args.no_carbamidomethyl_cys:
        if args.prep_alkylation not in {"unknown", "none"}:
            parser.error("--no-carbamidomethyl-cys conflicts with the selected alkylation reagent")
        args.prep_alkylation = "none"
    args.sample_prep = normalize_sample_prep({key: getattr(args, f"prep_{key}") for key in ALL_PREP_FIELDS})
    return args


def emit_progress(fraction: float, stage: str) -> None:
    """Emit a machine-readable progress update for the desktop worker."""
    print(f"LCMS_PROGRESS\t{max(0.0, min(1.0, fraction)):.4f}\t{stage}", flush=True)


def emit_timing(stage: str, started: float) -> None:
    """Emit a lightweight stage timer for performance diagnostics."""
    print(
        f"LCMS_TIMING\t{stage}\t{time.perf_counter() - started:.3f}",
        flush=True,
    )


def compact_psm(psm: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in psm.items() if key not in {"matched_fragments", "spectrum_peaks"}}


def component_search_excluded_feature_ids(
    annotations: list[dict[str, object]],
) -> set[str]:
    """Exclude established B/C links while leaving D-level links rescuable."""
    return {
        str(annotation.get("feature_or_candidate_id"))
        for annotation in annotations
        if annotation.get("feature_or_candidate_id")
        and annotation.get("confidence") in {
            "B_high_confidence_inferred",
            "C_tentative",
        }
    }


def conversion_metadata(mzml_path: Path) -> dict[str, object]:
    path = mzml_path.with_name(f"{mzml_path.stem}-metadata.json")
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    values = {
        str(item.get("name")): item.get("value")
        for section in payload.values()
        if isinstance(section, list)
        for item in section
        if isinstance(item, dict)
    }
    return {
        "instrument": values.get("Thermo Scientific instrument model"),
        "injection_volume": values.get("injection volume setting"),
        "acquisition_method": values.get("device acquisition method"),
        "metadata_file": str(path.resolve()),
    }


def search_ms2_sample(job: tuple[Path, str, list[object], float, float, float, float]) -> dict[str, object]:
    """Search one sample in a worker process.

    The parser and primary search are CPU-heavy Python work. Running one
    process per sample avoids the GIL while keeping the large decoded scan
    arrays on each worker's side of the process boundary. The parent reloads
    them from the local mzML cache for cross-sample rescue stages.
    """
    mzml_path, project_id, candidates, precursor_ppm, fragment_ppm, min_score, fdr = job
    raw_file, scans = read_mzml(mzml_path, project_id, ms_levels=(2,))
    searched = assign_q_values(
        search_scans(
            scans,
            candidates,
            precursor_tolerance_ppm=max(0.1, precursor_ppm),
            fragment_tolerance_ppm=max(0.1, fragment_ppm),
            min_score=max(0.0, min_score),
        )
    )
    accepted = [
        psm for psm in searched
        if not psm["is_decoy"] and float(psm["q_value"]) <= fdr
    ]
    return {
        "sample_summary": {
            "sample_id": raw_file.sample_id,
            "file": raw_file.file_name,
            "ms2_scans": len(scans),
            "searched_psms": len(searched),
            "accepted_psms": len(accepted),
            "decoy_winners": sum(bool(psm["is_decoy"]) for psm in searched),
            **conversion_metadata(mzml_path),
        },
        "accepted": accepted,
    }


def write_csv(path: Path, psms: list[dict[str, object]]) -> None:
    columns = [
        "sample_id", "scan_id", "rt", "precursor_mz", "precursor_charge", "precursor_error_ppm",
        "chain", "start", "end", "sequence", "proteolysis", "modification_text", "score", "q_value",
        "matched_ion_count", "fragment_coverage", "explained_intensity", "feature_group_id",
        "feature_difference", "feature_ranking", "feature_link_type", "feature_isotope_offset", "search_origin",
        "sequence_inference_level", "backbone_neutral_mass", "mass_delta",
        "localized_mass_offset_site", "sequence_tag_length",
        "complementary_ion_pair_count", "candidate_score_margin",
        "glycan_name", "glycan_composition", "glycan_site",
        "glycan_diagnostic_ion_count", "glycan_diagnostic_intensity_fraction",
        "glycan_core_y_ion_count", "glycan_core_y_type_count",
        "glycan_hexnac_fragment_count",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(psms)


def write_feature_evidence_csv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "rank", "feature_group_id", "parent_tic_peak_id", "representative_rt", "representative_mz",
        "true_peak_mz", "envelope_representative_mz",
        "difference_type", "max_fold_change", "ranking_score", "higher_abundance_sample", "lower_abundance_sample",
        "ms2_status", "sequence", "modification", "confidence", "feature_link_type", "feature_isotope_offset",
        "selected_precursor_charge", "peptide_likelihood", "peptide_likelihood_reason",
        "unidentified_reason", "candidate_mass_hypotheses",
        "sequence_region_candidates",
        "hypothesis", "normalized_area_by_sample", "coverage_by_sample",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            output = dict(row)
            output["normalized_area_by_sample"] = json.dumps(output.get("normalized_area_by_sample") or {}, ensure_ascii=False)
            output["coverage_by_sample"] = json.dumps(output.get("coverage_by_sample") or {}, ensure_ascii=False)
            output["candidate_mass_hypotheses"] = json.dumps(output.get("candidate_mass_hypotheses") or [], ensure_ascii=False)
            output["sequence_region_candidates"] = json.dumps(output.get("sequence_region_candidates") or [], ensure_ascii=False)
            writer.writerow(output)


def write_modified_findings_csv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "rank", "finding_id", "base_peptide_id", "chain", "start", "end", "sequence", "proteolysis",
        "modification", "modification_mass_delta", "neutral_mass", "confidence", "pairing_status",
        "evidence_scope", "linked_feature_ids", "difference_type", "max_fold_change",
        "higher_abundance_sample", "lower_abundance_sample", "normalized_area_by_sample", "sample_evidence",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            output = dict(row)
            output["linked_feature_ids"] = ";".join(row.get("linked_feature_ids") or [])
            output["normalized_area_by_sample"] = json.dumps(output.get("normalized_area_by_sample") or {}, ensure_ascii=False, separators=(",", ":"))
            output["sample_evidence"] = json.dumps(output.get("sample_evidence") or {}, ensure_ascii=False, separators=(",", ":"))
            writer.writerow(output)


def write_modification_pairs_csv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "rank", "pair_id", "base_peptide_id", "chain", "start", "end", "sequence", "modification",
        "modification_mass_delta", "linked_feature_ids", "modified_ms2_confidence", "reference_sample",
        "modified_status", "unmodified_status", "modified_consensus_rt", "unmodified_consensus_rt",
        "modified_area_by_sample", "unmodified_area_by_sample", "modified_unmodified_ratio_by_sample",
        "apparent_occupancy_by_sample", "modified_charge_states_by_sample", "unmodified_charge_states_by_sample",
        "pairwise_comparisons",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            modified = dict(row.get("modified_form") or {})
            unmodified = dict(row.get("unmodified_form") or {})
            output = {
                **row,
                "linked_feature_ids": ";".join(row.get("linked_feature_ids") or []),
                "modified_status": modified.get("status"),
                "unmodified_status": unmodified.get("status"),
                "modified_consensus_rt": modified.get("consensus_rt"),
                "unmodified_consensus_rt": unmodified.get("consensus_rt"),
                "modified_area_by_sample": modified.get("normalized_area_by_sample"),
                "unmodified_area_by_sample": unmodified.get("normalized_area_by_sample"),
                "modified_charge_states_by_sample": modified.get("charge_states_by_sample"),
                "unmodified_charge_states_by_sample": unmodified.get("charge_states_by_sample"),
            }
            for key in (
                "modified_area_by_sample", "unmodified_area_by_sample", "modified_unmodified_ratio_by_sample",
                "apparent_occupancy_by_sample", "modified_charge_states_by_sample",
                "unmodified_charge_states_by_sample", "pairwise_comparisons",
            ):
                output[key] = json.dumps(output.get(key), ensure_ascii=False, separators=(",", ":"))
            writer.writerow(output)


def write_modification_level_csv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "rank", "quantitation_id", "event_type", "event_name", "chain", "start", "end",
        "sequence", "quantitation_status", "denominator_scope", "relative_ms_response_only",
        "max_percentage_point_difference", "high_value_candidate", "linked_feature_ids",
        "forms", "total_area_by_sample", "pairwise_comparisons",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            output = dict(row)
            output["linked_feature_ids"] = ";".join(row.get("linked_feature_ids") or [])
            for key in ("forms", "total_area_by_sample", "pairwise_comparisons"):
                output[key] = json.dumps(output.get(key), ensure_ascii=False, separators=(",", ":"))
            writer.writerow(output)


def read_sqlite_spectra(path: Path, sample_id: str) -> list[dict[str, object]]:
    with sqlite3.connect(path) as connection:
        row = connection.execute(
            "SELECT payload_json FROM peak_first_artifacts WHERE artifact_key = ?",
            (f"spectra:{sample_id}",),
        ).fetchone()
    if row is None:
        return []
    payload = json.loads(str(row[0]))
    return payload if isinstance(payload, list) else []


def write_sqlite_artifact(path: Path, report: dict[str, object]) -> None:
    compact = dict(report)
    compact["psms"] = [compact_psm(psm) for psm in report["psms"]]
    compact["sequence_region_candidates"] = [
        compact_psm(candidate)
        for candidate in report.get("sequence_region_candidates") or []
    ]
    component_payload = {
        "schema_version": "LCMSMSComponentGroups/v1",
        "component_groups": report.get("ms1_component_groups") or [],
    }
    with sqlite3.connect(path) as connection:
        connection.execute(
            """
            INSERT INTO peak_first_artifacts (artifact_key, payload_json, updated_at)
            VALUES ('msms_identifications', ?, CURRENT_TIMESTAMP)
            ON CONFLICT(artifact_key) DO UPDATE SET
                payload_json = excluded.payload_json,
                updated_at = excluded.updated_at
            """,
            (json.dumps(compact, ensure_ascii=False, separators=(",", ":")),),
        )
        connection.execute(
            """
            INSERT INTO peak_first_artifacts (artifact_key, payload_json, updated_at)
            VALUES ('ms1_component_groups', ?, CURRENT_TIMESTAMP)
            ON CONFLICT(artifact_key) DO UPDATE SET
                payload_json = excluded.payload_json,
                updated_at = excluded.updated_at
            """,
            (json.dumps(component_payload, ensure_ascii=False, separators=(",", ":")),),
        )
        connection.commit()


def html_report(report: dict[str, object]) -> str:
    data = json.dumps(report, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>LC-MS/MS 差异组分鉴定</title><style>
body{{margin:0;font:14px Arial,"Microsoft YaHei",sans-serif;color:#172033;background:#f5f7fb}}header{{padding:16px 22px;background:#0f172a;color:white}}main{{padding:14px;display:grid;gap:12px}}section{{background:white;border:1px solid #dbe2ec;border-radius:8px;padding:12px}}h1,h2{{margin:0 0 8px}}.summary{{display:flex;gap:18px;flex-wrap:wrap}}input{{padding:6px 8px;border:1px solid #ccd5e1;border-radius:6px;width:420px}}.scroll{{overflow:auto;max-height:520px}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{padding:6px;border-bottom:1px solid #e8edf4;white-space:nowrap;text-align:right}}th:nth-child(-n+3),td:nth-child(-n+3),th:nth-child(7),td:nth-child(7){{text-align:left}}td.wrap{{white-space:normal;min-width:280px;text-align:left}}tr[data-index],tr[data-feature-index]{{cursor:pointer}}tr.selected{{background:#fff0c7}}canvas{{width:100%;height:330px;border:1px solid #e3e8f0;border-radius:6px}}.muted{{color:#667085}}.badge{{padding:2px 6px;border-radius:10px;background:#eef2ff}}.direct{{color:#067647;font-weight:600}}.isotope{{color:#175cd3;font-weight:600}}.tentative{{color:#b54708;font-weight:600}}.unresolved{{color:#667085}}
</style></head><body><header><h1>LC-MS/MS 差异组分鉴定与推测</h1><div id="summary" class="summary"></div></header><main>
<section><h2>MS1 显著差异 → MS2 定性证据</h2><input id="featureFilter" placeholder="筛选 Feature、序列、修饰、MS2 状态或高丰度样品"><div id="featureSummary" class="muted" style="margin:8px 0"></div><div class="scroll"><table id="featureTable"></table></div></section>
  <section><h2>修饰肽—未修饰母肽配对</h2><input id="pairFilter" placeholder="筛选序列、修饰或 Feature"><div id="pairSummary" class="muted" style="margin:8px 0"></div><div class="scroll"><table id="pairTable"></table></div><div class="muted" style="margin-top:8px">同一肽的不同电荷态和 M/M+1/M+2 同位素包络已合并。表观修饰占比用于样品间半定量比较，不等同于经过标准品校正的绝对占有率。</div></section>
  <section><h2>位点修饰水平与蛋白形式差异定量</h2><input id="levelFilter" placeholder="筛选序列、事件、修饰或 Feature"><div id="levelSummary" class="muted" style="margin:8px 0"></div><div class="scroll"><table id="levelTable"></table></div><div class="muted" style="margin-top:8px">仅 B/C 级 MS2 支持的形式可进入本表；B 级完整形式族进入正式相对定量，C 级单列为暂定。所有百分比均为相对 MS 响应构成。</div></section>
  <section><h2>鉴定证据</h2><input id="filter" placeholder="筛选样品、序列、修饰或 Feature"><div class="scroll"><table id="table"></table></div></section>
<section><h2>MS/MS 证据谱</h2><div id="detail" class="muted">点击上表查看谱图</div><canvas id="spectrum" width="1500" height="360"></canvas></section>
<section><h2>说明</h2><div id="warnings" style="color:#b42318;margin-bottom:6px"></div><div class="muted">B = target-decoy q≤1%、得分/覆盖度满足阈值的高可信推定；C = 部分鉴定；D = 低证据。A 级确认仍需标准品或正交方法。单蛋白小数据库的 FDR 仅作探索性控制。</div></section>
<section><h2>所有已定性修饰（不要求未修饰配对）</h2><input id="modFilter" placeholder="筛选序列、修饰、Feature 或配对状态"><div id="modSummary" class="muted" style="margin:8px 0"></div><div class="scroll"><table id="modTable"></table></div><div class="muted" style="margin-top:8px">MS2-only 表示修饰肽已获得序列证据，但尚未与显著 MS1 Feature 或未修饰肽完成定量配对；PSM 数和前体强度仅作为采集证据，不作为修饰占有率。</div></section>
</main><script id="data" type="application/json">{data}</script><script>
  const DATA=JSON.parse(document.getElementById('data').textContent),$=id=>document.getElementById(id);let selected=-1,selectedFeature=-1;
  $('summary').innerHTML=`<span>${{DATA.samples.length}} samples</span><span>${{DATA.total_ms2_scans}} MS2 scans</span><span>${{DATA.standard_accepted_psms??DATA.total_accepted_psms}} standard PSMs</span><span>${{DATA.component_consensus_accepted_psms||0}} component-consensus PSMs</span><span>${{DATA.component_sequence_tag_accepted_psms||0}} backbone/tag PSMs</span><span>${{DATA.sequence_region_candidate_count||0}} sequence-region candidates</span><span>${{DATA.consensus_feature_accepted_psms||0}} single-Feature consensus PSMs</span><span>${{DATA.feature_guided_accepted_psms||0}} feature-guided PSMs</span><span>${{DATA.feature_glycopeptide_form_count||0}} targeted glycopeptide forms / ${{DATA.feature_glycopeptide_feature_count||0}} Features</span><span>${{DATA.feature_open_mass_accepted_psms||0}} feature open-mass PSMs</span><span>${{DATA.feature_evidence?.length||0}} significant MS1 Features</span><span>${{DATA.modified_peptide_findings?.length||0}} modified findings</span><span>${{DATA.peptide_form_comparisons?.length||0}} modification pairs</span><span>${{DATA.formal_modification_level_count||0}} formal modification levels</span>`;
$('warnings').textContent=[...(DATA.warnings||[]),...(DATA.parameters?.sample_prep_model?.warnings||[])].join('；');
const nice=(v,n=3)=>v==null?'':Number(v).toFixed(n);function rows(){{const q=$('filter').value.toLowerCase();return DATA.psms.map((p,i)=>[p,i]).filter(([p])=>!q||[p.sample_id,p.sequence,p.modification_text,p.feature_group_id].join(' ').toLowerCase().includes(q));}}
const statusText={{identified_direct_precursor:'直接前体鉴定',identified_isotope_envelope:'同位素包络对应',identified_charge_state_envelope:'电荷态包络对应',tentative_component_consensus_identification:'组件多电荷/同位素联合推测',tentative_backbone_sequence_support:'未知质量偏移骨架序列推测',tentative_truncation_sequence_support:'截断骨架序列推测',tentative_sequence_region_candidate:'局部 b/y 序列区段候选',tentative_feature_consensus_identification:'单 Feature 共识谱推测',tentative_feature_guided_identification:'MS1 引导定向推测',tentative_feature_open_mass_identification:'单 Feature 开放质量搜索',tentative_feature_glycopeptide_identification:'差异 Feature 定向糖肽推测',low_evidence_sequence_candidate:'低证据序列候选（不替代 Unknown）',selected_precursor_unidentified:'已采集前体，未鉴定',coisolated_ms2_unresolved:'隔离窗覆盖，未定性',no_ms2_acquired:'未采集 MS2'}};
const statusClass=s=>s==='identified_direct_precursor'?'direct':['identified_isotope_envelope','identified_charge_state_envelope'].includes(s)?'isotope':String(s||'').startsWith('tentative_')?'tentative':'unresolved';
const effectiveFeatureStatus=f=>String(f?.confidence||'').startsWith('D_')?'low_evidence_sequence_candidate':f?.ms2_status;
function featureRows(){{const q=$('featureFilter').value.toLowerCase();return (DATA.feature_evidence||[]).map((f,i)=>[f,i]).filter(([f])=>!q||[f.feature_group_id,f.sequence,f.modification,f.ms2_status,f.higher_abundance_sample,f.hypothesis].join(' ').toLowerCase().includes(q));}}
  function coverageText(f){{return Object.entries(f.coverage_by_sample||{{}}).map(([sample,v])=>`${{sample}}: window ${{v.isolation_window_scan_count}}, selected ${{v.selected_precursor_scan_count}}, isotope ${{v.selected_isotope_scan_count}}`).join('；');}}
  function renderFeatures(){{const summary=DATA.feature_ms2_summary||{{}};$('featureSummary').textContent=Object.entries(summary).map(([k,v])=>`${{statusText[k]||k}} ${{v}}`).join('；');$('featureTable').innerHTML='<tr><th>rank</th><th>Feature</th><th>TIC peak</th><th>RT</th><th>true peak m/z</th><th>envelope m/z</th><th>MS1 type</th><th>fold</th><th>higher sample</th><th>MS2 status</th><th>sequence</th><th>modification / change</th><th>confidence</th><th>MS2 coverage</th><th>interpretation</th></tr>'+featureRows().map(([f,i])=>{{const s=effectiveFeatureStatus(f);return `<tr data-feature-index="${{i}}" class="${{i===selectedFeature?'selected':''}}"><td>${{f.rank}}</td><td>${{f.feature_group_id}}</td><td>${{f.parent_tic_peak_id||''}}</td><td>${{nice(f.representative_rt)}}</td><td>${{nice(f.true_peak_mz??f.representative_mz,5)}}</td><td>${{nice(f.envelope_representative_mz??f.true_peak_mz??f.representative_mz,5)}}</td><td>${{f.difference_type||''}}</td><td>${{nice(f.max_fold_change,2)}}</td><td>${{f.higher_abundance_sample||''}}</td><td class="${{statusClass(s)}}">${{statusText[s]||s}}</td><td>${{f.sequence||''}}</td><td>${{f.modification||''}}</td><td>${{f.confidence||''}}</td><td class="wrap">${{coverageText(f)}}</td><td class="wrap">${{f.hypothesis||''}}</td></tr>`;}}).join('');document.querySelectorAll('tr[data-feature-index]').forEach(tr=>tr.onclick=()=>{{selectedFeature=Number(tr.dataset.featureIndex);renderFeatures();const f=DATA.feature_evidence[selectedFeature],e=f.best_psm||f.candidate_psm||f.exploratory_psm||f.coverage_scan,s=[f.best_psm,f.candidate_psm,f.exploratory_psm,f.coverage_scan].find(item=>Array.isArray(item?.spectrum_peaks)&&item.spectrum_peaks.length)||e;if(e)draw({{...e,spectrum_peaks:s?.spectrum_peaks||[],spectrum_source_scan_id:s?.scan_id,feature_group_id:f.feature_group_id,hypothesis:f.exploratory_psm?.exploratory_note||f.hypothesis}});}});}}
  const formStatusText={{ms2_confirmed:'MS2确认',ms1_multi_charge_supported:'MS1多电荷态支持',ms1_candidate:'MS1候选',not_detected:'未检出'}};
  const directionText={{increased:'供试组修饰比升高',decreased:'供试组修饰比降低',stable:'修饰比基本不变',unavailable:'无法计算'}};
  function sampleMap(values,percent=false){{return Object.entries(values||{{}}).map(([sample,value])=>`${{sample}}: ${{value==null?'NA':percent?(Number(value)*100).toFixed(2)+'%':nice(value,4)}}`).join('；');}}
  function chargeMap(form){{return Object.entries(form?.charge_states_by_sample||{{}}).map(([sample,charges])=>`${{sample}}: ${{charges.length?charges.map(z=>'z'+z).join(','):'未检出'}}`).join('；');}}
  function pairRows(){{const q=$('pairFilter').value.toLowerCase();return (DATA.peptide_form_comparisons||[]).filter(p=>!q||[p.sequence,p.modification,...(p.linked_feature_ids||[])].join(' ').toLowerCase().includes(q));}}
  function renderPairs(){{const pairs=pairRows();$('pairSummary').textContent=`共 ${{DATA.peptide_form_comparisons?.length||0}} 个修饰肽组，已将同位素峰和不同电荷态归并；当前显示 ${{pairs.length}} 个。`;$('pairTable').innerHTML='<tr><th>rank</th><th>sequence</th><th>modification</th><th>Δmass</th><th>linked Features</th><th>modified form</th><th>unmodified form</th><th>modified/unmodified ratio</th><th>apparent occupancy</th><th>reference comparison</th><th>charge-trend support</th><th>MS2 confidence</th></tr>'+pairs.map(p=>{{const c=(p.pairwise_comparisons||[])[0]||{{}},mt=c.modified_charge_trend||{{}};return `<tr><td>${{p.rank}}</td><td>${{p.sequence}}</td><td>${{p.modification}}</td><td>${{nice(p.modification_mass_delta,4)}}</td><td class="wrap">${{(p.linked_feature_ids||[]).join(', ')}}</td><td class="wrap">${{formStatusText[p.modified_form?.status]||p.modified_form?.status}}; RT ${{nice(p.modified_form?.consensus_rt)}}; ${{chargeMap(p.modified_form)}}</td><td class="wrap">${{formStatusText[p.unmodified_form?.status]||p.unmodified_form?.status}}; RT ${{nice(p.unmodified_form?.consensus_rt)}}; ${{chargeMap(p.unmodified_form)}}</td><td class="wrap">${{sampleMap(p.modified_unmodified_ratio_by_sample)}}</td><td class="wrap">${{sampleMap(p.apparent_occupancy_by_sample,true)}}</td><td class="wrap">${{c.reference_sample||''}} → ${{c.test_sample||''}}: ${{directionText[c.modification_direction]||c.modification_direction||''}}; ratio fold ${{nice(c.modified_unmodified_ratio_fold,2)}}; occupancy Δ ${{c.occupancy_difference==null?'NA':(Number(c.occupancy_difference)*100).toFixed(2)+'%'}}</td><td>${{mt.status||'unavailable'}} (${{mt.agreeing_charge_count||0}}/${{mt.comparable_charge_count||0}})</td><td>${{p.modified_ms2_confidence||''}}</td></tr>`;}}).join('');}}
  function levelRows(){{const q=$('levelFilter').value.toLowerCase();return (DATA.modification_level_quantitation||[]).filter(r=>!q||[r.event_name,r.sequence,r.chain,...(r.linked_feature_ids||[]),...(r.forms||[]).map(f=>f.label)].join(' ').toLowerCase().includes(q));}}
  function renderLevels(){{const rows=levelRows(),formal=(DATA.modification_level_quantitation||[]).filter(r=>r.quantitation_status==='formal_relative_quantitation').length,high=(DATA.modification_level_quantitation||[]).filter(r=>r.high_value_candidate).length;$('levelSummary').textContent=`共 ${{DATA.modification_level_quantitation?.length||0}} 个形式族；正式相对定量 ${{formal}} 个；高价值候选 ${{high}} 个；当前显示 ${{rows.length}} 个。`;$('levelTable').innerHTML='<tr><th>rank</th><th>event</th><th>chain / site</th><th>sequence</th><th>form levels</th><th>largest difference</th><th>status</th><th>linked Features</th></tr>'+rows.map(r=>{{const forms=(r.forms||[]).filter(f=>f.included_in_denominator).map(f=>`${{f.label}}：${{sampleMap(f.relative_level_by_sample,true)}}`).join('<br>'),strongest=(r.pairwise_comparisons||[]).map(c=>c.strongest_change).filter(Boolean).sort((a,b)=>Math.abs(b.percentage_point_difference)-Math.abs(a.percentage_point_difference))[0];return `<tr><td>${{r.rank}}</td><td>${{r.event_name}}${{r.high_value_candidate?'；高价值候选':''}}</td><td>${{r.chain||''}} ${{r.start||''}}-${{r.end||''}}</td><td>${{r.sequence||''}}</td><td class="wrap">${{forms||'互补形式不足'}}</td><td class="wrap">${{strongest?strongest.label+'：'+(strongest.percentage_point_difference*100).toFixed(2)+' 个百分点':'NA'}}</td><td>${{r.quantitation_status}}</td><td class="wrap">${{(r.linked_feature_ids||[]).join(', ')}}</td></tr>`;}}).join('');}}
function modifiedRows(){{const q=$('modFilter').value.toLowerCase();return (DATA.modified_peptide_findings||[]).filter(f=>!q||[f.sequence,f.modification,f.pairing_status,f.evidence_scope,...(f.linked_feature_ids||[])].join(' ').toLowerCase().includes(q));}}
const originalDraw=draw;draw=function(p){{const fallback=(DATA.feature_evidence||[]).find(f=>f.feature_group_id===p.feature_group_id)?.coverage_scan?.spectrum_peaks||[];if((!p.spectrum_peaks||!p.spectrum_peaks.length)&&fallback.length)p={{...p,spectrum_peaks:fallback}};return originalDraw(p);}};
function sampleEvidenceText(f){{return Object.entries(f.sample_evidence||{{}}).map(([sample,v])=>`${{sample}}: PSM ${{v.psm_count}}, best ${{nice(v.best_score,1)}}, z ${{(v.charge_states||[]).join('/')}}`).join('；');}}
function renderModified(){{const findings=modifiedRows(),all=DATA.modified_peptide_findings||[];const linked=findings.filter(f=>f.evidence_scope==='differential_feature_linked').length,paired=findings.filter(f=>f.pairing_status==='paired_with_unmodified_form').length;$('modSummary').textContent=`共 ${{all.length}} 个修饰肽发现；当前显示 ${{findings.length}} 个，其中差异 Feature 关联 ${{linked}} 个、已配对 ${{paired}} 个。`;$('modTable').innerHTML='<tr><th>rank</th><th>sequence</th><th>modification</th><th>proteolysis</th><th>confidence</th><th>evidence</th><th>pairing</th><th>linked Features</th><th>MS1 direction</th><th>sample MS2 evidence</th></tr>'+findings.map(f=>`<tr data-mod-index="${{all.indexOf(f)}}"><td>${{f.rank}}</td><td>${{f.sequence}}</td><td>${{f.modification}}</td><td>${{f.proteolysis||''}}</td><td>${{f.confidence}}</td><td>${{f.evidence_scope}}</td><td>${{f.pairing_status}}</td><td class="wrap">${{(f.linked_feature_ids||[]).join(', ')}}</td><td>${{f.higher_abundance_sample?f.higher_abundance_sample+' high; fold '+nice(f.max_fold_change,2):''}}</td><td class="wrap">${{sampleEvidenceText(f)}}</td></tr>`).join('');document.querySelectorAll('tr[data-mod-index]').forEach(tr=>tr.onclick=()=>{{const f=all[Number(tr.dataset.modIndex)];if(f?.best_psm)draw(f.best_psm);}});}}
function render(){{$('table').innerHTML='<tr><th>sample</th><th>scan</th><th>Feature</th><th>RT</th><th>m/z</th><th>z</th><th>sequence</th><th>proteolysis</th><th>modification</th><th>ppm</th><th>ions</th><th>coverage</th><th>score</th><th>q</th></tr>'+rows().map(([p,i])=>`<tr data-index="${{i}}" class="${{i===selected?'selected':''}}"><td>${{p.sample_id}}</td><td>${{p.scan_id}}</td><td>${{p.feature_group_id||''}}</td><td>${{nice(p.rt)}}</td><td>${{nice(p.precursor_mz,5)}}</td><td>${{p.precursor_charge}}</td><td>${{p.sequence}}</td><td>${{p.proteolysis||''}}</td><td>${{p.modification_text}}</td><td>${{nice(p.precursor_error_ppm,2)}}</td><td>${{p.matched_ion_count}}</td><td>${{nice(p.fragment_coverage,2)}}</td><td>${{nice(p.score,1)}}</td><td>${{nice(p.q_value,4)}}</td></tr>`).join('');document.querySelectorAll('tr[data-index]').forEach(tr=>tr.onclick=()=>{{selected=Number(tr.dataset.index);render();draw(DATA.psms[selected]);}});}}
  function ionColor(label){{const text=String(label||'');if(/^b[0-9]/.test(text))return '#2563eb';if(/^y[0-9]/.test(text))return '#dc2626';return '#7c3aed';}}
  function draw(p){{const c=$('spectrum'),x=c.getContext('2d'),peaks=p.spectrum_peaks||[];x.clearRect(0,0,c.width,c.height);const sequence=p.sequence?(p.exploratory_only?`探索性候选：${{p.chain||''}}:${{p.start||''}}-${{p.end||''}} ${{p.sequence}} | ${{p.modification_text||''}}`:`${{p.chain||''}}:${{p.start||''}}-${{p.end||''}} ${{p.sequence}} | ${{p.modification_text||''}}`):'未获得合格序列鉴定',glycan=p.glycan_name?` | 糖链诊断离子 ${{p.glycan_diagnostic_ion_count||0}} | 核心 Y 离子 ${{p.glycan_core_y_ion_count||0}} | HexNAc 保留碎片 ${{p.glycan_hexnac_fragment_count||0}}`:'';const preview=p.exploratory_only?` | 仅供参考（${{p.exploratory_confidence||'low'}}） | b/y ${{p.matched_by_ion_count??0}}/${{nice(p.exploratory_by_coverage??p.fragment_coverage,3)}} | 解释强度 ${{nice(p.explained_intensity,3)}}`:'';$('detail').textContent=`${{p.feature_group_id||''}} | ${{p.sample_id||''}} | ${{p.scan_id||''}} | precursor ${{nice(p.precursor_mz,5)}} z${{p.precursor_charge||''}} | ${{sequence}} | score ${{nice(p.score,1)}} | q ${{nice(p.q_value,4)}}${{preview}}${{glycan}}${{p.hypothesis?' | '+p.hypothesis:''}}`;if(!peaks.length)return;const pad={{l:65,r:20,t:20,b:45}},maxMz=Math.max(...peaks.map(v=>v.mz)),minMz=Math.min(...peaks.map(v=>v.mz)),maxI=Math.max(...peaks.map(v=>v.intensity));x.strokeStyle='#94a3b8';x.beginPath();x.moveTo(pad.l,pad.t);x.lineTo(pad.l,c.height-pad.b);x.lineTo(c.width-pad.r,c.height-pad.b);x.stroke();peaks.forEach(v=>{{const px=pad.l+(v.mz-minMz)/Math.max(1e-9,maxMz-minMz)*(c.width-pad.l-pad.r),py=c.height-pad.b-v.intensity/maxI*(c.height-pad.t-pad.b);const color=v.label?ionColor(v.label):'#475569';x.strokeStyle=color;x.beginPath();x.moveTo(px,c.height-pad.b);x.lineTo(px,py);x.stroke();if(v.label){{x.fillStyle=color;x.font='11px Arial';x.fillText(v.label,px+2,Math.max(12,py-3));}}}});x.fillStyle='#475569';x.fillText(nice(minMz,1),pad.l,c.height-18);x.fillText(nice(maxMz,1),c.width-pad.r-35,c.height-18);}}
  $('featureFilter').oninput=renderFeatures;$('pairFilter').oninput=renderPairs;$('levelFilter').oninput=renderLevels;$('modFilter').oninput=renderModified;$('filter').oninput=render;renderFeatures();renderPairs();renderLevels();renderModified();render();
</script></body></html>"""


def main() -> None:
    args = parse_args()
    prep_model = sample_prep_model(args.sample_prep)
    overall_started = time.perf_counter()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    emit_progress(0.01, "validate_fasta")
    stage_started = time.perf_counter()
    chains = read_fasta(Path(args.fasta))
    candidates = generate_candidates(
        chains,
        max_missed_cleavages=args.max_missed_cleavages,
        min_length=args.min_peptide_length,
        max_length=args.max_peptide_length,
        carbamidomethyl_cys=prep_model["carbamidomethyl_cys"],
        max_variable_modifications=max(0, args.max_variable_modifications),
        semitryptic_max_trim=max(0, args.semitryptic_max_trim),
        enzyme=args.enzyme,
    )
    sequence_inference_candidates = generate_sequence_inference_candidates(
        chains,
        max_missed_cleavages=max(
            args.max_missed_cleavages,
            args.sequence_inference_max_missed_cleavages,
        ),
        min_length=args.min_peptide_length,
        max_length=args.max_peptide_length,
        carbamidomethyl_cys=prep_model["carbamidomethyl_cys"],
        max_terminal_trim=max(
            0,
            args.sequence_inference_max_terminal_trim,
        ),
        enzyme=args.enzyme,
    )
    candidates = filter_linear_candidates(candidates, args.sample_prep)
    sequence_inference_candidates = filter_linear_candidates(sequence_inference_candidates, args.sample_prep)
    emit_progress(0.08, "target_candidates")
    emit_timing("candidate_generation", stage_started)
    sample_summaries: list[dict[str, object]] = []
    accepted: list[dict[str, object]] = []
    sample_jobs = [
        (
            Path(mzml_text),
            args.project_id,
            candidates,
            args.precursor_ppm,
            args.fragment_ppm,
            args.min_score,
            args.fdr,
        )
        for mzml_text in args.mzml
    ]
    worker_count = max(1, min(
        len(sample_jobs),
        args.ms2_workers if args.ms2_workers > 0 else (os.cpu_count() or 1),
    ))
    stage_started = time.perf_counter()
    if worker_count > 1:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            sample_results = list(executor.map(search_ms2_sample, sample_jobs))
    else:
        sample_results = [search_ms2_sample(job) for job in sample_jobs]
    for result in sample_results:
        sample_summaries.append(dict(result["sample_summary"]))
        accepted.extend(list(result["accepted"]))
    emit_progress(0.44, "standard_ms2")
    emit_timing("standard_ms2", stage_started)

    # The worker processes have already populated the per-level mzML cache.
    # Reload MS2 once in the parent for the cross-sample feature/rescue stages.
    stage_started = time.perf_counter()
    all_ms2_scans = []
    for mzml_text in args.mzml:
        _, scans = read_mzml(Path(mzml_text), args.project_id, ms_levels=(2,))
        all_ms2_scans.extend(scans)
    total_ms2_scans = len(all_ms2_scans)
    emit_progress(0.50, "cross_sample_load")
    emit_timing("cross_sample_load", stage_started)
    feature_sqlite = Path(args.feature_sqlite).resolve() if args.feature_sqlite else None
    feature_payload = read_peak_first_payload(feature_sqlite) if feature_sqlite else None
    spectra_cache: dict[str, list[dict[str, object]]] = {}

    def load_cached_spectra(sample_id: str) -> list[dict[str, object]]:
        key = str(sample_id)
        if key not in spectra_cache:
            spectra_cache[key] = (
                read_sqlite_spectra(feature_sqlite, key)
                if feature_sqlite is not None
                else []
            )
        return spectra_cache[key]

    spectra_loader = load_cached_spectra if feature_sqlite is not None else None
    targeted: list[dict[str, object]] = []
    targeted_glycopeptides: list[dict[str, object]] = []
    open_mass: list[dict[str, object]] = []
    consensus: list[dict[str, object]] = []
    component_consensus: list[dict[str, object]] = []
    component_sequence_tags: list[dict[str, object]] = []
    sequence_region_candidates: list[dict[str, object]] = []
    provisional_component_groups: list[dict[str, object]] = []
    if feature_payload is not None:
        annotate_feature_groups(accepted, feature_payload)
        initial_annotations = differential_annotations(accepted)
        stage_started = time.perf_counter()
        provisional_component_groups = build_ms1_component_groups(
            feature_payload,
            initial_annotations,
            spectra_loader=spectra_loader,
        )
        emit_progress(0.56, "component_grouping")
        emit_timing("provisional_component_grouping", stage_started)
        # Keep low-evidence (D) single-spectrum links eligible for component
        # consensus rescue. Only Features already supported at B/C confidence
        # should be excluded from the joint charge/isotope search.
        linked_feature_ids = component_search_excluded_feature_ids(initial_annotations)
        stage_started = time.perf_counter()
        component_consensus = search_component_consensus_scans(
            feature_payload,
            all_ms2_scans,
            candidates,
            provisional_component_groups,
            exclude_feature_ids=linked_feature_ids,
            component_mass_tolerance_ppm=max(0.1, args.precursor_ppm * 2.0),
            fragment_tolerance_ppm=max(0.1, args.fragment_ppm),
            fdr_threshold=args.fdr,
        )
        emit_progress(0.63, "component_consensus")
        emit_timing("component_consensus", stage_started)
        linked_feature_ids.update(
            str(link.get("feature_group_id"))
            for psm in component_consensus
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        )
        stage_started = time.perf_counter()
        (
            component_sequence_tags,
            sequence_region_candidates,
        ) = search_component_sequence_tag_scans(
            feature_payload,
            all_ms2_scans,
            sequence_inference_candidates,
            provisional_component_groups,
            exclude_feature_ids=linked_feature_ids,
            fragment_tolerance_ppm=max(0.1, args.fragment_ppm),
            precursor_relation_ppm=max(0.1, args.precursor_ppm * 2.0),
            minimum_mass_delta_da=args.sequence_inference_min_delta_da,
            maximum_mass_delta_da=args.sequence_inference_max_delta_da,
            fdr_threshold=args.fdr,
        )
        emit_progress(0.70, "sequence_tag_search")
        emit_timing("sequence_tag_search", stage_started)
        linked_feature_ids.update(
            str(link.get("feature_group_id"))
            for psm in component_sequence_tags
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        )
        stage_started = time.perf_counter()
        consensus = search_selected_feature_consensus_scans(
            feature_payload,
            all_ms2_scans,
            candidates,
            exclude_feature_ids=linked_feature_ids,
            precursor_tolerance_ppm=max(0.1, args.precursor_ppm),
            fragment_tolerance_ppm=max(0.1, args.fragment_ppm),
            fdr_threshold=args.fdr,
        )
        emit_progress(0.76, "feature_consensus")
        emit_timing("feature_consensus", stage_started)
        linked_feature_ids.update(
            str(link.get("feature_group_id"))
            for psm in consensus
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        )
        stage_started = time.perf_counter()
        targeted = search_feature_guided_scans(
            feature_payload,
            all_ms2_scans,
            candidates,
            exclude_feature_ids=linked_feature_ids,
            precursor_tolerance_ppm=max(0.1, args.precursor_ppm),
            fragment_tolerance_ppm=max(0.1, args.fragment_ppm),
            fdr_threshold=args.fdr,
        )
        emit_progress(0.81, "feature_guided")
        emit_timing("feature_guided", stage_started)
        linked_feature_ids.update(
            str(link.get("feature_group_id"))
            for psm in targeted
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        )
        stage_started = time.perf_counter()
        targeted_glycopeptides = search_feature_glycopeptide_scans(
            feature_payload,
            all_ms2_scans,
            candidates,
            exclude_feature_ids=linked_feature_ids,
            precursor_tolerance_ppm=max(0.1, args.precursor_ppm * 2.0),
            fragment_tolerance_ppm=max(0.1, args.fragment_ppm),
            fdr_threshold=args.fdr,
            enzyme=args.enzyme,
        )
        emit_progress(0.86, "glycopeptide_search")
        emit_timing("glycopeptide_search", stage_started)
        linked_feature_ids.update(
            str(link.get("feature_group_id"))
            for psm in targeted_glycopeptides
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        )
        stage_started = time.perf_counter()
        open_mass = search_feature_open_mass_scans(
            feature_payload,
            all_ms2_scans,
            sequence_inference_candidates,
            exclude_feature_ids=linked_feature_ids,
            precursor_tolerance_ppm=max(0.1, args.precursor_ppm * 2.0),
            fragment_tolerance_ppm=max(0.1, args.fragment_ppm),
            minimum_mass_delta_da=args.sequence_inference_min_delta_da,
            maximum_mass_delta_da=args.sequence_inference_max_delta_da,
            fdr_threshold=args.fdr,
        )
        emit_progress(0.90, "open_mass_search")
        emit_timing("open_mass_search", stage_started)
    else:
        emit_progress(0.90, "no_feature_rescue")
    evidence = (
        accepted
        + component_consensus
        + component_sequence_tags
        + consensus
        + targeted
        + targeted_glycopeptides
        + open_mass
    )
    evidence.sort(key=lambda row: (float(row.get("feature_ranking") or 0.0), float(row["score"])), reverse=True)
    annotations = differential_annotations(evidence)
    feature_evidence: list[dict[str, object]] = []
    feature_ms2_summary: dict[str, int] = {}
    if feature_payload is not None:
        stage_started = time.perf_counter()
        feature_evidence, feature_ms2_summary = build_feature_ms2_evidence(
            feature_payload,
            all_ms2_scans,
            evidence,
            annotations,
            candidates=candidates,
            sequence_region_candidates=sequence_region_candidates,
            fragment_tolerance_ppm=args.fragment_ppm,
        )
        emit_timing("feature_evidence", stage_started)
    emit_progress(0.93, "feature_evidence")
    peptide_form_comparisons: list[dict[str, object]] = []
    if feature_payload is not None and feature_sqlite is not None:
        stage_started = time.perf_counter()
        peptide_form_comparisons = build_peptide_form_comparisons(
            feature_payload,
            evidence,
            annotations,
            candidates,
            spectra_loader,
            mz_tolerance_ppm=max(0.1, args.precursor_ppm * 2.0),
        )
        emit_timing("peptide_form_comparisons", stage_started)
    modification_level_quantitation = build_modification_level_quantitation(
        peptide_form_comparisons
    )
    proteolytic_modification_families = build_proteolytic_modification_families(
        modification_level_quantitation
    )
    modified_peptide_findings = build_modified_peptide_findings(
        evidence,
        annotations,
        feature_payload,
        peptide_form_comparisons,
    )
    emit_progress(0.96, "modification_summary")
    if feature_payload is not None:
        stage_started = time.perf_counter()
        ms1_component_groups = build_ms1_component_groups(
            feature_payload,
            annotations,
            spectra_loader=spectra_loader,
        )
        emit_timing("final_component_grouping", stage_started)
    else:
        ms1_component_groups = []
    unidentified_reason_summary: dict[str, int] = {}
    for row in feature_evidence:
        reason = str(row.get("unidentified_reason") or "")
        if reason:
            unidentified_reason_summary[reason] = unidentified_reason_summary.get(reason, 0) + 1
    injection_volumes = {str(sample.get("injection_volume")) for sample in sample_summaries if sample.get("injection_volume") not in (None, "")}
    warnings = [
        "样品进样体积不同；原始峰面积比不能直接解释为产品差异，请结合上样浓度和局部 TIC 归一化复核。"
    ] if len(injection_volumes) > 1 else []
    shown = evidence[:max(0, args.max_report_psms)]
    report = {
        "schema_version": "LCMSMSIdentification/v2",
        "project_id": args.project_id,
        "fasta": str(Path(args.fasta).resolve()),
        "chains": chains,
        "chain_lengths": {name: len(sequence) for name, sequence in chains.items()},
        "target_candidates": sum(not candidate.is_decoy for candidate in candidates),
        "decoy_candidates": sum(candidate.is_decoy for candidate in candidates),
        "parameters": {
            "enzyme": args.enzyme,
            "enzyme_label": ENZYMES[args.enzyme].label,
            "sample_prep": args.sample_prep,
            "sample_prep_model": prep_model,
            "max_missed_cleavages": args.max_missed_cleavages,
            "max_variable_modifications": args.max_variable_modifications,
            "semitryptic_max_trim": args.semitryptic_max_trim,
            "precursor_ppm": args.precursor_ppm,
            "fragment_ppm": args.fragment_ppm,
            "fdr": args.fdr,
            "fixed_modification": prep_model["fixed_modification"],
            "variable_modifications": [
                "Oxidation/Dioxidation@M/W", "Deamidation@N/Q", "Succinimide@N",
                "PyroGlu@peptide N-terminus", "Glycation@K/N-terminus",
                "G0F/G1F/G2F@N-glycosylation sequon", "C-terminal Lys clipping",
            ],
            "variable_modification_limit": args.max_variable_modifications,
            "targeted_n_glycopeptide_search": {
                "enabled": True,
                "differential_features_with_selected_ms2_only": True,
                "requires_hexnac_204_and_additional_diagnostic_ion": True,
                "requires_core_y_and_peptide_backbone_evidence": True,
                "glycan_panel": [
                    "G0/G0F", "G1/G1F", "G2/G2F", "M5-M9",
                    "G1F/G2F with one or two NeuAc",
                ],
            },
            "sequence_inference_max_missed_cleavages": (
                args.sequence_inference_max_missed_cleavages
            ),
            "sequence_inference_max_terminal_trim": (
                args.sequence_inference_max_terminal_trim
            ),
            "sequence_inference_delta_range_da": [
                args.sequence_inference_min_delta_da,
                args.sequence_inference_max_delta_da,
            ],
            "ms2_workers": worker_count,
            "mzml_parser_cache": True,
        },
        "samples": sample_summaries,
        "warnings": warnings,
        "total_ms2_scans": total_ms2_scans,
        "total_accepted_psms": len(evidence),
        "standard_accepted_psms": len(accepted),
        "feature_guided_accepted_psms": len(targeted),
        "feature_glycopeptide_accepted_psms": len(targeted_glycopeptides),
        "feature_glycopeptide_feature_count": len({
            str(link.get("feature_group_id") or "")
            for psm in targeted_glycopeptides
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        }),
        "feature_glycopeptide_form_count": len({
            (
                str(psm.get("base_peptide_id") or ""),
                str(psm.get("glycan_name") or psm.get("modification_text") or ""),
            )
            for psm in targeted_glycopeptides
        }),
        "feature_open_mass_accepted_psms": len(open_mass),
        "consensus_feature_accepted_psms": len(consensus),
        "component_consensus_accepted_psms": len(component_consensus),
        "component_sequence_tag_accepted_psms": len(
            component_sequence_tags
        ),
        "sequence_region_candidate_count": len(
            sequence_region_candidates
        ),
        "sequence_inference_target_candidates": sum(
            not candidate.is_decoy
            for candidate in sequence_inference_candidates
        ),
        "provisional_ms1_component_count": sum(
            bool(row.get("inferred_component"))
            for row in provisional_component_groups
        ),
        "reported_psms": len(shown),
        "psms": shown,
        "annotations": annotations,
        "ms1_component_groups": ms1_component_groups,
        "feature_ms2_summary": feature_ms2_summary,
        "unidentified_reason_summary": unidentified_reason_summary,
        "feature_evidence": feature_evidence,
        "sequence_region_candidates": sequence_region_candidates,
        "modified_peptide_findings": modified_peptide_findings,
        "peptide_form_comparisons": peptide_form_comparisons,
        "modification_level_quantitation": modification_level_quantitation,
        "modification_level_quantitation_count": len(modification_level_quantitation),
        "proteolytic_modification_families": proteolytic_modification_families,
        "proteolytic_modification_family_count": len(proteolytic_modification_families),
        "high_value_proteolytic_consensus_count": sum(
            bool(row.get("high_value_consensus"))
            for row in proteolytic_modification_families
        ),
        "formal_modification_level_count": sum(
            row.get("quantitation_status") == "formal_relative_quantitation"
            for row in modification_level_quantitation
        ),
        "high_value_modification_event_count": sum(
            bool(row.get("high_value_candidate"))
            for row in modification_level_quantitation
        ),
    }
    if not args.disable_unimod_rescue and feature_payload is not None:
        emit_progress(0.965, "unimod_second_pass")
        stage_started = time.perf_counter()
        report["unimod_rescue"] = search_unimod_rescue(
            feature_evidence, all_ms2_scans, sequence_inference_candidates, chains,
            precursor_ppm=max(0.1, args.precursor_ppm * 2.0),
            fragment_ppm=max(0.1, args.fragment_ppm),
            min_delta=args.sequence_inference_min_delta_da,
            max_delta=args.sequence_inference_max_delta_da,
            sample_prep=report["parameters"]["sample_prep"],
            progress=lambda done, total: emit_progress(0.965 + 0.014 * done / max(1, total), "unimod_second_pass"),
        )
        emit_timing("unimod_second_pass", stage_started)
    else:
        report["unimod_rescue"] = {"enabled": False}
    json_path = output_dir / "lcms_msms_identifications.json"
    csv_path = output_dir / "lcms_msms_identifications.csv"
    feature_csv_path = output_dir / "lcms_ms1_ms2_feature_evidence.csv"
    modified_findings_csv_path = output_dir / "lcms_modified_peptide_findings.csv"
    modification_pairs_csv_path = output_dir / "lcms_modification_pairs.csv"
    modification_level_csv_path = output_dir / "lcms_modification_level_quantitation.csv"
    html_path = output_dir / "lcms_msms_report.html"
    emit_progress(0.98, "write_report")
    stage_started = time.perf_counter()
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csv(csv_path, evidence)
    write_feature_evidence_csv(feature_csv_path, feature_evidence)
    write_modified_findings_csv(modified_findings_csv_path, modified_peptide_findings)
    write_modification_pairs_csv(modification_pairs_csv_path, peptide_form_comparisons)
    write_modification_level_csv(modification_level_csv_path, modification_level_quantitation)
    html_path.write_text(html_report(report), encoding="utf-8")
    if feature_sqlite:
        write_sqlite_artifact(feature_sqlite, report)
    emit_timing("write_report", stage_started)
    emit_timing("total_ms2_workflow", overall_started)
    emit_progress(1.0, "ms2_complete")
    print(f"MS2 scans: {total_ms2_scans}")
    print(f"MS2 sample workers: {worker_count}")
    print(f"Accepted target PSMs: {len(accepted)}")
    print(f"Feature-guided tentative PSMs: {len(targeted)}")
    print(f"Feature-targeted glycopeptide PSMs: {len(targeted_glycopeptides)}")
    print(f"Selected-feature consensus tentative PSMs: {len(consensus)}")
    print(f"Component consensus tentative PSMs: {len(component_consensus)}")
    print(
        "Mass-offset/sequence-tag backbone PSMs: "
        f"{len(component_sequence_tags)}"
    )
    print(
        "Sequence-region candidates: "
        f"{len(sequence_region_candidates)}"
    )
    print(f"Differential annotations: {len(annotations)}")
    print(f"MS1 Feature evidence: {feature_ms2_summary}")
    print(f"Unidentified reasons: {unidentified_reason_summary}")
    print(f"Modified peptide findings (paired or unpaired): {len(modified_peptide_findings)}")
    print(f"Modified/unmodified peptide pairs: {len(peptide_form_comparisons)}")
    print(f"Modification-level quantitation families: {len(modification_level_quantitation)}")
    print(f"Report: {html_path}")


if __name__ == "__main__":
    main()
