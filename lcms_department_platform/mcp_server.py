"""MCP adapter for curated LC-MS analytical evidence and review annotations.

The adapter deliberately exposes curated analytical evidence instead of the
task directory or arbitrary SQLite access.  Apart from append-only agent
annotations that are always pending human review, it is read-only.  It is
transport-agnostic: the HTTP portal uses :class:`MCPApplication` for JSON-RPC
requests, while a future stdio entry point can reuse the same class.
"""

from __future__ import annotations

import json
import re
import sqlite3
from contextlib import closing
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable
from urllib.parse import unquote, urlparse
from uuid import uuid4

try:
    from serve_peak_first_compare import read_artifact, read_bootstrap, read_msms_identifications
except ModuleNotFoundError:  # Support package-level imports in tests and integrations.
    from lcms_feature_mvp.serve_peak_first_compare import read_artifact, read_bootstrap, read_msms_identifications


MCP_PROTOCOL_VERSION = "2025-06-18"
MCP_SERVER_NAME = "lcms-department-platform"
MCP_SERVER_VERSION = "1.0.0"
MAX_TOOL_LIMIT = 500
MAX_SPECTRUM_PEAKS = 500
MAX_UNKNOWN_ANALYSIS_LIMIT = 100
MAX_ANNOTATIONS_PER_FEATURE = 100
MAX_ANNOTATION_TEXT_LENGTH = 2_000
MAX_ANNOTATION_RATIONALE_LENGTH = 4_000
MAX_AGENT_NAME_LENGTH = 120
TASK_ID_RE = re.compile(r"^[^\\/]+$")
UNKNOWN_FEATURE_STATUSES = {
    "selected_precursor_unidentified",
    "coisolated_ms2_unresolved",
    "low_evidence_sequence_candidate",
}


def _json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _as_dict(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: object) -> list[Any]:
    return value if isinstance(value, list) else []


def _limit(value: object, default: int = 100) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError):
        number = default
    return max(1, min(number, MAX_TOOL_LIMIT))


def _unknown_analysis_limit(value: object, default: int = 25) -> int:
    return min(_limit(value, default), MAX_UNKNOWN_ANALYSIS_LIMIT)


def _trim_text(value: object, limit: int, field: str, *, required: bool = False) -> str:
    text = str(value or "").strip()
    if required and not text:
        raise ValueError(f"{field} is required")
    if len(text) > limit:
        raise ValueError(f"{field} exceeds the {limit} character limit")
    return text


def _ensure_agent_annotation_table(connection: sqlite3.Connection) -> None:
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS agent_feature_annotations (
            annotation_id TEXT PRIMARY KEY,
            feature_group_id TEXT NOT NULL,
            suggestion TEXT NOT NULL,
            rationale TEXT NOT NULL DEFAULT '',
            agent_name TEXT NOT NULL DEFAULT 'MCP agent',
            review_status TEXT NOT NULL DEFAULT 'pending_human_review',
            created_at TEXT NOT NULL
        )
        """
    )
    connection.execute(
        """
        CREATE INDEX IF NOT EXISTS idx_agent_feature_annotations_feature
        ON agent_feature_annotations (feature_group_id, created_at)
        """
    )


def read_agent_annotations(db_path: Path, feature_id: str | None = None) -> list[dict[str, Any]]:
    """Read append-only annotations without creating a table for legacy tasks."""

    try:
        with closing(sqlite3.connect(db_path)) as connection:
            if feature_id:
                rows = connection.execute(
                    """
                    SELECT annotation_id, feature_group_id, suggestion, rationale, agent_name, review_status, created_at
                    FROM agent_feature_annotations
                    WHERE feature_group_id = ?
                    ORDER BY created_at DESC, annotation_id DESC
                    """,
                    (feature_id,),
                ).fetchall()
            else:
                rows = connection.execute(
                    """
                    SELECT annotation_id, feature_group_id, suggestion, rationale, agent_name, review_status, created_at
                    FROM agent_feature_annotations
                    ORDER BY created_at DESC, annotation_id DESC
                    """
                ).fetchall()
    except sqlite3.OperationalError:
        return []
    fields = ("annotation_id", "feature_group_id", "suggestion", "rationale", "agent_name", "review_status", "created_at")
    return [dict(zip(fields, row, strict=True)) for row in rows]


def save_agent_annotation(
    db_path: Path,
    feature_id: str,
    suggestion: object,
    rationale: object = "",
    agent_name: object = "MCP agent",
) -> dict[str, Any]:
    """Append one review-only agent annotation; platform identifications stay untouched."""

    suggestion_text = _trim_text(suggestion, MAX_ANNOTATION_TEXT_LENGTH, "suggestion", required=True)
    rationale_text = _trim_text(rationale, MAX_ANNOTATION_RATIONALE_LENGTH, "rationale")
    agent_text = _trim_text(agent_name, MAX_AGENT_NAME_LENGTH, "agent_name") or "MCP agent"
    annotation_id = f"agent_{uuid4().hex}"
    created_at = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    with closing(sqlite3.connect(db_path)) as connection:
        _ensure_agent_annotation_table(connection)
        current_count = connection.execute(
            "SELECT COUNT(*) FROM agent_feature_annotations WHERE feature_group_id = ?",
            (feature_id,),
        ).fetchone()
        if int(current_count[0] if current_count else 0) >= MAX_ANNOTATIONS_PER_FEATURE:
            raise ValueError(f"feature already has the maximum of {MAX_ANNOTATIONS_PER_FEATURE} agent annotations")
        connection.execute(
            """
            INSERT INTO agent_feature_annotations
                (annotation_id, feature_group_id, suggestion, rationale, agent_name, review_status, created_at)
            VALUES (?, ?, ?, ?, ?, 'pending_human_review', ?)
            """,
            (annotation_id, feature_id, suggestion_text, rationale_text, agent_text, created_at),
        )
        connection.commit()
    return {
        "annotation_id": annotation_id,
        "feature_group_id": feature_id,
        "suggestion": suggestion_text,
        "rationale": rationale_text,
        "agent_name": agent_text,
        "review_status": "pending_human_review",
        "created_at": created_at,
    }


def _selected_task(task: dict[str, Any]) -> dict[str, Any]:
    """Return task metadata without exposing local file paths or references."""

    fields = (
        "task_id",
        "task_name",
        "status",
        "stage",
        "progress",
        "created_at",
        "updated_at",
        "error",
        "result_summary",
    )
    return {field: task.get(field) for field in fields if field in task}


def _selected_component(component: dict[str, Any]) -> dict[str, Any]:
    fields = (
        "component_group_id",
        "feature_group_id",
        "component_primary_feature_id",
        "component_label",
        "component_confidence",
        "component_neutral_mass",
        "component_representative_mz",
        "component_representative_charge",
        "component_charge_states",
        "component_isotope_offsets",
        "component_relation_types",
        "component_mass_evidence",
        "component_inference_score",
        "component_xic_shape_score",
        "identified_component",
        "inferred_component",
        "members",
        "area_by_sample",
        "normalized_area_by_sample",
        "presence_by_sample",
        "difference_type",
        "max_fold_change",
        "representative_rt",
    )
    result = {field: component.get(field) for field in fields if field in component}
    members = result.get("members")
    if isinstance(members, list):
        result["members"] = members[:100]
    return result


def _selected_evidence(evidence: dict[str, Any]) -> dict[str, Any]:
    fields = (
        "rank",
        "feature_group_id",
        "parent_tic_peak_id",
        "representative_rt",
        "representative_mz",
        "true_peak_mz",
        "envelope_representative_mz",
        "difference_type",
        "max_fold_change",
        "ranking_score",
        "higher_abundance_sample",
        "lower_abundance_sample",
        "normalized_area_by_sample",
        "ms2_status",
        "sequence",
        "modification",
        "confidence",
        "feature_link_type",
        "feature_isotope_offset",
        "selected_precursor_charge",
        "peptide_likelihood",
        "peptide_likelihood_reason",
        "unidentified_reason",
        "candidate_mass_hypotheses",
        "sequence_region_candidates",
        "coverage_by_sample",
        "best_psm",
        "candidate_psm",
        "coverage_scan",
    )
    return {field: evidence.get(field) for field in fields if field in evidence}


def _spectrum_payload(scan: object, max_peaks: int = MAX_SPECTRUM_PEAKS) -> dict[str, Any] | None:
    if not isinstance(scan, dict):
        return None
    peaks = _as_list(scan.get("spectrum_peaks"))
    peaks = [peak for peak in peaks if isinstance(peak, dict)]
    peaks = peaks[: max(1, min(int(max_peaks), MAX_SPECTRUM_PEAKS))]
    result = {
        key: scan.get(key)
        for key in (
            "sample_id",
            "scan_id",
            "rt",
            "precursor_mz",
            "precursor_charge",
            "activation_method",
            "collision_energy",
            "relation_type",
            "isotope_offset",
        )
        if key in scan
    }
    result["spectrum_peaks"] = peaks
    result["peak_count"] = len(peaks)
    result["truncated"] = len(_as_list(scan.get("spectrum_peaks"))) > len(peaks)
    return result


class MCPApplication:
    """Stateless MCP JSON-RPC application with one append-only review tool."""

    def __init__(self, config: Any, store: Any, audit: Callable[[str], None] | None = None) -> None:
        self.config = config
        self.store = store
        self.audit = audit or (lambda _message: None)

    def _task(self, task_id: str) -> dict[str, Any]:
        task_id = str(task_id or "")
        if not task_id or not TASK_ID_RE.fullmatch(task_id) or task_id in {".", ".."}:
            raise ValueError("invalid task_id")
        task = self.store.get(task_id)
        if not task:
            raise ValueError("task not found")
        task_root = (Path(self.config.jobs_dir) / task_id).resolve()
        jobs_root = Path(self.config.jobs_dir).resolve()
        if task_root.parent != jobs_root:
            raise ValueError("task path is outside the configured jobs directory")
        return task

    def _db(self, task_id: str) -> Path:
        self._task(task_id)
        path = (Path(self.config.jobs_dir) / task_id / "output" / "lcms_peak_first_compare.sqlite").resolve()
        jobs_root = Path(self.config.jobs_dir).resolve()
        if path.parent.parent.parent != jobs_root or not path.exists():
            raise ValueError("task result database is not available")
        return path

    @staticmethod
    def _component_groups(db_path: Path, msms: dict[str, Any] | None = None) -> list[dict[str, Any]]:
        if isinstance(msms, dict):
            groups = msms.get("ms1_component_groups")
            if isinstance(groups, list) and groups:
                return [item for item in groups if isinstance(item, dict)]
        for key in ("ms1_component_groups", "bootstrap"):
            try:
                payload = read_artifact(db_path, key)
            except KeyError:
                continue
            if key == "ms1_component_groups" and isinstance(payload, dict):
                groups = payload.get("component_groups")
            elif isinstance(payload, dict):
                groups = payload.get("ms1_component_groups")
            else:
                groups = None
            if isinstance(groups, list):
                return [item for item in groups if isinstance(item, dict)]
        return []

    @staticmethod
    def _feature_evidence(msms: dict[str, Any]) -> list[dict[str, Any]]:
        return [item for item in _as_list(msms.get("feature_evidence")) if isinstance(item, dict)]

    def _summary(self, task_id: str) -> dict[str, Any]:
        task = self._task(task_id)
        result: dict[str, Any] = {"task_id": task_id, "task": _selected_task(task)}
        if task.get("status") != "finished":
            result["analysis_available"] = False
            return result
        db_path = self._db(task_id)
        bootstrap = _as_dict(read_bootstrap(db_path))
        msms: dict[str, Any] = {}
        try:
            msms = _as_dict(read_msms_identifications(db_path))
        except KeyError:
            pass
        groups = self._component_groups(db_path, msms)
        evidence = self._feature_evidence(msms)
        status_counts: dict[str, int] = {}
        for item in evidence:
            status = str(item.get("ms2_status") or "unknown")
            status_counts[status] = status_counts.get(status, 0) + 1
        result.update(
            {
                "analysis_available": True,
                "project_id": bootstrap.get("project_id"),
                "sample_ids": bootstrap.get("sample_ids", []),
                "component_group_count": len(groups),
                "feature_evidence_count": len(evidence),
                "pending_agent_annotation_count": len(read_agent_annotations(db_path)),
                "ms2_status_counts": status_counts,
                "ms2_metrics": {
                    key: msms.get(key)
                    for key in (
                        "total_ms2_scans",
                        "total_accepted_psms",
                        "reported_psms",
                        "feature_guided_accepted_psms",
                        "feature_open_mass_accepted_psms",
                        "consensus_feature_accepted_psms",
                        "component_consensus_accepted_psms",
                        "component_sequence_tag_accepted_psms",
                        "feature_glycopeptide_accepted_psms",
                        "feature_glycopeptide_feature_count",
                        "feature_glycopeptide_form_count",
                        "modification_level_quantitation_count",
                        "formal_modification_level_count",
                        "high_value_modification_event_count",
                    )
                    if key in msms
                },
            }
        )
        return result

    def _list_features(self, arguments: dict[str, Any]) -> dict[str, Any]:
        task_id = str(arguments.get("task_id") or "")
        db_path = self._db(task_id)
        try:
            msms = _as_dict(read_msms_identifications(db_path))
        except KeyError:
            msms = {}
        evidence = self._feature_evidence(msms)
        include_no_ms2 = bool(arguments.get("include_no_ms2", False))
        status = str(arguments.get("ms2_status") or "").strip().lower()
        query = str(arguments.get("query") or "").strip().lower()
        filtered: list[dict[str, Any]] = []
        for item in evidence:
            item_status = str(item.get("ms2_status") or "unknown")
            if not include_no_ms2 and item_status.startswith("no_ms2"):
                continue
            if status and item_status.lower() != status:
                continue
            haystack = " ".join(
                str(item.get(field) or "")
                for field in ("feature_group_id", "sequence", "modification", "hypothesis", "unidentified_reason")
            ).lower()
            if query and query not in haystack:
                continue
            filtered.append(_selected_evidence(item))
        offset = max(0, int(arguments.get("offset") or 0))
        limit = _limit(arguments.get("limit"), 100)
        return {
            "task_id": task_id,
            "ms2_only_default": True,
            "returned_count": len(filtered[offset : offset + limit]),
            "total_matching": len(filtered),
            "offset": offset,
            "limit": limit,
            "features": filtered[offset : offset + limit],
        }

    def _feature_evidence_payload(self, arguments: dict[str, Any]) -> dict[str, Any]:
        task_id = str(arguments.get("task_id") or "")
        feature_id = str(arguments.get("feature_id") or "")
        if not feature_id:
            raise ValueError("feature_id is required")
        db_path = self._db(task_id)
        try:
            msms = _as_dict(read_msms_identifications(db_path))
        except KeyError:
            msms = {}
        evidence = next(
            (item for item in self._feature_evidence(msms) if str(item.get("feature_group_id")) == feature_id),
            None,
        )
        groups = self._component_groups(db_path, msms)
        component = next(
            (
                item
                for item in groups
                if str(item.get("feature_group_id")) == feature_id
                or str(item.get("component_primary_feature_id")) == feature_id
                or feature_id in {str(value) for value in _as_list(item.get("merged_feature_group_ids"))}
            ),
            None,
        )
        if evidence is None and component is None:
            raise ValueError("feature not found")
        coverage_scan = evidence.get("coverage_scan") if evidence else None
        ms2_available = isinstance(coverage_scan, dict) and bool(coverage_scan.get("spectrum_peaks"))
        result = {
            "task_id": task_id,
            "feature_group_id": feature_id,
            "ms2_available": ms2_available,
            "analysis_eligible": ms2_available,
            "analysis_policy": "no downstream analysis for features without stored MS2" if not ms2_available else "MS2 evidence available",
            "evidence": _selected_evidence(evidence) if evidence else None,
            "component": _selected_component(component) if component else None,
            "agent_annotations": read_agent_annotations(db_path, feature_id),
        }
        if arguments.get("include_spectrum", True) and ms2_available:
            result["spectrum"] = _spectrum_payload(coverage_scan, arguments.get("max_peaks", MAX_SPECTRUM_PEAKS))
        return result

    def _spectrum(self, arguments: dict[str, Any]) -> dict[str, Any]:
        arguments = dict(arguments)
        arguments["include_spectrum"] = True
        result = self._feature_evidence_payload(arguments)
        if not result.get("ms2_available"):
            return {
                "task_id": result["task_id"],
                "feature_group_id": result["feature_group_id"],
                "ms2_available": False,
                "spectrum": None,
                "message": "No displayable MS/MS spectrum is stored for this feature.",
            }
        return {
            "task_id": result["task_id"],
            "feature_group_id": result["feature_group_id"],
            "ms2_available": True,
            "spectrum": result.get("spectrum"),
        }

    @staticmethod
    def _analysis_scan(evidence: dict[str, Any]) -> dict[str, Any] | None:
        for key in ("best_psm", "candidate_psm", "coverage_scan"):
            candidate = evidence.get(key)
            if isinstance(candidate, dict) and _as_list(candidate.get("spectrum_peaks")):
                return candidate
        return None

    @staticmethod
    def _unknown_triage(evidence: dict[str, Any]) -> dict[str, str]:
        status = str(evidence.get("ms2_status") or "unknown")
        has_sequence = bool(str(evidence.get("sequence") or "").strip())
        if status == "low_evidence_sequence_candidate":
            return {
                "priority": "high",
                "next_step": "Review the competing sequence/modification candidate against the stored fragment spectrum; do not promote it without human confirmation.",
            }
        if status == "coisolated_ms2_unresolved":
            return {
                "priority": "high",
                "next_step": "The isolation window contains signal but attribution is unresolved; inspect co-isolation and consider a narrower-isolation reacquisition.",
            }
        if status == "selected_precursor_unidentified":
            return {
                "priority": "medium",
                "next_step": "Review the stored MS/MS spectrum, diagnostic ions and mass hypotheses; a new targeted acquisition may be needed if no coherent fragment series is present.",
            }
        return {
            "priority": "low" if has_sequence else "medium",
            "next_step": "Review the stored evidence before proposing a tentative annotation.",
        }

    def _analyze_unknown_features(self, arguments: dict[str, Any]) -> dict[str, Any]:
        """Return a bounded, evidence-only review queue; it does not rerun search algorithms."""

        task_id = str(arguments.get("task_id") or "")
        db_path = self._db(task_id)
        try:
            msms = _as_dict(read_msms_identifications(db_path))
        except KeyError:
            msms = {}
        requested_status = str(arguments.get("ms2_status") or "").strip().lower()
        include_no_ms2 = bool(arguments.get("include_no_ms2", False))
        include_spectrum = bool(arguments.get("include_spectrum", False))
        max_peaks = min(_limit(arguments.get("max_peaks"), 100), MAX_SPECTRUM_PEAKS)
        candidates: list[dict[str, Any]] = []
        for evidence in self._feature_evidence(msms):
            status = str(evidence.get("ms2_status") or "unknown")
            if requested_status:
                if status.lower() != requested_status:
                    continue
            elif status not in UNKNOWN_FEATURE_STATUSES:
                continue
            if not include_no_ms2 and status.startswith("no_ms2"):
                continue
            scan = self._analysis_scan(evidence)
            item = _selected_evidence(evidence)
            item["analysis_eligible"] = bool(scan)
            item["analysis_policy"] = "MS2 evidence available" if scan else "no downstream analysis for features without stored MS2"
            item["triage"] = self._unknown_triage(evidence)
            item["pending_agent_annotation_count"] = len(read_agent_annotations(db_path, str(evidence.get("feature_group_id") or "")))
            if include_spectrum and scan:
                item["spectrum"] = _spectrum_payload(scan, max_peaks)
            candidates.append(item)
        candidates.sort(
            key=lambda item: (
                {"high": 0, "medium": 1, "low": 2}.get(str(_as_dict(item.get("triage")).get("priority")), 3),
                -float(item.get("ranking_score") or 0.0),
                int(item.get("rank") or 10**9),
            )
        )
        offset = max(0, int(arguments.get("offset") or 0))
        limit = _unknown_analysis_limit(arguments.get("limit"), 25)
        return {
            "task_id": task_id,
            "analysis_kind": "bounded_unknown_feature_triage",
            "algorithm_execution": "none; this tool only curates existing MS1/MS2 evidence for agent review",
            "returned_count": len(candidates[offset : offset + limit]),
            "total_matching": len(candidates),
            "offset": offset,
            "limit": limit,
            "include_spectrum": include_spectrum,
            "features": candidates[offset : offset + limit],
        }

    def _save_agent_annotation(self, arguments: dict[str, Any]) -> dict[str, Any]:
        task_id = str(arguments.get("task_id") or "")
        feature_id = str(arguments.get("feature_id") or "").strip()
        if not feature_id:
            raise ValueError("feature_id is required")
        db_path = self._db(task_id)
        # Validate the feature before any database write.  An annotation can be
        # attached to a known component or a Feature evidence record only.
        self._feature_evidence_payload({"task_id": task_id, "feature_id": feature_id, "include_spectrum": False})
        annotation = save_agent_annotation(
            db_path,
            feature_id,
            arguments.get("suggestion"),
            arguments.get("rationale"),
            arguments.get("agent_name"),
        )
        return {
            "task_id": task_id,
            "annotation": annotation,
            "write_policy": "append_only_pending_human_review; platform identification fields were not changed",
        }

    def call_tool(self, name: str, arguments: dict[str, Any]) -> object:
        self.audit(f"tool={name} task_id={arguments.get('task_id', '')}")
        if name == "list_tasks":
            return {"tasks": [_selected_task(task) for task in self.store.load()]}
        if name == "get_task_status":
            task_id = str(arguments.get("task_id") or "")
            return {"task_id": task_id, "task": _selected_task(self._task(task_id))}
        if name == "get_task_summary":
            return self._summary(str(arguments.get("task_id") or ""))
        if name in {"list_features", "search_features"}:
            return self._list_features(arguments)
        if name == "get_feature_evidence":
            return self._feature_evidence_payload(arguments)
        if name == "get_msms_spectrum":
            return self._spectrum(arguments)
        if name == "analyze_unknown_features":
            return self._analyze_unknown_features(arguments)
        if name == "save_agent_annotation":
            return self._save_agent_annotation(arguments)
        raise ValueError(f"unknown tool: {name}")

    @staticmethod
    def tool_definitions() -> list[dict[str, Any]]:
        task = {"type": "string", "description": "平台任务 ID"}
        read_annotations = {"readOnlyHint": True, "destructiveHint": False, "idempotentHint": True, "openWorldHint": False}
        write_annotations = {"readOnlyHint": False, "destructiveHint": False, "idempotentHint": False, "openWorldHint": False}
        return [
            {
                "name": "list_tasks",
                "description": "列出任务及其状态。只返回任务元数据，不暴露本地文件路径。",
                "annotations": read_annotations,
                "inputSchema": {"type": "object", "properties": {}, "additionalProperties": False},
            },
            {
                "name": "get_task_status",
                "description": "读取一个任务的运行状态。",
                "annotations": read_annotations,
                "inputSchema": {"type": "object", "properties": {"task_id": task}, "required": ["task_id"]},
            },
            {
                "name": "get_task_summary",
                "description": "读取任务的样本、feature、MS2状态和鉴定统计摘要。",
                "annotations": read_annotations,
                "inputSchema": {"type": "object", "properties": {"task_id": task}, "required": ["task_id"]},
            },
            {
                "name": "list_features",
                "description": "检索任务中的 feature 证据；默认排除没有采集到 MS2 的 feature。",
                "annotations": read_annotations,
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "task_id": task,
                        "query": {"type": "string"},
                        "ms2_status": {"type": "string"},
                        "include_no_ms2": {"type": "boolean", "default": False},
                        "offset": {"type": "integer", "minimum": 0, "default": 0},
                        "limit": {"type": "integer", "minimum": 1, "maximum": MAX_TOOL_LIMIT, "default": 100},
                    },
                    "required": ["task_id"],
                },
            },
            {
                "name": "search_features",
                "description": "按序列、修饰、feature ID 或未定性原因搜索 MS2 feature。",
                "annotations": read_annotations,
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "task_id": task,
                        "query": {"type": "string"},
                        "ms2_status": {"type": "string"},
                        "include_no_ms2": {"type": "boolean", "default": False},
                        "offset": {"type": "integer", "minimum": 0, "default": 0},
                        "limit": {"type": "integer", "minimum": 1, "maximum": MAX_TOOL_LIMIT, "default": 100},
                    },
                    "required": ["task_id"],
                },
            },
            {
                "name": "get_feature_evidence",
                "description": "读取单个 feature 的序列、修饰、评分、MS2证据和成分信息。",
                "annotations": read_annotations,
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "task_id": task,
                        "feature_id": {"type": "string"},
                        "include_spectrum": {"type": "boolean", "default": True},
                        "max_peaks": {"type": "integer", "minimum": 1, "maximum": MAX_SPECTRUM_PEAKS, "default": MAX_SPECTRUM_PEAKS},
                    },
                    "required": ["task_id", "feature_id"],
                },
            },
            {
                "name": "get_msms_spectrum",
                "description": "读取单个 feature 保存的 MS/MS 谱图；没有 MS2 时返回空谱图和原因。",
                "annotations": read_annotations,
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "task_id": task,
                        "feature_id": {"type": "string"},
                        "max_peaks": {"type": "integer", "minimum": 1, "maximum": MAX_SPECTRUM_PEAKS, "default": MAX_SPECTRUM_PEAKS},
                    },
                    "required": ["task_id", "feature_id"],
                },
            },
            {
                "name": "analyze_unknown_features",
                "description": "批量整理已有 MS1/MS2 证据中的未知或低证据 Feature，提供受控的人工复核队列；不会重跑搜索算法，也不会修改平台结果。",
                "annotations": read_annotations,
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "task_id": task,
                        "ms2_status": {"type": "string", "description": "可选：只分析指定 MS2 状态"},
                        "include_no_ms2": {"type": "boolean", "default": False},
                        "include_spectrum": {"type": "boolean", "default": False},
                        "max_peaks": {"type": "integer", "minimum": 1, "maximum": MAX_SPECTRUM_PEAKS, "default": 100},
                        "offset": {"type": "integer", "minimum": 0, "default": 0},
                        "limit": {"type": "integer", "minimum": 1, "maximum": MAX_UNKNOWN_ANALYSIS_LIMIT, "default": 25},
                    },
                    "required": ["task_id"],
                },
            },
            {
                "name": "save_agent_annotation",
                "description": "向指定 Feature 追加 Agent 建议。每条记录固定为“待人工确认”，不会覆盖序列、修饰、置信等级或其他平台定性字段。",
                "annotations": write_annotations,
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "task_id": task,
                        "feature_id": {"type": "string"},
                        "suggestion": {"type": "string", "maxLength": MAX_ANNOTATION_TEXT_LENGTH},
                        "rationale": {"type": "string", "maxLength": MAX_ANNOTATION_RATIONALE_LENGTH},
                        "agent_name": {"type": "string", "maxLength": MAX_AGENT_NAME_LENGTH, "default": "MCP agent"},
                    },
                    "required": ["task_id", "feature_id", "suggestion"],
                },
            },
        ]

    def _resource_payload(self, uri: str) -> object:
        parsed = urlparse(uri)
        if parsed.scheme != "lcms" or parsed.netloc != "task":
            raise ValueError("unsupported resource URI")
        parts = [unquote(part) for part in parsed.path.strip("/").split("/") if part]
        if len(parts) < 2:
            raise ValueError("resource URI must be lcms://task/{task_id}/{resource}")
        task_id, resource = parts[0], parts[1]
        if resource == "summary":
            return self._summary(task_id)
        if resource == "features":
            return self._list_features({"task_id": task_id, "limit": MAX_TOOL_LIMIT, "include_no_ms2": False})
        if resource == "feature" and len(parts) >= 3:
            return self._feature_evidence_payload({"task_id": task_id, "feature_id": parts[2], "include_spectrum": True})
        if resource == "spectrum" and len(parts) >= 3:
            return self._spectrum({"task_id": task_id, "feature_id": parts[2]})
        raise ValueError("unsupported resource")

    def dispatch(self, request: dict[str, Any]) -> dict[str, Any] | None:
        request_id = request.get("id")
        method = str(request.get("method") or "")
        params = request.get("params")
        params = params if isinstance(params, dict) else {}
        if request_id is None and method.startswith("notifications/"):
            return None
        try:
            if method == "initialize":
                requested = str(params.get("protocolVersion") or MCP_PROTOCOL_VERSION)
                protocol = requested if requested == MCP_PROTOCOL_VERSION else MCP_PROTOCOL_VERSION
                result = {
                    "protocolVersion": protocol,
                    "capabilities": {
                        "tools": {},
                        "resources": {"subscribe": False, "listChanged": False},
                    },
                    "serverInfo": {"name": MCP_SERVER_NAME, "version": MCP_SERVER_VERSION},
                    "instructions": "LC-MS analysis server. Features without stored MS2 are excluded from downstream analysis by default. The only write operation is append-only Agent annotations, which are always pending human review and never alter platform identifications.",
                }
            elif method == "ping":
                result = {}
            elif method == "tools/list":
                result = {"tools": self.tool_definitions()}
            elif method == "tools/call":
                name = str(params.get("name") or "")
                arguments = params.get("arguments")
                arguments = arguments if isinstance(arguments, dict) else {}
                try:
                    payload = self.call_tool(name, arguments)
                    result = {"content": [{"type": "text", "text": _json(payload)}], "isError": False}
                except Exception as exc:
                    result = {"content": [{"type": "text", "text": _json({"error": str(exc)})}], "isError": True}
            elif method == "resources/list":
                result = {"resources": []}
            elif method == "resources/templates/list":
                result = {
                    "resourceTemplates": [
                        {"uriTemplate": "lcms://task/{task_id}/summary", "name": "Task summary", "mimeType": "application/json"},
                        {"uriTemplate": "lcms://task/{task_id}/features", "name": "MS2 feature evidence", "mimeType": "application/json"},
                        {"uriTemplate": "lcms://task/{task_id}/feature/{feature_id}", "name": "Feature evidence", "mimeType": "application/json"},
                        {"uriTemplate": "lcms://task/{task_id}/spectrum/{feature_id}", "name": "Feature MS/MS spectrum", "mimeType": "application/json"},
                    ]
                }
            elif method == "resources/read":
                uri = str(params.get("uri") or "")
                payload = self._resource_payload(uri)
                result = {"contents": [{"uri": uri, "mimeType": "application/json", "text": _json(payload)}]}
            elif method == "prompts/list":
                result = {"prompts": []}
            else:
                return {"jsonrpc": "2.0", "id": request_id, "error": {"code": -32601, "message": f"method not found: {method}"}}
            return {"jsonrpc": "2.0", "id": request_id, "result": result}
        except Exception as exc:
            return {"jsonrpc": "2.0", "id": request_id, "error": {"code": -32602, "message": str(exc)}}

    def handle_json(self, body: bytes) -> dict[str, Any] | list[dict[str, Any]] | None:
        payload = json.loads(body.decode("utf-8"))
        if isinstance(payload, list):
            responses = [self.dispatch(item) for item in payload if isinstance(item, dict)]
            return [item for item in responses if item is not None]
        if not isinstance(payload, dict):
            raise ValueError("MCP request must be a JSON-RPC object or batch")
        return self.dispatch(payload)
