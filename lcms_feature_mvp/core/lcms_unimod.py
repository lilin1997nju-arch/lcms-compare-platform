"""Offline Unimod rules and bounded, non-identifying second-pass MS2 search.

This module never adds its candidates to formal PSMs, FDR or quantitation.
Only one new modification on an otherwise unmodified (optionally fixed-CAM)
backbone is considered. Neutral-loss rules are retained but not scored.
"""
from bisect import bisect_left, bisect_right
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import time
import xml.etree.ElementTree as ET

from .lcms_msms import (
    PROTON, _candidate, _prepared_observed_spectrum,
    _screen_mass_offset_fragments, _labeled_spectrum_peaks,
    match_fragments, theoretical_fragments,
)
from .lcms_sample_prep import PREP_FILTER_VERSION, normalize_sample_prep, filter_reagent_rules, sample_prep_model, filter_linear_candidates

DATA = Path(__file__).resolve().parent / "data"
NS = {"u": "http://www.unimod.org/xmlns/schema/unimod_2"}
# These require different experimental/search models, not a single PTM.
EXCLUDED_CLASSES = {"AA substitution", "Isotopic label", "Cross-link"}
NOTE = ("Unimod 定向搜索，仅供参考：未经独立 FDR 校准；质量与 b/y 匹配不等于"
        "序列、修饰身份或位点已确认，不参与正式鉴定和修饰定量。")


class UnimodCatalog:
    def __init__(self, xml_path=DATA / "unimod.xml"):
        data = Path(xml_path).read_bytes()
        if b"<!DOCTYPE" in data.upper():
            raise ValueError("Unexpected DTD in Unimod snapshot")
        root = ET.fromstring(data)
        self.records = []
        for mod in root.findall("u:modifications/u:mod", NS):
            delta = mod.find("u:delta", NS)
            if delta is None:
                continue
            rules = []
            for rule in mod.findall("u:specificity", NS):
                rules.append({**rule.attrib, "neutral_losses": [
                    dict(loss.attrib) for loss in rule.findall("u:NeutralLoss", NS)
                ], "notes": rule.findtext("u:misc_notes", default="", namespaces=NS)})
            self.records.append({
                "id": int(mod.attrib["record_id"]), "title": mod.attrib["title"],
                "full_name": mod.attrib.get("full_name", ""),
                "mass": float(delta.attrib["mono_mass"]),
                "composition": delta.attrib.get("composition", ""),
                "approved": mod.attrib.get("approved") == "1", "rules": rules,
            })
        if not self.records:
            raise ValueError("Unimod snapshot contains no modifications")
        self.records.sort(key=lambda r: (r["mass"], r["id"]))
        self.masses = [r["mass"] for r in self.records]
        self.metadata = {"record_count": len(self.records),
                         "sha256": hashlib.sha256(data).hexdigest(),
                         "source": "https://www.unimod.org/xml/unimod.xml"}
        manifest = Path(xml_path).with_name("unimod-snapshot.json")
        if manifest.exists():
            info = json.loads(manifest.read_text(encoding="utf-8"))
            if info["sha256"] != self.metadata["sha256"]:
                raise ValueError("Unimod snapshot checksum mismatch")
            self.metadata.update(info)

    def match_mass(self, delta, tolerance):
        if not math.isfinite(delta) or not math.isfinite(tolerance) or tolerance < 0:
            return []
        return self.records[bisect_left(self.masses, delta-tolerance):
                            bisect_right(self.masses, delta+tolerance)]


@lru_cache(maxsize=1)
def load_catalog():
    return UnimodCatalog()


def applicable_sites(record, candidate, chain_length=None):
    """Return zero-based mass-placement indices, preserving terminal semantics.

    The fragmentation engine attaches terminal delta to the first/last residue;
    the original chemical site and position remain in each rule for display.
    Protein termini are accepted only with known FASTA coordinates/length.
    """
    result = []
    sequence = candidate.sequence
    # Monolinks may be classified as Chemical derivative, and virtual ion
    # entries as Other. Neither is a generic intact-peptide PTM model.
    title = str(record.get("title") or "")
    if title.startswith("Xlink:") or title.endswith("-type-ion"):
        return result
    for rule in record["rules"]:
        classification = rule.get("classification", "")
        if classification in EXCLUDED_CLASSES or "cross-link" in classification.lower():
            continue
        # Internal amino-acid deletion needs a different backbone, not a
        # residue with zero mass. Retain real protein-terminal clipping rules.
        if title.endswith("-loss") and rule.get("position") == "Anywhere":
            continue
        site, position = rule.get("site"), rule.get("position")
        if site == "N-term":
            indices = [0]
        elif site == "C-term":
            indices = [len(sequence)-1]
        elif isinstance(site, str) and len(site) == 1 and site in "ACDEFGHIKLMNPQRSTVWY":
            indices = [i for i, aa in enumerate(sequence) if aa == site]
        else:
            continue
        for index in indices:
            if position in {"Any N-term", "Protein N-term"} and index != 0:
                continue
            if position in {"Any C-term", "Protein C-term"} and index != len(sequence)-1:
                continue
            if position == "Protein N-term" and candidate.start != 1:
                continue
            if position == "Protein C-term" and (chain_length is None or candidate.end != chain_length):
                continue
            if position not in {"Anywhere", "Any N-term", "Any C-term", "Protein N-term", "Protein C-term"}:
                continue
            # Unimod deltas are relative to unmodified residues, not CAM-Cys.
            if site == "C" and candidate.carbamidomethyl_cys:
                continue
            if any(i == index for i, _, _ in candidate.modifications):
                continue
            result.append((index, rule))
    return result


def search_unimod_rescue(feature_evidence, scans, backbones, chains, *,
                         precursor_ppm=20.0, fragment_ppm=20.0,
                         min_delta=-250.0, max_delta=2500.0,
                         max_backbones_per_scan=12, max_results_per_sample=3,
                         catalog=None, progress=None, sample_prep=None):
    """Attach a separate candidate list; leave all first-pass evidence intact.

    Reuse the first-pass selected scan per sample (not isolation-only scans).
    Screen mass/site-compatible backbones, then score every allowed position
    for the top twelve backbones. No candidate site is silently declared unique.
    Limits and discarded backbone counts are recorded in the search summary.
    """
    started = time.perf_counter()
    catalog = catalog or load_catalog()
    sample_prep = normalize_sample_prep(sample_prep)
    permitted = {}
    excluded_ids = []
    removed_rules = 0
    for record in catalog.records:
        rules = filter_reagent_rules(record, sample_prep)
        removed_rules += len(record["rules"]) - len(rules)
        if rules:
            permitted[record["id"]] = {**record, "rules": rules}
        else:
            excluded_ids.append(record["id"])
    targets = sorted((c for c in filter_linear_candidates(backbones, sample_prep) if not c.is_decoy and not c.modifications),
                     key=lambda c: (c.neutral_mass, c.candidate_id))
    masses = [c.neutral_mass for c in targets]
    scan_index = {(s.sample_id, s.scan_id): s for s in scans}
    summary = {"enabled": True, "catalog": catalog.metadata, "note": NOTE,
               "sample_prep": sample_prep, "prep_filter_version": PREP_FILTER_VERSION,
               "sample_prep_model": sample_prep_model(sample_prep),
               "prep_excluded_record_ids": sorted(excluded_ids),
               "prep_excluded_rule_count": removed_rules,
               "max_new_modifications": 1, "max_backbones_per_scan": max_backbones_per_scan,
               "max_results_per_sample": max_results_per_sample,
               "precursor_ppm": precursor_ppm, "fragment_ppm": fragment_ppm,
               "delta_range_da": [min_delta, max_delta],
               "excluded_classes": sorted(EXCLUDED_CLASSES),
               "excluded_models": ["Xlink:*", "*-type-ion", "internal *-loss"],
               "neutral_loss_scoring": False, "independent_fdr": False,
               "eligible_features": 0, "selected_scans": 0,
               "mass_site_backbones": 0, "screened_backbones": 0,
               "discarded_backbones": 0, "scored_site_candidates": 0,
               "features_with_candidates": 0, "reported_candidates": 0}
    for feature_index, feature in enumerate(feature_evidence):
        if progress and feature_index % 25 == 0:
            progress(feature_index, len(feature_evidence))
        # Existing B/C-level interpretations are not re-searched or overwritten.
        confidence = str(feature.get("confidence") or "")
        if confidence.startswith(("B_", "C_")):
            continue
        if feature.get("difference_type") not in {"presence_absence", "area_changed", "moderate_difference"}:
            continue
        feature["unimod_rescue_candidates"] = []
        feature["unimod_rescue_note"] = NOTE
        coverages = feature.get("coverage_scan_by_sample") or {}
        if not coverages and feature.get("coverage_scan"):
            coverages = {feature["coverage_scan"]["sample_id"]: feature["coverage_scan"]}
        selected = [r for r in coverages.values() if r and r.get("relation_type") in
                    {"selected_precursor", "selected_isotope_envelope"}]
        if not selected:
            continue
        summary["eligible_features"] += 1
        for coverage in selected:
            # The shared legacy fragment cache is unbounded. Keep the new
            # exploratory site enumeration from retaining every trial peptide.
            if theoretical_fragments.cache_info().currsize > 16000:
                theoretical_fragments.cache_clear()
            scan = scan_index.get((coverage["sample_id"], coverage["scan_id"]))
            charge = int(coverage.get("precursor_charge") or 0)
            if scan is None or charge <= 0 or not scan.precursor_mz:
                continue
            summary["selected_scans"] += 1
            # Same selected precursor is used for both mass filtering and scoring.
            # The coverage offset describes feature vs selected isotope; it must
            # not be subtracted from the selected scan mass a second time.
            neutral_mass = (scan.precursor_mz - PROTON) * charge
            tolerance = neutral_mass * precursor_ppm / 1e6
            mz, intensities, _ = _prepared_observed_spectrum(scan)
            screened = []
            for backbone in targets[bisect_left(masses, neutral_mass-max_delta-tolerance):
                                    bisect_right(masses, neutral_mass-min_delta+tolerance)]:
                delta = neutral_mass - backbone.neutral_mass
                if abs(delta) < 0.5:
                    continue
                options = []
                for record in catalog.match_mass(delta, tolerance):
                    record = permitted.get(record["id"])
                    if record is None:
                        continue
                    sites = applicable_sites(record, backbone, len(chains[backbone.chain]) if backbone.chain in chains else None)
                    if sites:
                        options.append((record, sites))
                if not options:
                    continue
                summary["mass_site_backbones"] += 1
                screen = _screen_mass_offset_fragments(backbone, delta, mz, intensities, fragment_ppm)
                matches = screen.get("matched_fragments") or []
                if len(matches) < 4:
                    continue
                rank = (len(matches), float(screen.get("explained_intensity") or 0))
                screened.append((rank, backbone, options))
            screened.sort(key=lambda item: (item[0], item[1].candidate_id), reverse=True)
            summary["screened_backbones"] += len(screened)
            summary["discarded_backbones"] += max(0, len(screened)-max_backbones_per_scan)
            results = []
            for _, backbone, options in screened[:max_backbones_per_scan]:
                for record, sites in options:
                    scores_by_index = {}
                    for index, rule in sites:
                        if index in scores_by_index:
                            continue
                        candidate = _candidate(backbone.chain, backbone.start, backbone.end, backbone.sequence,
                                               ((index, f"UNIMOD:{record['id']}:{record['title']}", record["mass"]),),
                                               carbamidomethyl_cys=backbone.carbamidomethyl_cys,
                                               proteolysis=backbone.proteolysis)
                        error = (neutral_mass-candidate.neutral_mass)/candidate.neutral_mass*1e6
                        if abs(error) > precursor_ppm:
                            continue
                        evidence = match_fragments(scan, candidate, error, precursor_ppm, fragment_ppm,
                                                   include_extended_fragments=False)
                        summary["scored_site_candidates"] += 1
                        scores_by_index[index] = (candidate, error, evidence)
                    if not scores_by_index:
                        continue
                    index, (candidate, error, evidence) = max(scores_by_index.items(), key=lambda item: item[1][2]["score"])
                    # Reporting gate only. This is NOT an identification/FDR cutoff.
                    if evidence["matched_ion_count"] < 4 or evidence["fragment_coverage"] < 0.10:
                        continue
                    fragments = evidence["matched_fragments"]
                    alternatives = [{"position": i+1, "site": rule["site"],
                                     "specificity": rule["position"], "classification": rule.get("classification"),
                                     "hidden": rule.get("hidden") == "1",
                                     "neutral_losses": rule.get("neutral_losses", []),
                                     "score": scores_by_index[i][2]["score"]}
                                    for i, rule in sites if i in scores_by_index]
                    results.append({
                        "feature_group_id": feature["feature_group_id"],
                        "sample_id": scan.sample_id, "scan_id": scan.scan_id, "rt": scan.rt,
                        "precursor_mz": scan.precursor_mz, "precursor_charge": charge,
                        "precursor_error_ppm": error, "chain": backbone.chain,
                        "start": backbone.start, "end": backbone.end, "sequence": backbone.sequence,
                        "proteolysis": backbone.proteolysis, "modification_text": candidate.modification_text,
                        "mass_delta": record["mass"], "unimod_id": record["id"],
                        "unimod_title": record["title"], "unimod_approved": record["approved"],
                        "unimod_composition": record["composition"],
                        "display_site": index+1, "site_candidates": alternatives,
                        "matched_by_ion_count": evidence["matched_ion_count"],
                        "localization_status": "not_validated", "q_value": None,
                        "exploratory_only": True, "exploratory_note": NOTE,
                        "search_origin": "unimod_second_pass",
                        **{k: v for k, v in evidence.items() if k != "matched_fragments"},
                        "spectrum_peaks": _labeled_spectrum_peaks(scan, fragments),
                    })
            results.sort(key=lambda r: (-r["score"], abs(r["precursor_error_ppm"]), r["sequence"], r["unimod_id"]))
            feature["unimod_rescue_candidates"].extend(results[:max_results_per_sample])
        hits = feature["unimod_rescue_candidates"]
        if hits:
            hits.sort(key=lambda r: (-r["score"], r["sample_id"], r["unimod_id"]))
            summary["features_with_candidates"] += 1
            summary["reported_candidates"] += len(hits)
    summary["elapsed_seconds"] = time.perf_counter()-started
    return summary
