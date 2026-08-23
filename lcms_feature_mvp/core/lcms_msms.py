#!/usr/bin/env python3
"""Sequence-driven peptide MS/MS search using only the Python standard library."""

from __future__ import annotations

import json
import math
import re
import sqlite3
import statistics
from bisect import bisect_left, bisect_right
from collections import Counter
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Callable

from .lcms_models import LCMSSpectrumScan, XICTrace
from .lcms_peak_detection import detect_xic_peaks


PROTON = 1.007276466621
WATER = 18.010564684
CARBAMIDOMETHYL = 57.021463735
ISOTOPE_MASS_DIFF = 1.00335483507
AVERAGINE_UNIT_MASS = 111.1254


MODIFICATION_SITE_PATTERN = re.compile(r"^(?P<name>.+)@(?P<position>\d+)$")
DISPLAY_MODIFICATION_SITE_PATTERN = re.compile(
    r"^(?P<name>.+)@(?P<residue>[A-Z]?)(?P<position>\d+)$"
)


def modification_token_with_residue(token: str, sequence: str) -> str:
    """Add the peptide residue letter to a localized modification label."""
    text = str(token or "").strip()
    match = MODIFICATION_SITE_PATTERN.match(text)
    if not match:
        return text
    name = match.group("name").strip()
    position = int(match.group("position"))
    if name.lower().startswith(("n-terminal", "c-terminal")):
        return f"{name}@{position}"
    residue = sequence[position - 1] if 1 <= position <= len(sequence) else ""
    if not residue or name.endswith(f"-{residue}"):
        return f"{name}@{position}"
    return f"{name}@{residue}{position}"


def format_modification_alternatives(alternatives: list[str], sequence: str) -> str:
    """Render isobaric localization alternatives without implying coexistence."""
    clean = [str(item or "").strip() for item in alternatives if str(item or "").strip()]
    if not clean:
        return ""
    raw_token_lists = [
        [token.strip() for token in item.split(";") if token.strip()]
        for item in clean
    ]
    if len(raw_token_lists) == 1:
        return "；".join(raw_token_lists[0])
    token_lists = [
        [modification_token_with_residue(token, sequence) for token in item.split(";") if token.strip()]
        for item in clean
    ]

    common = Counter(token_lists[0])
    for tokens in token_lists[1:]:
        common &= Counter(tokens)
    common_tokens: list[str] = []
    remaining_common = Counter(common)
    for token in token_lists[0]:
        if remaining_common[token] > 0:
            common_tokens.append(token)
            remaining_common[token] -= 1

    residual_lists: list[list[str]] = []
    for tokens in token_lists:
        residual_common = Counter(common)
        residual: list[str] = []
        for token in tokens:
            if residual_common[token] > 0:
                residual_common[token] -= 1
            else:
                residual.append(token)
        residual_lists.append(residual)

    if all(len(tokens) == 1 for tokens in residual_lists):
        parsed_sites = [DISPLAY_MODIFICATION_SITE_PATTERN.match(tokens[0]) for tokens in residual_lists]
        if all(parsed_sites):
            names = [match.group("name") for match in parsed_sites if match]
            if len(set(names)) == 1:
                site_labels = []
                for match in parsed_sites:
                    assert match is not None
                    position = int(match.group("position"))
                    residue = match.group("residue") or (
                        sequence[position - 1] if 1 <= position <= len(sequence) else ""
                    )
                    site_labels.append(f"{residue}{position}" if residue else str(position))
                site_labels = sorted(set(site_labels), key=lambda value: int(re.search(r"\d+", value).group()))
                ambiguous = f"{names[0]}@{'/'.join(site_labels)}（位点未区分）"
                return "；".join([*common_tokens, ambiguous])

    option_labels = []
    for index, tokens in enumerate(token_lists):
        option_name = chr(ord("A") + index) if index < 26 else str(index + 1)
        option_labels.append(f"方案{option_name}：{'；'.join(tokens)}")
    return "｜".join(option_labels) + "（方案未区分）"
AVERAGINE_COMPOSITION = {
    "C": 4.9384,
    "H": 7.7583,
    "N": 1.3577,
    "O": 1.4773,
    "S": 0.0417,
}
# Nominal neutron shifts and natural abundances.  The truncated polynomial is
# sufficient for scoring the first isotopes normally visible in peptide MS1.
NATURAL_ISOTOPE_POLYNOMIALS = {
    "C": ((0, 0.9893), (1, 0.0107)),
    "H": ((0, 0.999885), (1, 0.000115)),
    "N": ((0, 0.99636), (1, 0.00364)),
    "O": ((0, 0.99757), (1, 0.00038), (2, 0.00205)),
    "S": ((0, 0.9499), (1, 0.0075), (2, 0.0425), (4, 0.0001)),
}
AA_MASS = {
    "A": 71.037113805, "R": 156.10111105, "N": 114.04292747, "D": 115.026943065,
    "C": 103.009184505, "E": 129.042593135, "Q": 128.05857754, "G": 57.021463735,
    "H": 137.058911875, "I": 113.084063975, "L": 113.084063975, "K": 128.094963015,
    "M": 131.040484645, "F": 147.068413945, "P": 97.052763875, "S": 87.032028435,
    "T": 101.047678505, "W": 186.07931298, "Y": 163.063328575, "V": 99.068413945,
}
COMMON_MODIFICATIONS = (
    ("Oxidation", frozenset("MW"), 15.99491462, "residue"),
    ("Dioxidation", frozenset("MW"), 31.98982924, "residue"),
    ("Deamidation", frozenset("NQ"), 0.984015585, "residue"),
    ("Succinimide", frozenset("N"), -17.026549105, "residue"),
    ("PyroGlu-Q", frozenset("Q"), -17.026549105, "n_term"),
    ("PyroGlu-E", frozenset("E"), -18.010564684, "n_term"),
    ("Glycation", frozenset("K"), 162.05282343, "residue"),
    ("N-terminal glycation", frozenset(AA_MASS), 162.05282343, "n_term"),
    ("G0F N-glycan", frozenset("N"), 1444.533870, "n_glycan_sequon"),
    ("G1F N-glycan", frozenset("N"), 1606.586694, "n_glycan_sequon"),
    ("G2F N-glycan", frozenset("N"), 1768.639517, "n_glycan_sequon"),
)

GLYCAN_RESIDUE_MASS = {
    "Hex": 162.0528234315,
    "HexNAc": 203.0793725330,
    "Fuc": 146.0579088090,
    "NeuAc": 291.0954165270,
}


def _glycan_composition_mass(composition: dict[str, int]) -> float:
    return sum(
        GLYCAN_RESIDUE_MASS[name] * max(0, int(count))
        for name, count in composition.items()
    )


# A deliberately bounded panel for therapeutic-antibody peptide mapping.  It
# is queried only for differential Features with a selected MS2 precursor and
# glycan-diagnostic ions; this is not a dataset-wide glycoproteomics search.
TARGETED_N_GLYCANS = tuple(
    {
        "name": name,
        "composition": composition,
        "mass": _glycan_composition_mass(composition),
    }
    for name, composition in (
        ("G0 N-glycan", {"Hex": 3, "HexNAc": 4}),
        ("G0F N-glycan", {"Hex": 3, "HexNAc": 4, "Fuc": 1}),
        ("G1 N-glycan", {"Hex": 4, "HexNAc": 4}),
        ("G1F N-glycan", {"Hex": 4, "HexNAc": 4, "Fuc": 1}),
        ("G2 N-glycan", {"Hex": 5, "HexNAc": 4}),
        ("G2F N-glycan", {"Hex": 5, "HexNAc": 4, "Fuc": 1}),
        ("M5 N-glycan", {"Hex": 5, "HexNAc": 2}),
        ("M6 N-glycan", {"Hex": 6, "HexNAc": 2}),
        ("M7 N-glycan", {"Hex": 7, "HexNAc": 2}),
        ("M8 N-glycan", {"Hex": 8, "HexNAc": 2}),
        ("M9 N-glycan", {"Hex": 9, "HexNAc": 2}),
        ("G1F+NeuAc1 N-glycan", {"Hex": 4, "HexNAc": 4, "Fuc": 1, "NeuAc": 1}),
        ("G2F+NeuAc1 N-glycan", {"Hex": 5, "HexNAc": 4, "Fuc": 1, "NeuAc": 1}),
        ("G2F+NeuAc2 N-glycan", {"Hex": 5, "HexNAc": 4, "Fuc": 1, "NeuAc": 2}),
    )
)


def _convolve_isotope_polynomials(
    first: tuple[float, ...],
    second: tuple[float, ...],
    isotope_count: int,
) -> tuple[float, ...]:
    output = [0.0] * isotope_count
    for first_index, first_value in enumerate(first):
        if first_value <= 0.0:
            continue
        for second_index, second_value in enumerate(second[: isotope_count - first_index]):
            if second_value > 0.0:
                output[first_index + second_index] += first_value * second_value
    return tuple(output)


def _isotope_polynomial_power(
    base: tuple[float, ...],
    exponent: int,
    isotope_count: int,
) -> tuple[float, ...]:
    result = (1.0,) + (0.0,) * (isotope_count - 1)
    factor = base
    remaining = max(0, int(exponent))
    while remaining:
        if remaining & 1:
            result = _convolve_isotope_polynomials(result, factor, isotope_count)
        remaining >>= 1
        if remaining:
            factor = _convolve_isotope_polynomials(factor, factor, isotope_count)
    return result


@lru_cache(maxsize=4096)
def _cached_averagine_isotope_distribution(
    rounded_neutral_mass: int,
    isotope_count: int,
) -> tuple[float, ...]:
    count = max(3, min(16, int(isotope_count)))
    scale = max(1.0, float(rounded_neutral_mass)) / AVERAGINE_UNIT_MASS
    distribution = (1.0,) + (0.0,) * (count - 1)
    for element, average_count in AVERAGINE_COMPOSITION.items():
        atom_count = max(0, int(round(average_count * scale)))
        base = [0.0] * count
        for shift, abundance in NATURAL_ISOTOPE_POLYNOMIALS[element]:
            if shift < count:
                base[shift] = abundance
        element_distribution = _isotope_polynomial_power(tuple(base), atom_count, count)
        distribution = _convolve_isotope_polynomials(distribution, element_distribution, count)
    total = sum(distribution)
    if total <= 0.0:
        return (1.0,) + (0.0,) * (count - 1)
    return tuple(value / total for value in distribution)


def theoretical_averagine_isotope_distribution(
    neutral_mass: float,
    isotope_count: int = 10,
) -> list[float]:
    """Return a normalized nominal-isotope envelope for an average peptide."""
    return list(_cached_averagine_isotope_distribution(int(round(max(1.0, neutral_mass))), isotope_count))


def fit_theoretical_isotope_envelope(
    observed_intensities: list[float],
    neutral_mass_observed_first: float,
    max_monoisotopic_offset: int = 4,
) -> dict[str, object]:
    """Fit an observed envelope to shifted averagine isotope distributions.

    The shift is the number of theoretical isotopes preceding the first
    observed peak.  This directly detects cases where M+1 or M+2 was picked as
    the first observable isotope and supplies a confidence score for the full
    relative-intensity pattern rather than isotope spacing alone.
    """
    observed = [max(0.0, float(value)) for value in observed_intensities]
    observed_norm = math.sqrt(sum(value * value for value in observed))
    if not observed or observed_norm <= 0.0:
        return {
            "monoisotopic_offset": 0,
            "cosine_score": 0.0,
            "explained_theoretical_fraction": 0.0,
            "fit_score": 0.0,
            "theoretical_intensity": [],
        }
    maximum_offset = max(0, min(6, int(max_monoisotopic_offset)))
    isotope_count = min(16, len(observed) + maximum_offset + 4)
    best: tuple[float, int, dict[str, object]] | None = None
    for offset in range(maximum_offset + 1):
        neutral_mass = max(
            1.0,
            float(neutral_mass_observed_first) - offset * ISOTOPE_MASS_DIFF,
        )
        theoretical = theoretical_averagine_isotope_distribution(neutral_mass, isotope_count)
        expected = theoretical[offset: offset + len(observed)]
        if len(expected) != len(observed):
            continue
        expected_norm = math.sqrt(sum(value * value for value in expected))
        cosine = (
            sum(first * second for first, second in zip(observed, expected))
            / max(observed_norm * expected_norm, 1e-12)
        )
        explained_fraction = sum(expected) / max(sum(theoretical), 1e-12)
        # Shape is the primary evidence for a missing M/M+1 peak.  Coverage is
        # deliberately a weak term because the true monoisotopic peak may be
        # below the detection threshold for larger peptides.
        fit_score = max(0.0, min(1.0, cosine)) * (
            0.96 + 0.04 * math.sqrt(max(0.0, min(1.0, explained_fraction)))
        ) - 0.01 * offset
        fit_score = max(0.0, min(1.0, fit_score))
        result = {
            "monoisotopic_offset": offset,
            "cosine_score": max(0.0, min(1.0, cosine)),
            "explained_theoretical_fraction": explained_fraction,
            "fit_score": fit_score,
            "theoretical_intensity": expected,
        }
        rank = (fit_score, -offset)
        if best is None or rank > (best[0], best[1]):
            best = (fit_score, -offset, result)
    return best[2] if best is not None else {
        "monoisotopic_offset": 0,
        "cosine_score": 0.0,
        "explained_theoretical_fraction": 0.0,
        "fit_score": 0.0,
        "theoretical_intensity": [],
    }


@dataclass(frozen=True)
class PeptideCandidate:
    candidate_id: str
    base_peptide_id: str
    chain: str
    start: int
    end: int
    sequence: str
    neutral_mass: float
    modifications: tuple[tuple[int, str, float], ...] = ()
    carbamidomethyl_cys: bool = True
    proteolysis: str = "fully_tryptic"
    is_decoy: bool = False

    @property
    def modification_text(self) -> str:
        return "; ".join(f"{name}@{position + 1}" for position, name, _ in self.modifications) or "Unmodified"


def read_fasta(path: Path) -> dict[str, str]:
    chains: dict[str, str] = {}
    name = ""
    parts: list[str] = []

    def commit_record() -> None:
        if not name:
            return
        sequence = "".join(parts)
        if not sequence:
            raise ValueError(f"FASTA record has no amino acid sequence: {name}")
        chains[name] = sequence

    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            commit_record()
            name = line[1:].split()[0] or f"chain_{len(chains) + 1}"
            parts = []
        else:
            parts.append(line.upper().replace(" ", ""))
    commit_record()
    if not chains:
        raise ValueError(f"No FASTA sequences found: {path}")
    empty = sorted(name for name, sequence in chains.items() if not sequence)
    if empty:
        raise ValueError(f"FASTA contains empty sequence record(s): {empty}")
    invalid = sorted({aa for sequence in chains.values() for aa in sequence if aa not in AA_MASS})
    if invalid:
        raise ValueError(f"Unsupported amino acids in FASTA: {invalid}")
    return chains


def digest_trypsin(
    chains: dict[str, str],
    max_missed_cleavages: int = 2,
    min_length: int = 6,
    max_length: int = 60,
) -> list[tuple[str, int, int, str]]:
    peptides: list[tuple[str, int, int, str]] = []
    for chain, sequence in chains.items():
        cuts = [0]
        cuts.extend(
            index + 1
            for index, aa in enumerate(sequence)
            if aa in "KR" and (index + 1 == len(sequence) or sequence[index + 1] != "P")
        )
        if cuts[-1] != len(sequence):
            cuts.append(len(sequence))
        for left in range(len(cuts) - 1):
            for missed in range(max(0, max_missed_cleavages) + 1):
                right = left + missed + 1
                if right >= len(cuts):
                    break
                start, end = cuts[left], cuts[right]
                if min_length <= end - start <= max_length:
                    peptides.append((chain, start + 1, end, sequence[start:end]))
    return peptides


def _residue_masses(
    sequence: str,
    modifications: tuple[tuple[int, str, float], ...],
    carbamidomethyl_cys: bool = True,
) -> list[float]:
    deltas = {position: delta for position, _, delta in modifications}
    return [
        AA_MASS[aa] + (CARBAMIDOMETHYL if carbamidomethyl_cys and aa == "C" else 0.0) + deltas.get(index, 0.0)
        for index, aa in enumerate(sequence)
    ]


def _candidate(
    chain: str,
    start: int,
    end: int,
    sequence: str,
    modifications: tuple[tuple[int, str, float], ...] = (),
    carbamidomethyl_cys: bool = True,
    proteolysis: str = "fully_tryptic",
    is_decoy: bool = False,
    base_peptide_id: str | None = None,
) -> PeptideCandidate:
    base_id = base_peptide_id or f"{chain}:{start}-{end}:{sequence}"
    mod_text = ";".join(f"{name}@{position + 1}" for position, name, _ in modifications) or "Unmodified"
    prefix = "DECOY:" if is_decoy else ""
    return PeptideCandidate(
        candidate_id=f"{prefix}{base_id}:{mod_text}",
        base_peptide_id=base_id,
        chain=chain,
        start=start,
        end=end,
        sequence=sequence,
        neutral_mass=sum(_residue_masses(sequence, modifications, carbamidomethyl_cys)) + WATER,
        modifications=modifications,
        carbamidomethyl_cys=carbamidomethyl_cys,
        proteolysis=proteolysis,
        is_decoy=is_decoy,
    )


def _decoy(target: PeptideCandidate) -> PeptideCandidate:
    tail = target.sequence[-1] if target.sequence[-1] in "KR" else ""
    core = target.sequence[:-1] if tail else target.sequence
    decoy_sequence = core[::-1] + tail
    if decoy_sequence == target.sequence and len(core) > 1:
        decoy_sequence = core[1:] + core[:1] + tail
    core_length = len(core)
    modifications = tuple(
        (core_length - 1 - position if position < core_length else position, name, delta)
        for position, name, delta in target.modifications
    )
    return _candidate(
        target.chain,
        target.start,
        target.end,
        decoy_sequence,
        modifications,
        carbamidomethyl_cys=target.carbamidomethyl_cys,
        proteolysis=target.proteolysis,
        is_decoy=True,
        base_peptide_id=target.base_peptide_id,
    )


def _modification_sites(sequence: str) -> list[tuple[int, str, float]]:
    sites: list[tuple[int, str, float]] = []
    for position, aa in enumerate(sequence):
        for name, residues, delta, location in COMMON_MODIFICATIONS:
            location_matches = (
                location == "residue"
                or (location == "n_term" and position == 0)
                or (
                    location == "n_glycan_sequon"
                    and position + 2 < len(sequence)
                    and sequence[position + 1] != "P"
                    and sequence[position + 2] in "ST"
                )
            )
            if aa in residues and location_matches:
                sites.append((position, name, delta))
    return sites


def generate_candidates(
    chains: dict[str, str],
    max_missed_cleavages: int = 2,
    min_length: int = 6,
    max_length: int = 60,
    carbamidomethyl_cys: bool = True,
    max_variable_modifications: int = 1,
    semitryptic_max_trim: int = 0,
    include_decoys: bool = True,
) -> list[PeptideCandidate]:
    targets: list[PeptideCandidate] = []
    fully_tryptic = digest_trypsin(chains, max_missed_cleavages, min_length, max_length)
    digested: dict[tuple[str, int, int, str], str] = {
        row: "fully_tryptic" for row in fully_tryptic
    }
    for chain, start, end, sequence in fully_tryptic:
        for trim in range(1, max(0, semitryptic_max_trim) + 1):
            if len(sequence) - trim < min_length:
                break
            digested.setdefault((chain, start + trim, end, sequence[trim:]), "semi_tryptic")
            digested.setdefault((chain, start, end - trim, sequence[:-trim]), "semi_tryptic")

    for (chain, start, end, sequence), proteolysis in sorted(digested.items()):
        targets.append(_candidate(
            chain,
            start,
            end,
            sequence,
            carbamidomethyl_cys=carbamidomethyl_cys,
            proteolysis=proteolysis,
        ))
        sites = _modification_sites(sequence)
        for modification_count in range(1, max(0, max_variable_modifications) + 1):
            for selected in combinations(sites, modification_count):
                if len({position for position, _, __ in selected}) != modification_count:
                    continue
                targets.append(_candidate(
                    chain,
                    start,
                    end,
                    sequence,
                    tuple(selected),
                    carbamidomethyl_cys=carbamidomethyl_cys,
                    proteolysis=proteolysis,
                ))
        if end == len(chains[chain]) and sequence.endswith("K") and len(sequence) - 1 >= min_length:
            clipped = sequence[:-1]
            targets.append(_candidate(
                chain,
                start,
                end - 1,
                clipped,
                ((len(clipped) - 1, "C-terminal Lys clipping", 0.0),),
                carbamidomethyl_cys=carbamidomethyl_cys,
                proteolysis=proteolysis,
            ))
    return targets + ([_decoy(candidate) for candidate in targets] if include_decoys else [])


def generate_sequence_inference_candidates(
    chains: dict[str, str],
    max_missed_cleavages: int = 4,
    min_length: int = 6,
    max_length: int = 60,
    carbamidomethyl_cys: bool = True,
    max_terminal_trim: int = 20,
    include_decoys: bool = True,
) -> list[PeptideCandidate]:
    """Build unmodified backbone candidates including directional truncations."""
    digested: dict[tuple[str, int, int, str], str] = {
        row: "fully_tryptic"
        for row in digest_trypsin(
            chains,
            max_missed_cleavages,
            min_length,
            max_length,
        )
    }
    for chain, start, end, sequence in list(digested):
        for trim in range(1, max(0, max_terminal_trim) + 1):
            if len(sequence) - trim < min_length:
                break
            digested.setdefault(
                (chain, start + trim, end, sequence[trim:]),
                "n_terminal_truncation_candidate",
            )
            digested.setdefault(
                (chain, start, end - trim, sequence[:-trim]),
                "c_terminal_truncation_candidate",
            )
    targets = [
        _candidate(
            chain,
            start,
            end,
            sequence,
            carbamidomethyl_cys=carbamidomethyl_cys,
            proteolysis=proteolysis,
        )
        for (chain, start, end, sequence), proteolysis
        in sorted(digested.items())
    ]
    return targets + (
        [_decoy(candidate) for candidate in targets]
        if include_decoys else []
    )


@lru_cache(maxsize=None)
def theoretical_fragments(
    candidate: PeptideCandidate,
    max_charge: int = 2,
    include_extended: bool = False,
) -> tuple[tuple[str, float, str, int], ...]:
    masses = _residue_masses(candidate.sequence, candidate.modifications, candidate.carbamidomethyl_cys)
    total = sum(masses)
    prefix = 0.0
    ions: list[tuple[str, float, str, int]] = []
    for ordinal in range(1, len(masses)):
        prefix += masses[ordinal - 1]
        suffix = total - prefix
        for charge in range(1, max(1, max_charge) + 1):
            suffix_charge = "" if charge == 1 else f"^{charge}"
            b_mz = (prefix + charge * PROTON) / charge
            y_ordinal = len(masses) - ordinal
            y_mz = (suffix + WATER + charge * PROTON) / charge
            ions.append((f"b{ordinal}{suffix_charge}", b_mz, "b", ordinal))
            ions.append((f"y{y_ordinal}{suffix_charge}", y_mz, "y", y_ordinal))
            if include_extended:
                # HCD commonly produces a-ions and neutral-loss variants.  They
                # are scored as supporting evidence, while sequence coverage is
                # still calculated from the underlying b/y positions.
                if b_mz > 0.0:
                    ions.append((f"a{ordinal}{suffix_charge}", b_mz - 27.99491462 / charge, "a", ordinal))
                if b_mz > 18.010564684 / charge:
                    ions.append((f"b{ordinal}-H2O{suffix_charge}", b_mz - 18.010564684 / charge, "b", ordinal))
                if b_mz > 17.026549101 / charge:
                    ions.append((f"b{ordinal}-NH3{suffix_charge}", b_mz - 17.026549101 / charge, "b", ordinal))
                if y_mz > 18.010564684 / charge:
                    ions.append((f"y{y_ordinal}-H2O{suffix_charge}", y_mz - 18.010564684 / charge, "y", y_ordinal))
                if y_mz > 17.026549101 / charge:
                    ions.append((f"y{y_ordinal}-NH3{suffix_charge}", y_mz - 17.026549101 / charge, "y", y_ordinal))
    return tuple(ions)


def _longest_run(values: set[int]) -> int:
    longest = current = 0
    previous = None
    for value in sorted(values):
        current = current + 1 if previous is not None and value == previous + 1 else 1
        longest = max(longest, current)
        previous = value
    return longest


def match_fragments(
    scan: LCMSSpectrumScan,
    candidate: PeptideCandidate,
    precursor_error_ppm: float,
    precursor_tolerance_ppm: float,
    fragment_tolerance_ppm: float,
    max_fragment_charge: int = 2,
    include_extended_fragments: bool = True,
) -> dict[str, object]:
    pairs = sorted(zip(scan.mz_array, scan.intensity_array), key=lambda pair: pair[0])
    observed_mz = [pair[0] for pair in pairs]
    observed_intensity = [max(0.0, pair[1]) for pair in pairs]
    used: set[int] = set()
    matches: list[dict[str, object]] = []
    for label, theoretical_mz, series, ordinal in theoretical_fragments(
        candidate,
        max_charge=max_fragment_charge,
        include_extended=include_extended_fragments,
    ):
        tolerance = theoretical_mz * fragment_tolerance_ppm / 1_000_000.0
        left = bisect_left(observed_mz, theoretical_mz - tolerance)
        right = bisect_right(observed_mz, theoretical_mz + tolerance)
        choices = [index for index in range(left, right) if index not in used]
        if not choices:
            continue
        index = max(choices, key=observed_intensity.__getitem__)
        used.add(index)
        matches.append({
            "label": label,
            "series": series,
            "ordinal": ordinal,
            "theoretical_mz": theoretical_mz,
            "observed_mz": observed_mz[index],
            "error_ppm": (observed_mz[index] - theoretical_mz) / theoretical_mz * 1_000_000.0,
            "intensity": observed_intensity[index],
            "observed_index": index,
        })
    # Neutral-loss/a-ion support should not inflate sequence coverage.  Coverage
    # remains based on the canonical b/y backbone positions.
    positions = {
        (str(match["series"]), int(match["ordinal"]))
        for match in matches
        if str(match["series"]) in {"b", "y"}
    }
    denominator = max(1, 2 * (len(candidate.sequence) - 1))
    coverage = len(positions) / denominator
    total_intensity = sum(observed_intensity)
    explained_intensity = sum(observed_intensity[index] for index in used) / total_intensity if total_intensity else 0.0
    continuity = max(
        _longest_run({ordinal for series, ordinal in positions if series == "b"}),
        _longest_run({ordinal for series, ordinal in positions if series == "y"}),
    ) / max(1, len(candidate.sequence) - 1)
    mass_score = max(0.0, 1.0 - abs(precursor_error_ppm) / max(precursor_tolerance_ppm, 1e-9))
    score = 30.0 * mass_score + 35.0 * coverage + 25.0 * explained_intensity + 10.0 * continuity
    return {
        "score": score,
        "matched_ion_count": len(matches),
        "fragment_coverage": coverage,
        "explained_intensity": explained_intensity,
        "ion_continuity": continuity,
        "matched_fragments": matches,
    }


def search_scans(
    scans: list[LCMSSpectrumScan],
    candidates: list[PeptideCandidate],
    precursor_tolerance_ppm: float = 10.0,
    fragment_tolerance_ppm: float = 20.0,
    min_score: float = 20.0,
) -> list[dict[str, object]]:
    ordered = sorted(candidates, key=lambda candidate: candidate.neutral_mass)
    masses = [candidate.neutral_mass for candidate in ordered]
    psms: list[dict[str, object]] = []
    for scan in scans:
        if scan.ms_level != 2 or scan.precursor_mz is None or not scan.mz_array:
            continue
        best: dict[str, object] | None = None
        for charge in ([scan.precursor_charge] if scan.precursor_charge else range(1, 8)):
            if charge is None or charge <= 0:
                continue
            neutral_mass = (scan.precursor_mz - PROTON) * charge
            tolerance = neutral_mass * precursor_tolerance_ppm / 1_000_000.0
            for candidate in ordered[bisect_left(masses, neutral_mass - tolerance):bisect_right(masses, neutral_mass + tolerance)]:
                error_ppm = (neutral_mass - candidate.neutral_mass) / candidate.neutral_mass * 1_000_000.0
                evidence = match_fragments(scan, candidate, error_ppm, precursor_tolerance_ppm, fragment_tolerance_ppm)
                if float(evidence["score"]) < min_score:
                    continue
                row = {
                    "sample_id": scan.sample_id,
                    "scan_id": scan.scan_id,
                    "rt": scan.rt,
                    "precursor_mz": scan.precursor_mz,
                    "precursor_charge": charge,
                    "precursor_intensity": scan.precursor_intensity or 0.0,
                    "precursor_error_ppm": error_ppm,
                    "activation_method": scan.activation_method,
                    "collision_energy": scan.collision_energy,
                    **asdict(candidate),
                    "modification_text": candidate.modification_text,
                    **evidence,
                }
                if best is None or float(row["score"]) > float(best["score"]):
                    best = row
        if best is not None:
            # Keep the observed spectrum for every target candidate that passes
            # the search minimum.  Confidence is assigned later; discarding the
            # peaks below the C-level score threshold made D-level sequence
            # candidates look identified while leaving the report spectrum blank.
            if not best.get("is_decoy"):
                pairs = sorted(zip(scan.mz_array, scan.intensity_array), key=lambda pair: pair[0])
                labels = {
                    int(match["observed_index"]): str(match["label"])
                    for match in best["matched_fragments"]
                }
                top_indices = sorted(range(len(pairs)), key=lambda index: pairs[index][1], reverse=True)[:60]
                best["spectrum_peaks"] = [
                    {"mz": pairs[index][0], "intensity": pairs[index][1], "label": labels.get(index, "")}
                    for index in sorted(top_indices, key=lambda index: pairs[index][0])
                ]
            else:
                best["spectrum_peaks"] = []
            psms.append(best)
    return psms


def assign_q_values(psms: list[dict[str, object]]) -> list[dict[str, object]]:
    ordered = sorted(psms, key=lambda row: float(row["score"]), reverse=True)
    targets = decoys = 0
    fdrs: list[float] = []
    for row in ordered:
        if row.get("is_decoy"):
            decoys += 1
        else:
            targets += 1
        fdrs.append(decoys / max(1, targets))
    q_value = 1.0
    for row, fdr in zip(reversed(ordered), reversed(fdrs)):
        q_value = min(q_value, fdr)
        row["q_value"] = q_value
    return ordered


def read_peak_first_payload(path: Path) -> dict[str, object]:
    with sqlite3.connect(path) as connection:
        row = connection.execute(
            "SELECT payload_json FROM peak_first_artifacts WHERE artifact_key = 'bootstrap'"
        ).fetchone()
    if row is None:
        raise ValueError(f"Peak-first bootstrap is missing: {path}")
    return json.loads(row[0])


def annotate_feature_groups(
    psms: list[dict[str, object]],
    payload: dict[str, object],
    mz_tolerance_ppm: float = 20.0,
    rt_tolerance_min: float = 0.5,
    min_mz_tolerance_da: float = 0.02,
    min_isotope_offset: int = -1,
    max_isotope_offset: int = 4,
    min_feature_charge: int = 1,
    max_feature_charge: int = 7,
) -> None:
    rows = sorted(payload.get("global_feature_groups") or [], key=lambda row: float(row.get("representative_mz") or 0.0))
    mz_values = [float(row.get("representative_mz") or 0.0) for row in rows]
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    source_groups = {
        str(group.get("feature_group_id")): group
        for peak in payload.get("peak_results") or []
        for group in peak.get("feature_groups") or []
        if group.get("feature_group_id")
    }

    def coordinates(row: dict[str, object], sample_id: str) -> list[tuple[float, float, float]]:
        values: list[tuple[float, float, float]] = [(
            float(row.get("representative_mz") or 0.0),
            float(row.get("representative_rt") or 0.0),
            float(dict(row.get("rt_correction_by_sample") or {}).get(sample_id) or 0.0),
        )]
        source_ids = row.get("source_feature_group_ids") or [row.get("feature_group_id")]
        for source_id in source_ids:
            source = source_groups.get(str(source_id))
            if not source:
                continue
            sample_feature = dict(source.get("features_by_sample") or {}).get(sample_id) or {}
            feature_mz = float(sample_feature.get("mz") or source.get("representative_mz") or row.get("representative_mz") or 0.0)
            feature_rt = float(sample_feature.get("aligned_rt_apex") or source.get("representative_rt") or row.get("representative_rt") or 0.0)
            local_shift = float(dict(source.get("rt_correction_by_sample") or {}).get(sample_id) or 0.0)
            values.append((feature_mz, feature_rt, local_shift))
        return values

    for psm in psms:
        mz = float(psm["precursor_mz"])
        charge = int(psm.get("precursor_charge") or 0)
        sample_id = str(psm["sample_id"])
        if charge <= 0:
            continue
        neutral_mass = float(psm.get("neutral_mass") or ((mz - PROTON) * charge))
        links_by_feature: dict[str, dict[str, object]] = {}
        for feature_charge in range(min_feature_charge, max_feature_charge + 1):
            for isotope_offset in range(min_isotope_offset, max_isotope_offset + 1):
                expected_mz = (neutral_mass + isotope_offset * ISOTOPE_MASS_DIFF) / feature_charge + PROTON
                tolerance = max(expected_mz * mz_tolerance_ppm / 1_000_000.0, min_mz_tolerance_da)
                left = bisect_left(mz_values, expected_mz - tolerance - min_mz_tolerance_da)
                right = bisect_right(mz_values, expected_mz + tolerance + min_mz_tolerance_da)
                for row in rows[left:right]:
                    best_coordinate: tuple[float, float, float, float] | None = None
                    for feature_mz, feature_rt, local_shift in coordinates(row, sample_id):
                        mz_error = feature_mz - expected_mz
                        if abs(mz_error) > tolerance:
                            continue
                        aligned_rt = float(psm["rt"]) + float(shifts.get(sample_id) or 0.0) + local_shift
                        rt_error = aligned_rt - feature_rt
                        if abs(rt_error) > rt_tolerance_min:
                            continue
                        distance = (
                            abs(rt_error) / max(rt_tolerance_min, 1e-9)
                            + abs(mz_error) / max(tolerance, 1e-9)
                            + 0.02 * abs(isotope_offset)
                            + 0.01 * abs(feature_charge - charge)
                        )
                        if best_coordinate is None or distance < best_coordinate[0]:
                            best_coordinate = (distance, feature_mz, feature_rt, mz_error)
                    if best_coordinate is None:
                        continue
                    distance, sample_feature_mz, sample_feature_rt, mz_error = best_coordinate
                    feature_id = str(row.get("feature_group_id"))
                    if feature_charge != charge:
                        link_type = "charge_state_envelope"
                    elif isotope_offset != 0:
                        link_type = "isotope_envelope"
                    else:
                        link_type = "direct_precursor"
                    link = {
                        "feature_group_id": row.get("feature_group_id"),
                        "feature_mz": row.get("representative_mz"),
                        "feature_sample_mz": sample_feature_mz,
                        "feature_rt": row.get("representative_rt"),
                        "feature_sample_rt": sample_feature_rt,
                        "feature_similarity": row.get("similarity_score"),
                        "feature_difference": row.get("difference_score"),
                        "feature_ranking": row.get("ranking_score"),
                        "feature_difference_type": row.get("difference_type"),
                        "feature_area_by_sample": row.get("area_by_sample") or {},
                        "feature_normalized_area_by_sample": row.get("normalized_area_by_sample") or {},
                        "feature_link_type": link_type,
                        "feature_charge": feature_charge,
                        "feature_isotope_offset": isotope_offset,
                        "feature_mz_error_da": mz_error,
                        "feature_mz_error_ppm": mz_error / max(expected_mz, 1e-12) * 1_000_000.0,
                        "feature_link_distance": distance,
                    }
                    previous = links_by_feature.get(feature_id)
                    if previous is None or distance < float(previous["feature_link_distance"]):
                        links_by_feature[feature_id] = link
        links = sorted(links_by_feature.values(), key=lambda link: float(link["feature_link_distance"]))
        if links:
            psm["feature_links"] = links
            psm.update({
                key: value
                for key, value in links[0].items()
                if key != "feature_link_distance"
            })


def psm_confidence(psm: dict[str, object]) -> str:
    open_confidence = str(psm.get("open_modification_confidence") or "")
    if open_confidence == "B_open_modification_supported":
        return "B_high_confidence_inferred"
    if open_confidence == "C_open_mass_candidate":
        return "C_tentative"
    if (
        float(psm.get("q_value", 1.0)) <= 0.01
        and float(psm.get("score") or 0.0) >= 60.0
        and int(psm.get("matched_ion_count") or 0) >= 6
        and float(psm.get("fragment_coverage") or 0.0) >= 0.15
    ):
        return "B_high_confidence_inferred"
    if float(psm.get("score") or 0.0) >= 40.0 and int(psm.get("matched_ion_count") or 0) >= 4:
        return "C_tentative"
    return "D_low_evidence"


def differential_annotations(psms: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    for psm in psms:
        links = psm.get("feature_links") or []
        if links:
            for link in links:
                linked_psm = dict(psm)
                linked_psm.update(link)
                group_id = str(link.get("feature_group_id") or psm["candidate_id"])
                grouped.setdefault((group_id, str(psm["candidate_id"])), []).append(linked_psm)
        else:
            group_id = str(psm.get("feature_group_id") or psm["candidate_id"])
            grouped.setdefault((group_id, str(psm["candidate_id"])), []).append(psm)
    annotations: list[dict[str, object]] = []
    for (group_id, candidate_id), evidence in grouped.items():
        best_by_sample: dict[str, dict[str, object]] = {}
        for psm in evidence:
            sample_id = str(psm["sample_id"])
            if sample_id not in best_by_sample or float(psm["score"]) > float(best_by_sample[sample_id]["score"]):
                best_by_sample[sample_id] = psm
        best = max(evidence, key=lambda row: float(row["score"]))
        modification = str(best["modification_text"])
        interpretation = f"修饰肽候选：{modification}" if modification != "Unmodified" else "序列肽段候选"
        link_type = str(best.get("feature_link_type") or "unlinked_sequence")
        isotope_offset = int(best.get("feature_isotope_offset") or 0)
        if link_type == "isotope_envelope":
            isotope_label = f"M+{isotope_offset}" if isotope_offset >= 0 else f"M{isotope_offset}"
            interpretation = f"同一肽的 {isotope_label} 同位素峰，不是独立修饰；{interpretation}"
        elif link_type == "charge_state_envelope":
            isotope_label = f"M+{isotope_offset}" if isotope_offset > 0 else (f"M{isotope_offset}" if isotope_offset < 0 else "M")
            interpretation = f"同一肽的 z={best.get('feature_charge')}、{isotope_label} 电荷/同位素峰，不是独立修饰；{interpretation}"
        elif link_type == "isolation_window_targeted_search":
            interpretation = f"隔离窗内 MS1 引导定向重搜：{interpretation}；存在共隔离/嵌合谱歧义"
        elif link_type == "selected_precursor_consensus_search":
            interpretation = f"重复采集的选中前体 MS2 共识谱推测：{interpretation}；需用单次高质量谱或靶向采集复核"
        elif link_type == "component_charge_isotope_consensus_search":
            interpretation = (
                "同一 MS1 组件的多电荷态/同位素成员联合 b/y 证据推测："
                f"{interpretation}；组件中性质量用于候选筛选，需以独立高质量谱复核"
            )
        elif link_type == "component_mass_offset_sequence_tag_search":
            inference_level = str(
                best.get("sequence_inference_level")
                or "backbone_sequence_supported"
            )
            if inference_level == "truncation_sequence_supported":
                interpretation = (
                    "已知轻/重链序列约束的截断候选：组件中性质量、"
                    "连续 b/y 标签及同组多谱共同支持该截断骨架；"
                    f"{interpretation}，仍需独立高质量谱或肽段标准品复核"
                )
            else:
                interpretation = (
                    "质量偏移容忍的骨架序列推测：组件中性质量用于计算"
                    "未知 ΔMass，并由连续 b/y 标签和同组多谱选择候选；"
                    f"{interpretation}，修饰类别/位点可能仍未确定"
                )
        elif link_type == "feature_open_mass_search":
            interpretation = (
                f"单 Feature 开放质量搜索：{interpretation}；"
                "未知 ΔMass 由 b/y 及扩展碎片支持，需结合标准品或进一步采集复核"
            )
        elif link_type == "feature_targeted_glycopeptide_search":
            interpretation = (
                f"差异 Feature 定向糖肽搜索：{interpretation}；"
                f"检测到 {int(best.get('glycan_diagnostic_ion_count') or 0)} 个糖链诊断离子、"
                f"{int(best.get('glycan_core_y_ion_count') or 0)} 个核心 Y 离子及"
                f"{int(best.get('matched_ion_count') or 0)} 个肽骨架位点证据；"
                "糖型和位点仍建议结合标准品或正交糖分析复核"
            )
        elif link_type == "global_open_modification_search":
            interpretation = (
                f"全局开放修饰搜索：{interpretation}；该 ΔMass 经重复谱聚类、"
                "目标-诱饵控制及限制性二次重搜支持，仍建议结合正交方法确认具体化学结构"
            )
        if len(best_by_sample) == 1:
            interpretation += f"；仅 {next(iter(best_by_sample))} 获得合格 MS/MS"
        else:
            interpretation += "；两样品均获得 MS/MS"
        confidence = psm_confidence(best)
        if link_type in {
            "isolation_window_targeted_search",
            "selected_precursor_consensus_search",
            "component_charge_isotope_consensus_search",
            "component_mass_offset_sequence_tag_search",
            "feature_open_mass_search",
            "feature_targeted_glycopeptide_search",
            "global_open_modification_search",
        } and confidence.startswith("B_"):
            if link_type != "global_open_modification_search":
                confidence = "C_tentative"
        annotations.append({
            "feature_or_candidate_id": group_id,
            "candidate_id": candidate_id,
            "base_peptide_id": best["base_peptide_id"],
            "neutral_mass": best.get("neutral_mass"),
            "chain": best["chain"],
            "start": best["start"],
            "end": best["end"],
            "sequence": best["sequence"],
            "modification": modification,
            "proteolysis": best.get("proteolysis", "fully_tryptic"),
            "confidence": confidence,
            "interpretation": interpretation,
            "feature_mz": best.get("feature_mz", best.get("precursor_mz")),
            "feature_rt": best.get("feature_rt", best.get("rt")),
            "feature_difference": best.get("feature_difference"),
            "feature_ranking": best.get("feature_ranking"),
            "area_by_sample": best.get("feature_area_by_sample") or {},
            "normalized_area_by_sample": best.get("feature_normalized_area_by_sample") or {},
            "feature_link_type": link_type,
            "feature_isotope_offset": isotope_offset,
            "feature_charge": best.get("feature_charge"),
            "feature_mz_error_da": best.get("feature_mz_error_da"),
            "feature_mz_error_ppm": best.get("feature_mz_error_ppm"),
            "search_origin": best.get("search_origin", "standard_sequence_search"),
            "sequence_inference_level": best.get("sequence_inference_level"),
            "mass_delta": best.get("mass_delta"),
            "backbone_neutral_mass": best.get("backbone_neutral_mass"),
            "sequence_tag_length": best.get("sequence_tag_length"),
            "complementary_ion_pair_count": best.get(
                "complementary_ion_pair_count"
            ),
            "candidate_score_margin": best.get("candidate_score_margin"),
            "sample_evidence": {
                sample_id: {
                    "scan_id": row["scan_id"],
                    "rt": row["rt"],
                    "score": row["score"],
                    "q_value": row.get("q_value"),
                    "matched_ion_count": row["matched_ion_count"],
                    "fragment_coverage": row["fragment_coverage"],
                }
                for sample_id, row in best_by_sample.items()
            },
        })
    return sorted(annotations, key=lambda row: (float(row.get("feature_ranking") or 0.0), max(float(item["score"]) for item in row["sample_evidence"].values())), reverse=True)


def _infer_unidentified_ms1_groups(
    payload: dict[str, object],
    feature_rows: dict[str, dict[str, object]],
    excluded_feature_ids: set[str],
    spectra_loader: Callable[[str], list[dict[str, object]]] | None,
    mass_tolerance_ppm: float = 25.0,
    rt_tolerance_min: float = 0.18,
    min_xic_shape_score: float = 0.82,
    max_xic_apex_delta_min: float = 0.10,
    min_inference_score: float = 0.82,
) -> tuple[list[dict[str, object]], set[str]]:
    """Infer unknown MS1 components from mass, coelution, and abundance evidence.

    The result deliberately remains an MS1-level inference. It is never promoted
    to a peptide identification and only forms a group when a narrow-mass
    charge/isotope relationship is corroborated by a coeluting raw XIC pair.
    """
    rows = [
        row for feature_id, row in feature_rows.items()
        if feature_id not in excluded_feature_ids
        and float(row.get("representative_mz") or 0.0) > PROTON
        and float(row.get("representative_rt") or 0.0) > 0.0
    ]
    if len(rows) < 2:
        return [], set()

    source_groups = {
        str(group.get("feature_group_id") or ""): group
        for peak in payload.get("peak_results") or []
        for group in peak.get("feature_groups") or []
        if group.get("feature_group_id")
    }
    sample_ids = [str(sample_id) for sample_id in payload.get("sample_ids") or []]

    def source_observations(row: dict[str, object]) -> dict[str, dict[str, float]]:
        output: dict[str, dict[str, float]] = {}
        source_ids = (
            row.get("merged_feature_group_ids")
            or row.get("source_feature_group_ids")
            or [row.get("feature_group_id")]
        )
        for source_id in source_ids:
            source = source_groups.get(str(source_id))
            if not source:
                continue
            for sample_id, feature in dict(source.get("features_by_sample") or {}).items():
                if not isinstance(feature, dict):
                    continue
                observation = {
                    "rt_apex": float(feature.get("aligned_rt_apex") or feature.get("rt_apex") or row.get("representative_rt") or 0.0),
                    "rt_start": float(feature.get("rt_start") or row.get("representative_rt") or 0.0),
                    "rt_end": float(feature.get("rt_end") or row.get("representative_rt") or 0.0),
                    "area": float(feature.get("normalized_area") or feature.get("area") or 0.0),
                }
                previous = output.get(str(sample_id))
                if previous is None or observation["area"] > previous["area"]:
                    output[str(sample_id)] = observation
        if not output:
            values = dict(row.get("normalized_area_by_sample") or row.get("area_by_sample") or {})
            for sample_id, value in values.items():
                if float(value or 0.0) <= 0.0:
                    continue
                rt = float(row.get("representative_rt") or 0.0)
                output[str(sample_id)] = {"rt_apex": rt, "rt_start": rt - 0.06, "rt_end": rt + 0.06, "area": float(value)}
        return output

    observations = {
        str(row.get("feature_group_id")): source_observations(row)
        for row in rows
    }

    def centered_log_profile(row: dict[str, object]) -> list[float]:
        values = dict(row.get("normalized_area_by_sample") or row.get("area_by_sample") or {})
        ordered = [max(0.0, float(values.get(sample_id) or 0.0)) for sample_id in sample_ids]
        maximum = max(ordered, default=0.0)
        if maximum <= 0.0:
            return [0.0 for _ in ordered]
        floor = max(maximum * 0.005, 1e-12)
        logs = [math.log(max(value, floor)) for value in ordered]
        center = statistics.mean(logs)
        return [value - center for value in logs]

    log_profiles = {
        str(row.get("feature_group_id")): centered_log_profile(row)
        for row in rows
    }

    def trend_score(first_id: str, second_id: str) -> float:
        first = log_profiles[first_id]
        second = log_profiles[second_id]
        if len(first) < 2 or len(second) != len(first):
            return 0.75
        first_norm = math.sqrt(sum(value * value for value in first))
        second_norm = math.sqrt(sum(value * value for value in second))
        if first_norm < 0.10 and second_norm < 0.10:
            return 1.0
        if first_norm < 0.10 or second_norm < 0.10:
            return 0.55
        direction = sum(a * b for a, b in zip(first, second)) / max(first_norm * second_norm, 1e-12)
        if direction <= 0.0:
            return 0.0
        rms_delta = math.sqrt(statistics.mean((a - b) ** 2 for a, b in zip(first, second)))
        return max(0.0, min(1.0, direction * math.exp(-rms_delta / 1.5)))

    def rt_evidence(
        first: dict[str, object],
        second: dict[str, object],
    ) -> tuple[float, dict[str, float]] | None:
        first_id = str(first.get("feature_group_id"))
        second_id = str(second.get("feature_group_id"))
        representative_delta = abs(float(first.get("representative_rt") or 0.0) - float(second.get("representative_rt") or 0.0))
        shared = set(observations[first_id]) & set(observations[second_id])
        centers: dict[str, float] = {}
        apex_deltas: list[float] = []
        overlap_scores: list[float] = []
        for sample_id in shared:
            left = observations[first_id][sample_id]
            right = observations[second_id][sample_id]
            delta = abs(left["rt_apex"] - right["rt_apex"])
            apex_deltas.append(delta)
            centers[sample_id] = (left["rt_apex"] + right["rt_apex"]) * 0.5
            intersection = max(0.0, min(left["rt_end"], right["rt_end"]) - max(left["rt_start"], right["rt_start"]))
            smaller_width = max(1e-9, min(left["rt_end"] - left["rt_start"], right["rt_end"] - right["rt_start"]))
            overlap_scores.append(min(1.0, intersection / smaller_width))
        best_apex_delta = min(apex_deltas, default=float("inf"))
        if representative_delta > rt_tolerance_min and best_apex_delta > 0.12:
            return None
        if not shared:
            centers = {
                sample_id: statistics.mean([
                    float(first.get("representative_rt") or 0.0),
                    float(second.get("representative_rt") or 0.0),
                ])
                for sample_id in sample_ids
            }
        apex_score = math.exp(-min(apex_deltas, default=representative_delta) / 0.08)
        overlap_score = max(overlap_scores, default=apex_score)
        score = max(0.0, min(1.0, 0.55 * apex_score + 0.30 * overlap_score + 0.15 * math.exp(-representative_delta / 0.12)))
        return score, centers

    refined_envelopes: dict[str, dict[str, object]] = {}

    def resolve_isotope_envelope(
        spectra: list[dict[str, object]],
        scan_rts: list[float],
        row: dict[str, object],
        center_rt: float,
        sample_id: str,
    ) -> dict[str, object] | None:
        target_mz = float(row.get("representative_mz") or 0.0)
        left = bisect_left(scan_rts, center_rt - 0.10)
        right = bisect_right(scan_rts, center_rt + 0.10)
        local_scans = spectra[left:right]
        if not local_scans:
            return None
        apex_scan = max(
            local_scans,
            key=lambda scan: sum(
                float(intensity)
                for mz, intensity in zip(scan.get("mz") or [], scan.get("intensity") or [])
                if abs(float(mz) - target_mz) <= 0.65
            ),
        )
        points = sorted(
            (
                (float(mz), float(intensity))
                for mz, intensity in zip(apex_scan.get("mz") or [], apex_scan.get("intensity") or [])
                if abs(float(mz) - target_mz) <= 2.2 and float(intensity) > 0.0
            ),
            key=lambda item: item[0],
        )
        if len(points) < 2:
            return None
        maximum = max(intensity for _, intensity in points)
        points = [item for item in points if item[1] >= maximum * 0.005]
        mz_values = [item[0] for item in points]
        best: tuple[float, dict[str, object]] | None = None

        def nearest_index(expected: float, tolerance: float) -> int | None:
            index = bisect_left(mz_values, expected)
            choices = [candidate for candidate in (index - 1, index) if 0 <= candidate < len(mz_values)]
            if not choices:
                return None
            match = min(choices, key=lambda candidate: abs(mz_values[candidate] - expected))
            return match if abs(mz_values[match] - expected) <= tolerance else None

        for charge in range(1, 7):
            spacing = ISOTOPE_MASS_DIFF / charge
            tolerance = max(0.008, target_mz * 20.0 / 1_000_000.0)
            for anchor_index, (anchor_mz, _) in enumerate(points):
                if abs(anchor_mz - target_mz) > 0.70:
                    continue
                sequence = {anchor_index}
                cursor_mz = anchor_mz
                for _ in range(4):
                    match = nearest_index(cursor_mz - spacing, tolerance)
                    if match is None or match in sequence:
                        break
                    sequence.add(match)
                    cursor_mz = mz_values[match]
                cursor_mz = anchor_mz
                for _ in range(4):
                    match = nearest_index(cursor_mz + spacing, tolerance)
                    if match is None or match in sequence:
                        break
                    sequence.add(match)
                    cursor_mz = mz_values[match]
                ordered_indices = sorted(sequence)
                if len(ordered_indices) < 2:
                    continue
                envelope = [points[index] for index in ordered_indices]
                spacing_errors = [
                    abs((envelope[index][0] - envelope[index - 1][0]) - spacing)
                    for index in range(1, len(envelope))
                ]
                mean_error = statistics.mean(spacing_errors)
                intensity_fraction = sum(intensity for _, intensity in envelope) / max(
                    sum(intensity for _, intensity in points),
                    1e-12,
                )
                observed_first_mass = (envelope[0][0] - PROTON) * charge
                isotope_fit = fit_theoretical_isotope_envelope(
                    [intensity for _, intensity in envelope],
                    observed_first_mass,
                )
                monoisotopic_offset = int(isotope_fit["monoisotopic_offset"])
                monoisotopic_mz = envelope[0][0] - monoisotopic_offset * spacing
                charge_prior = 1.0 if charge in {2, 3} else (0.92 if charge == 4 else 0.78)
                score = (
                    len(envelope)
                    + 0.45 * intensity_fraction
                    + 0.15 * charge_prior
                    + 1.10 * float(isotope_fit["fit_score"])
                    - 0.50 * mean_error / max(tolerance, 1e-12)
                )
                candidate = {
                    "charge": charge,
                    "observed_first_isotope_mz": envelope[0][0],
                    "monoisotopic_mz": monoisotopic_mz,
                    "monoisotopic_offset": monoisotopic_offset,
                    "target_isotope_index": min(
                        range(len(envelope)),
                        key=lambda index: abs(envelope[index][0] - target_mz),
                    ),
                    "neutral_mass_observed_first": observed_first_mass,
                    "neutral_mass_monoisotopic": (monoisotopic_mz - PROTON) * charge,
                    "envelope_mz": [mz for mz, _ in envelope],
                    "envelope_intensity": [intensity for _, intensity in envelope],
                    "isotope_peak_count": len(envelope),
                    "mean_spacing_error_da": mean_error,
                    "isotope_fit_score": isotope_fit["fit_score"],
                    "isotope_fit_cosine": isotope_fit["cosine_score"],
                    "isotope_fit_explained_fraction": isotope_fit["explained_theoretical_fraction"],
                    "theoretical_isotope_intensity": isotope_fit["theoretical_intensity"],
                    "source_sample_id": sample_id,
                    "source_scan_id": apex_scan.get("scan_id"),
                    "source_rt": float(apex_scan.get("aligned_rt") or apex_scan.get("rt") or center_rt),
                    "envelope_score": score,
                }
                if best is None or score > best[0]:
                    best = (score, candidate)
        return best[1] if best is not None else None

    if spectra_loader is not None:
        rows_by_strongest_sample: dict[str, list[dict[str, object]]] = {}
        for row in rows:
            feature_id = str(row.get("feature_group_id"))
            values = dict(row.get("normalized_area_by_sample") or row.get("area_by_sample") or {})
            candidate_samples = set(observations[feature_id]) or set(sample_ids)
            if not candidate_samples:
                continue
            strongest_sample = max(
                candidate_samples,
                key=lambda sample_id: float(values.get(sample_id) or observations[feature_id].get(sample_id, {}).get("area") or 0.0),
            )
            rows_by_strongest_sample.setdefault(str(strongest_sample), []).append(row)
        for sample_id, sample_rows in rows_by_strongest_sample.items():
            spectra = sorted(
                spectra_loader(sample_id),
                key=lambda scan: float(scan.get("aligned_rt") or scan.get("rt") or 0.0),
            )
            scan_rts = [float(scan.get("aligned_rt") or scan.get("rt") or 0.0) for scan in spectra]
            for row in sample_rows:
                feature_id = str(row.get("feature_group_id"))
                center_rt = float(
                    observations[feature_id].get(sample_id, {}).get("rt_apex")
                    or row.get("representative_rt")
                    or 0.0
                )
                envelope = resolve_isotope_envelope(spectra, scan_rts, row, center_rt, sample_id)
                if envelope is not None:
                    refined_envelopes[feature_id] = envelope

    hypotheses: dict[str, list[tuple[int, int, float]]] = {}
    for row in rows:
        feature_id = str(row.get("feature_group_id"))
        mz = float(row.get("representative_mz") or 0.0)
        refined = refined_envelopes.get(feature_id)
        if refined is not None:
            charge = int(refined["charge"])
            observed_mass = float(refined["neutral_mass_observed_first"])
            best_offset = int(refined.get("monoisotopic_offset") or 0)
            fit_score = float(refined.get("isotope_fit_score") or 0.0)
            offsets = (
                [best_offset]
                if fit_score >= 0.72
                else sorted(range(0, 5), key=lambda offset: (abs(offset - best_offset), offset))
            )
            hypotheses[feature_id] = [
                (charge, isotope_offset, observed_mass - isotope_offset * ISOTOPE_MASS_DIFF)
                for isotope_offset in offsets
                if 500.0 <= observed_mass - isotope_offset * ISOTOPE_MASS_DIFF <= 15_000.0
            ]
        else:
            hypotheses[feature_id] = [
                (charge, isotope_offset, (mz - PROTON) * charge - isotope_offset * ISOTOPE_MASS_DIFF)
                for charge in range(1, 7)
                for isotope_offset in range(0, 5)
                if 500.0 <= (mz - PROTON) * charge - isotope_offset * ISOTOPE_MASS_DIFF <= 15_000.0
            ]

    def mass_relation(first: dict[str, object], second: dict[str, object]) -> dict[str, object] | None:
        first_id = str(first.get("feature_group_id"))
        second_id = str(second.get("feature_group_id"))
        first_mz = float(first.get("representative_mz") or 0.0)
        second_mz = float(second.get("representative_mz") or 0.0)
        best: tuple[float, dict[str, object]] | None = None

        low_id, high_id = (first_id, second_id) if first_mz <= second_mz else (second_id, first_id)
        low_mz, high_mz = sorted((first_mz, second_mz))
        delta = high_mz - low_mz
        for charge in range(1, 7):
            for isotope_step in range(1, 5):
                expected_delta = isotope_step * ISOTOPE_MASS_DIFF / charge
                mz_error = abs(delta - expected_delta)
                mz_tolerance = max(0.008, high_mz * 20.0 / 1_000_000.0)
                if mz_error > mz_tolerance:
                    continue
                for low_offset in range(0, 5 - isotope_step):
                    high_offset = low_offset + isotope_step
                    low_mass = (low_mz - PROTON) * charge - low_offset * ISOTOPE_MASS_DIFF
                    high_mass = (high_mz - PROTON) * charge - high_offset * ISOTOPE_MASS_DIFF
                    neutral_mass = (low_mass + high_mass) * 0.5
                    if not 500.0 <= neutral_mass <= 15_000.0:
                        continue
                    ppm_error = abs(low_mass - high_mass) / neutral_mass * 1_000_000.0
                    if ppm_error > mass_tolerance_ppm:
                        continue
                    charge_prior = 1.0 if charge in {2, 3} else (0.92 if charge == 4 else 0.80)
                    score = 0.78 * (1.0 - ppm_error / mass_tolerance_ppm) + 0.14 * charge_prior + 0.08 * (1.0 - low_offset / 4.0)
                    relation = {
                        "relation_type": "same_charge_isotope",
                        "neutral_mass": neutral_mass,
                        "mass_error_ppm": ppm_error,
                        "mass_score": max(0.0, min(1.0, score)),
                        "mass_evidence": "coarse_feature_isotope_spacing",
                        "assignments": {
                            low_id: {"charge": charge, "isotope_offset": low_offset},
                            high_id: {"charge": charge, "isotope_offset": high_offset},
                        },
                    }
                    if best is None or score > best[0]:
                        best = (score, relation)

        for first_charge, first_offset, first_mass in hypotheses[first_id]:
            for second_charge, second_offset, second_mass in hypotheses[second_id]:
                refined_pair = first_id in refined_envelopes and second_id in refined_envelopes
                if first_charge == second_charge and not refined_pair:
                    continue
                neutral_mass = (first_mass + second_mass) * 0.5
                ppm_error = abs(first_mass - second_mass) / max(neutral_mass, 1e-12) * 1_000_000.0
                if ppm_error > mass_tolerance_ppm:
                    continue
                charge_prior = statistics.mean(
                    1.0 if charge in {2, 3} else (0.92 if charge == 4 else 0.78)
                    for charge in (first_charge, second_charge)
                )
                offset_prior = 1.0 - (first_offset + second_offset) / 8.0
                score = 0.80 * (1.0 - ppm_error / mass_tolerance_ppm) + 0.14 * charge_prior + 0.06 * offset_prior
                relation = {
                    "relation_type": (
                        "same_charge_isotope"
                        if first_charge == second_charge
                        else "different_charge_state"
                    ),
                    "neutral_mass": neutral_mass,
                    "mass_error_ppm": ppm_error,
                    "mass_score": max(0.0, min(1.0, score)),
                    "mass_evidence": (
                        "resolved_isotope_envelopes"
                        if refined_pair else "coarse_feature_mz"
                    ),
                    "assignments": {
                        first_id: {
                            "charge": first_charge,
                            "isotope_offset": first_offset + int(
                                dict(refined_envelopes.get(first_id) or {}).get("target_isotope_index") or 0
                            ),
                        },
                        second_id: {
                            "charge": second_charge,
                            "isotope_offset": second_offset + int(
                                dict(refined_envelopes.get(second_id) or {}).get("target_isotope_index") or 0
                            ),
                        },
                    },
                }
                if best is None or score > best[0]:
                    best = (score, relation)
        return best[1] if best is not None else None

    ordered_rows = sorted(rows, key=lambda row: float(row.get("representative_rt") or 0.0))
    candidate_edges: list[dict[str, object]] = []
    for left_index, first in enumerate(ordered_rows):
        first_rt = float(first.get("representative_rt") or 0.0)
        for second in ordered_rows[left_index + 1:]:
            second_rt = float(second.get("representative_rt") or 0.0)
            if second_rt - first_rt > max(rt_tolerance_min, 0.40):
                break
            rt_result = rt_evidence(first, second)
            if rt_result is None:
                continue
            relation = mass_relation(first, second)
            if relation is None:
                continue
            first_id = str(first.get("feature_group_id"))
            second_id = str(second.get("feature_group_id"))
            consistency = trend_score(first_id, second_id)
            if consistency < 0.50:
                continue
            rt_score, centers = rt_result
            shared_samples = list(centers)
            if not shared_samples:
                continue
            area_first = dict(first.get("normalized_area_by_sample") or first.get("area_by_sample") or {})
            area_second = dict(second.get("normalized_area_by_sample") or second.get("area_by_sample") or {})
            best_sample = max(
                shared_samples,
                key=lambda sample_id: min(float(area_first.get(sample_id) or 0.0), float(area_second.get(sample_id) or 0.0)),
            )
            candidate_edges.append({
                **relation,
                "first_id": first_id,
                "second_id": second_id,
                "rt_score": rt_score,
                "trend_score": consistency,
                "xic_sample_id": best_sample,
                "xic_center_rt": centers[best_sample],
            })

    def xic_evidence(
        spectra: list[dict[str, object]],
        scan_rts: list[float],
        first_targets: list[float],
        second_targets: list[float],
        center_rt: float,
    ) -> tuple[float, float] | None:
        left = bisect_left(scan_rts, center_rt - 0.40)
        right = bisect_right(scan_rts, center_rt + 0.40)
        local = spectra[left:right]
        if len(local) < 5:
            return None
        traces = ([], [])
        for scan in local:
            mz_values = scan.get("mz") or []
            intensity_values = scan.get("intensity") or []
            for trace, targets in zip(traces, (first_targets, second_targets)):
                trace.append(sum(
                    float(intensity)
                    for mz, intensity in zip(mz_values, intensity_values)
                    if any(
                        abs(float(mz) - target_mz)
                        <= max(0.02, target_mz * 20.0 / 1_000_000.0)
                        for target_mz in targets
                    )
                ))
        if any(max(trace, default=0.0) <= 0.0 or sum(value > 0.0 for value in trace) < 3 for trace in traces):
            return None
        transformed = [[math.sqrt(max(0.0, value)) for value in trace] for trace in traces]
        dot = sum(a * b for a, b in zip(*transformed))
        norm = math.sqrt(sum(value * value for value in transformed[0]) * sum(value * value for value in transformed[1]))
        cosine_score = dot / max(norm, 1e-12)
        active_sets = []
        for trace in traces:
            threshold = max(trace) * 0.02
            active_sets.append({index for index, value in enumerate(trace) if value >= threshold})
        union = active_sets[0] | active_sets[1]
        overlap = len(active_sets[0] & active_sets[1]) / max(1, len(union))
        apex_indices = [max(range(len(trace)), key=trace.__getitem__) for trace in traces]
        apex_delta = abs(float(local[apex_indices[0]].get("aligned_rt") or local[apex_indices[0]].get("rt") or 0.0) - float(local[apex_indices[1]].get("aligned_rt") or local[apex_indices[1]].get("rt") or 0.0))
        return max(0.0, min(1.0, 0.82 * cosine_score + 0.18 * overlap)), apex_delta

    accepted_edges: list[dict[str, object]] = []
    if spectra_loader is not None and candidate_edges:
        row_by_id = {str(row.get("feature_group_id")): row for row in rows}
        for sample_id in sorted({str(edge["xic_sample_id"]) for edge in candidate_edges}):
            spectra = sorted(
                spectra_loader(sample_id),
                key=lambda scan: float(scan.get("aligned_rt") or scan.get("rt") or 0.0),
            )
            scan_rts = [float(scan.get("aligned_rt") or scan.get("rt") or 0.0) for scan in spectra]
            for edge in (item for item in candidate_edges if str(item["xic_sample_id"]) == sample_id):
                first = row_by_id[str(edge["first_id"])]
                second = row_by_id[str(edge["second_id"])]
                evidence = xic_evidence(
                    spectra,
                    scan_rts,
                    [
                        float(value)
                        for value in dict(refined_envelopes.get(str(edge["first_id"])) or {}).get(
                            "envelope_mz",
                            [float(first.get("representative_mz") or 0.0)],
                        )
                    ],
                    [
                        float(value)
                        for value in dict(refined_envelopes.get(str(edge["second_id"])) or {}).get(
                            "envelope_mz",
                            [float(second.get("representative_mz") or 0.0)],
                        )
                    ],
                    float(edge["xic_center_rt"]),
                )
                if evidence is None:
                    continue
                shape_score, apex_delta = evidence
                if shape_score < min_xic_shape_score or apex_delta > max_xic_apex_delta_min:
                    continue
                edge["xic_shape_score"] = shape_score
                edge["xic_apex_delta_min"] = apex_delta
                inference_score = (
                    0.30 * float(edge["mass_score"])
                    + 0.20 * float(edge["rt_score"])
                    + 0.30 * shape_score
                    + 0.20 * float(edge["trend_score"])
                )
                if inference_score < min_inference_score:
                    continue
                edge["inference_score"] = inference_score
                accepted_edges.append(edge)
    accepted_edges.sort(key=lambda edge: float(edge.get("inference_score") or 0.0), reverse=True)

    clusters: list[dict[str, object]] = []
    cluster_by_feature: dict[str, dict[str, object]] = {}

    def assignment_compatible(cluster: dict[str, object], edge: dict[str, object]) -> bool:
        assignments = dict(cluster["assignments"])
        for feature_id, assignment in dict(edge["assignments"]).items():
            previous = assignments.get(feature_id)
            if previous and (
                int(previous["charge"]) != int(assignment["charge"])
                or int(previous["isotope_offset"]) != int(assignment["isotope_offset"])
            ):
                return False
        neutral_mass = float(cluster["neutral_mass"])
        edge_mass = float(edge["neutral_mass"])
        return abs(neutral_mass - edge_mass) / max(neutral_mass, edge_mass, 1e-12) * 1_000_000.0 <= mass_tolerance_ppm

    for edge in accepted_edges:
        first_id = str(edge["first_id"])
        second_id = str(edge["second_id"])
        first_cluster = cluster_by_feature.get(first_id)
        second_cluster = cluster_by_feature.get(second_id)
        if first_cluster is not None and second_cluster is first_cluster:
            if assignment_compatible(first_cluster, edge):
                first_cluster["edges"].append(edge)
            continue
        if first_cluster is None and second_cluster is None:
            cluster = {
                "members": {first_id, second_id},
                "assignments": dict(edge["assignments"]),
                "neutral_mass": float(edge["neutral_mass"]),
                "edges": [edge],
            }
            clusters.append(cluster)
            cluster_by_feature[first_id] = cluster
            cluster_by_feature[second_id] = cluster
            continue
        target = first_cluster or second_cluster
        source = second_cluster if first_cluster is not None else first_cluster
        assert target is not None
        if len(set(target["members"]) | {first_id, second_id}) > 8:
            continue
        if not assignment_compatible(target, edge):
            continue
        if source is not None:
            if len(set(target["members"]) | set(source["members"])) > 8:
                continue
            if not assignment_compatible(source, edge):
                continue
            if any(
                feature_id in target["assignments"]
                and (
                    int(target["assignments"][feature_id]["charge"]) != int(assignment["charge"])
                    or int(target["assignments"][feature_id]["isotope_offset"]) != int(assignment["isotope_offset"])
                )
                for feature_id, assignment in dict(source["assignments"]).items()
            ):
                continue
            target["members"].update(source["members"])
            target["assignments"].update(source["assignments"])
            target["edges"].extend(source["edges"])
            target["neutral_mass"] = statistics.mean([float(target["neutral_mass"]), float(source["neutral_mass"])])
            for feature_id in source["members"]:
                cluster_by_feature[str(feature_id)] = target
            if source in clusters:
                clusters.remove(source)
        target["members"].update({first_id, second_id})
        target["assignments"].update(dict(edge["assignments"]))
        target["edges"].append(edge)
        target["neutral_mass"] = statistics.mean([float(target["neutral_mass"]), float(edge["neutral_mass"])])
        cluster_by_feature[first_id] = target
        cluster_by_feature[second_id] = target

    difference_order = {
        "presence_absence": 5,
        "area_changed": 4,
        "moderate_difference": 3,
        "common_feature": 2,
        "low_confidence": 1,
    }

    def summed_sample_map(entries: list[dict[str, object]], key: str) -> dict[str, float]:
        output: dict[str, float] = {}
        for feature in entries:
            for sample_id, value in dict(feature.get(key) or {}).items():
                output[str(sample_id)] = output.get(str(sample_id), 0.0) + max(0.0, float(value or 0.0))
        return output

    component_rows: list[dict[str, object]] = []
    inferred_feature_ids: set[str] = set()
    for cluster in clusters:
        member_ids = sorted(str(feature_id) for feature_id in cluster["members"])
        if len(member_ids) < 2:
            continue
        entries = [feature_rows[feature_id] for feature_id in member_ids]
        primary = max(entries, key=lambda row: float(row.get("ranking_score") or 0.0))
        relations_by_feature: dict[str, set[str]] = {feature_id: set() for feature_id in member_ids}
        for edge in cluster["edges"]:
            relation_type = str(edge["relation_type"])
            relations_by_feature[str(edge["first_id"])].add(relation_type)
            relations_by_feature[str(edge["second_id"])].add(relation_type)
        members: list[dict[str, object]] = []
        for feature in entries:
            feature_id = str(feature.get("feature_group_id"))
            assignment = dict(cluster["assignments"])[feature_id]
            relations = relations_by_feature[feature_id]
            envelope = refined_envelopes.get(feature_id) or {}
            member = dict(feature)
            member.update({
                "component_charge": int(assignment["charge"]),
                "component_isotope_offset": int(assignment["isotope_offset"]),
                "component_link_type": (
                    "ms1_inferred_isotope"
                    if "same_charge_isotope" in relations
                    else "ms1_inferred_charge"
                ),
                "component_confidence": "MS1_inferred",
                "true_peak_mz": feature.get("true_peak_mz", feature.get("representative_mz")),
                "component_observed_first_isotope_mz": envelope.get("observed_first_isotope_mz"),
                "component_monoisotopic_mz": envelope.get("monoisotopic_mz"),
                "component_monoisotopic_offset": envelope.get("monoisotopic_offset"),
                "envelope_representative_mz": (
                    envelope.get("monoisotopic_mz")
                    or envelope.get("observed_first_isotope_mz")
                    or feature.get("envelope_representative_mz")
                    or feature.get("representative_mz")
                ),
                "component_isotope_peak_count": envelope.get("isotope_peak_count"),
                "component_envelope_mz": envelope.get("envelope_mz") or [],
                "component_envelope_intensity": envelope.get("envelope_intensity") or [],
                "component_theoretical_isotope_intensity": envelope.get("theoretical_isotope_intensity") or [],
                "component_isotope_fit_score": envelope.get("isotope_fit_score"),
                "component_isotope_fit_cosine": envelope.get("isotope_fit_cosine"),
                "component_isotope_fit_explained_fraction": envelope.get("isotope_fit_explained_fraction"),
            })
            members.append(member)
            inferred_feature_ids.add(feature_id)
        members.sort(key=lambda member: (
            int(member.get("component_charge") or 0),
            int(member.get("component_isotope_offset") or 0),
            float(member.get("representative_mz") or 0.0),
        ))
        raw_area = summed_sample_map(entries, "area_by_sample")
        normalized_area = summed_sample_map(entries, "normalized_area_by_sample")
        relation_types = sorted({str(edge["relation_type"]) for edge in cluster["edges"]})
        mean_shape = statistics.mean(float(edge.get("xic_shape_score") or 0.0) for edge in cluster["edges"])
        mean_score = statistics.mean(float(edge.get("inference_score") or 0.0) for edge in cluster["edges"])
        isotope_fit_scores = [
            float(member.get("component_isotope_fit_score") or 0.0)
            for member in members
            if member.get("component_isotope_fit_score") is not None
        ]
        mean_isotope_fit = statistics.mean(isotope_fit_scores) if isotope_fit_scores else None
        component_mass_evidence = {
            str(edge.get("mass_evidence") or "")
            for edge in cluster["edges"]
            if edge.get("mass_evidence")
        }
        if isotope_fit_scores:
            component_mass_evidence.add("averagine_full_isotope_envelope")
        high_confidence = (
            mean_shape >= 0.90
            and mean_score >= 0.86
            and (len(member_ids) >= 3 or len(relation_types) >= 2)
        )
        confidence = "MS1_high" if high_confidence else "MS1_medium"
        component = dict(primary)
        if raw_area:
            component["area_by_sample"] = raw_area
            component["max_area"] = max(raw_area.values(), default=0.0)
        if normalized_area:
            component["normalized_area_by_sample"] = normalized_area
        neutral_mass = float(cluster["neutral_mass"])
        component.update({
            "identified_component": False,
            "inferred_component": True,
            "component_label": f"Unknown {neutral_mass:.4f} Da",
            "component_confidence": confidence,
            "component_member_count": len(members),
            "component_charge_states": sorted({int(member["component_charge"]) for member in members}),
            "component_isotope_offsets": sorted({int(member["component_isotope_offset"]) for member in members}),
            "component_neutral_mass": neutral_mass,
            "component_primary_feature_id": primary.get("feature_group_id"),
            "component_relation_types": relation_types,
            "component_mass_evidence": sorted(component_mass_evidence),
            "component_inference_score": mean_score,
            "component_xic_shape_score": mean_shape,
            "component_isotope_fit_score": mean_isotope_fit,
            "component_isotope_fit_min_score": min(isotope_fit_scores) if isotope_fit_scores else None,
            "component_monoisotope_corrected": any(
                int(member.get("component_monoisotopic_offset") or 0) > 0
                for member in members
            ),
            "members": members,
            "ranking_score": max(float(feature.get("ranking_score") or 0.0) for feature in entries),
            "difference_type": max(
                (str(feature.get("difference_type") or "") for feature in entries),
                key=lambda value: difference_order.get(value, 0),
            ),
            "merged_feature_count": sum(int(feature.get("merged_feature_count") or 1) for feature in entries),
        })
        primary_member = next(
            (
                member for member in members
                if str(member.get("feature_group_id")) == str(primary.get("feature_group_id"))
            ),
            members[0],
        )
        component["true_peak_mz"] = primary_member.get("true_peak_mz", primary_member.get("representative_mz"))
        component["envelope_representative_mz"] = primary_member.get(
            "envelope_representative_mz",
            primary_member.get("representative_mz"),
        )
        component["component_representative_mz"] = component["envelope_representative_mz"]
        component["component_representative_charge"] = primary_member.get("component_charge")
        component_rows.append(component)
    return component_rows, inferred_feature_ids


def build_ms1_component_groups(
    payload: dict[str, object],
    annotations: list[dict[str, object]],
    rt_tolerance_min: float = 0.75,
    include_low_evidence: bool = False,
    spectra_loader: Callable[[str], list[dict[str, object]]] | None = None,
) -> list[dict[str, object]]:
    """Collapse identified charge states and isotopes into component-level rows.

    A peptide form is keyed by base peptide plus neutral mass and is additionally
    split by retention time. This avoids merging separate chromatographic species
    that happen to receive the same sequence assignment. Every MS1 Feature is
    assigned to at most one component; ambiguous assignments retain the strongest
    sequence evidence. Unassigned Features remain as singleton rows.
    """
    feature_rows = {
        str(row.get("feature_group_id") or ""): row
        for row in payload.get("global_feature_groups") or []
        if row.get("feature_group_id")
    }
    confidence_order = {"B_high_confidence_inferred": 3, "C_tentative": 2, "D_low_evidence": 1}
    link_order = {
        "direct_precursor": 6,
        "charge_state_envelope": 5,
        "isotope_envelope": 4,
        "component_charge_isotope_consensus_search": 3,
        "component_mass_offset_sequence_tag_search": 3,
        "selected_precursor_consensus_search": 3,
        "feature_targeted_glycopeptide_search": 3,
        "global_open_modification_search": 3,
        "isolation_window_targeted_search": 2,
        "unlinked_sequence": 0,
    }

    def annotation_score(annotation: dict[str, object]) -> tuple[float, ...]:
        sample_scores = [
            float(item.get("score") or 0.0)
            for item in dict(annotation.get("sample_evidence") or {}).values()
            if isinstance(item, dict)
        ]
        return (
            float(confidence_order.get(str(annotation.get("confidence") or ""), 0)),
            float(link_order.get(str(annotation.get("feature_link_type") or ""), 0)),
            max(sample_scores, default=0.0),
            -abs(float(annotation.get("feature_mz_error_ppm") or 0.0)),
            float(annotation.get("feature_ranking") or 0.0),
        )

    assignments: dict[str, dict[str, object]] = {}
    low_evidence_by_feature: dict[str, dict[str, object]] = {}
    for annotation in annotations:
        feature_id = str(annotation.get("feature_or_candidate_id") or "")
        confidence = str(annotation.get("confidence") or "")
        if feature_id not in feature_rows or confidence not in confidence_order:
            continue
        if confidence == "D_low_evidence" and not include_low_evidence:
            previous = low_evidence_by_feature.get(feature_id)
            if previous is None or annotation_score(annotation) > annotation_score(previous):
                low_evidence_by_feature[feature_id] = annotation
            continue
        previous = assignments.get(feature_id)
        if previous is None or annotation_score(annotation) > annotation_score(previous):
            assignments[feature_id] = annotation

    identity_buckets: dict[tuple[str, object], list[tuple[dict[str, object], dict[str, object]]]] = {}
    for feature_id, annotation in assignments.items():
        neutral_mass = annotation.get("neutral_mass")
        if neutral_mass is not None and float(neutral_mass) > 0.0:
            form_identity: object = round(float(neutral_mass), 4)
        else:
            form_identity = str(annotation.get("candidate_id") or "")
        identity = (str(annotation.get("base_peptide_id") or annotation.get("sequence") or ""), form_identity)
        identity_buckets.setdefault(identity, []).append((feature_rows[feature_id], annotation))

    clustered: list[list[tuple[dict[str, object], dict[str, object]]]] = []
    for entries in identity_buckets.values():
        entries.sort(key=lambda item: float(item[0].get("representative_rt") or 0.0))
        local_clusters: list[list[tuple[dict[str, object], dict[str, object]]]] = []
        for entry in entries:
            rt = float(entry[0].get("representative_rt") or 0.0)
            if not local_clusters:
                local_clusters.append([entry])
                continue
            center = statistics.mean(float(item[0].get("representative_rt") or 0.0) for item in local_clusters[-1])
            if abs(rt - center) <= max(0.0, rt_tolerance_min):
                local_clusters[-1].append(entry)
            else:
                local_clusters.append([entry])
        clustered.extend(local_clusters)

    def summed_sample_map(entries: list[tuple[dict[str, object], dict[str, object]]], key: str) -> dict[str, float]:
        output: dict[str, float] = {}
        for feature, _ in entries:
            for sample_id, value in dict(feature.get(key) or {}).items():
                output[str(sample_id)] = output.get(str(sample_id), 0.0) + max(0.0, float(value or 0.0))
        return output

    def primary_score(entry: tuple[dict[str, object], dict[str, object]]) -> tuple[float, ...]:
        feature, annotation = entry
        link_type = str(annotation.get("feature_link_type") or "")
        isotope_offset = int(annotation.get("feature_isotope_offset") or 0)
        return (
            float(link_type == "direct_precursor"),
            float(isotope_offset == 0),
            float(confidence_order.get(str(annotation.get("confidence") or ""), 0)),
            float(link_order.get(link_type, 0)),
            float(feature.get("ranking_score") or 0.0),
        )

    difference_order = {
        "presence_absence": 5,
        "area_changed": 4,
        "moderate_difference": 3,
        "common_feature": 2,
        "low_confidence": 1,
    }
    component_rows: list[dict[str, object]] = []
    assigned_feature_ids: set[str] = set()
    for entries in clustered:
        primary_feature, primary_annotation = max(entries, key=primary_score)
        best_annotation = max((annotation for _, annotation in entries), key=annotation_score)
        raw_area = summed_sample_map(entries, "area_by_sample")
        normalized_area = summed_sample_map(entries, "normalized_area_by_sample")
        members: list[dict[str, object]] = []
        for feature, annotation in sorted(
            entries,
            key=lambda item: (
                int(item[1].get("feature_charge") or 0),
                int(item[1].get("feature_isotope_offset") or 0),
                float(item[0].get("representative_mz") or 0.0),
            ),
        ):
            feature_id = str(feature.get("feature_group_id") or "")
            assigned_feature_ids.add(feature_id)
            member = dict(feature)
            feature_charge = int(annotation.get("feature_charge") or 0)
            neutral_mass = float(annotation.get("neutral_mass") or 0.0)
            envelope_mz = (
                (neutral_mass / feature_charge + PROTON)
                if feature_charge > 0 and neutral_mass > 0.0
                else float(feature.get("representative_mz") or 0.0)
            )
            member.update({
                "component_charge": annotation.get("feature_charge"),
                "component_isotope_offset": int(annotation.get("feature_isotope_offset") or 0),
                "component_link_type": annotation.get("feature_link_type"),
                "component_confidence": annotation.get("confidence"),
                "component_candidate_id": annotation.get("candidate_id"),
                "true_peak_mz": feature.get("true_peak_mz", feature.get("representative_mz")),
                "envelope_representative_mz": envelope_mz,
            })
            members.append(member)
        modifications = sorted({str(annotation.get("modification") or "Unmodified") for _, annotation in entries})
        component_modification = format_modification_alternatives(
            modifications,
            str(best_annotation.get("sequence") or ""),
        )
        component = dict(primary_feature)
        if raw_area:
            component["area_by_sample"] = raw_area
            component["max_area"] = max(raw_area.values(), default=0.0)
        if normalized_area:
            component["normalized_area_by_sample"] = normalized_area
        component.update({
            "identified_component": True,
            "component_label": f"{best_annotation.get('sequence') or ''} | {component_modification}",
            "component_sequence": best_annotation.get("sequence"),
            "component_modification": component_modification,
            "component_chain": best_annotation.get("chain"),
            "component_start": best_annotation.get("start"),
            "component_end": best_annotation.get("end"),
            "component_confidence": best_annotation.get("confidence"),
            "component_member_count": len(members),
            "component_charge_states": sorted({int(member["component_charge"]) for member in members if member.get("component_charge")}),
            "component_isotope_offsets": sorted({int(member["component_isotope_offset"]) for member in members}),
            "component_candidate_ids": sorted({str(annotation.get("candidate_id") or "") for _, annotation in entries}),
            "component_neutral_mass": best_annotation.get("neutral_mass"),
            "component_primary_feature_id": primary_feature.get("feature_group_id"),
            "members": members,
            "ranking_score": max(float(feature.get("ranking_score") or 0.0) for feature, _ in entries),
            "difference_type": max(
                (str(feature.get("difference_type") or "") for feature, _ in entries),
                key=lambda value: difference_order.get(value, 0),
            ),
            "merged_feature_count": sum(int(feature.get("merged_feature_count") or 1) for feature, _ in entries),
        })
        primary_member = next(
            (
                member for member in members
                if str(member.get("feature_group_id")) == str(primary_feature.get("feature_group_id"))
            ),
            members[0],
        )
        component["true_peak_mz"] = primary_member.get("true_peak_mz", primary_member.get("representative_mz"))
        component["envelope_representative_mz"] = primary_member.get(
            "envelope_representative_mz",
            primary_member.get("representative_mz"),
        )
        component["component_representative_mz"] = component["envelope_representative_mz"]
        component["component_representative_charge"] = primary_member.get("component_charge")
        component_rows.append(component)

    inferred_components, inferred_feature_ids = _infer_unidentified_ms1_groups(
        payload,
        feature_rows,
        assigned_feature_ids,
        spectra_loader,
    )
    component_rows.extend(inferred_components)
    assigned_feature_ids.update(inferred_feature_ids)

    for feature_id, feature in feature_rows.items():
        if feature_id in assigned_feature_ids:
            continue
        singleton = dict(feature)
        singleton.update({
            "identified_component": False,
            "component_label": feature_id,
            "component_member_count": 1,
            "component_primary_feature_id": feature_id,
            "true_peak_mz": feature.get("true_peak_mz", feature.get("representative_mz")),
            "envelope_representative_mz": feature.get(
                "envelope_representative_mz",
                feature.get("representative_mz"),
            ),
            "component_representative_mz": feature.get(
                "envelope_representative_mz",
                feature.get("representative_mz"),
            ),
            "members": [dict(feature)],
        })
        component_rows.append(singleton)

    # A D-level match is useful review evidence, but it must not turn an MS1
    # component into an identified peptide.  Attach the best candidate to the
    # unresolved/inferred component so every consumer can show one consistent
    # label while preserving the Unknown identity and confidence boundary.
    for row in component_rows:
        if row.get("identified_component"):
            continue
        feature_ids = {
            str(row.get("feature_group_id") or ""),
            str(row.get("component_primary_feature_id") or ""),
        }
        feature_ids.update(str(value) for value in row.get("source_feature_group_ids") or [])
        feature_ids.update(str(value) for value in row.get("merged_feature_group_ids") or [])
        for member in row.get("members") or []:
            if not isinstance(member, dict):
                continue
            feature_ids.add(str(member.get("feature_group_id") or ""))
            feature_ids.update(str(value) for value in member.get("source_feature_group_ids") or [])
            feature_ids.update(str(value) for value in member.get("merged_feature_group_ids") or [])
        candidates = [
            low_evidence_by_feature[feature_id]
            for feature_id in feature_ids
            if feature_id in low_evidence_by_feature
        ]
        if not candidates:
            continue
        candidate = max(candidates, key=annotation_score)
        sequence = str(candidate.get("sequence") or "").strip()
        modification = str(candidate.get("modification") or "Unmodified").strip()
        base_label = str(row.get("component_label") or row.get("feature_group_id") or "Unknown")
        candidate_text = f"candidate {sequence}" if sequence else "low-evidence sequence candidate"
        if sequence and sequence not in base_label:
            row["component_label"] = f"{base_label} | {candidate_text}"
        row.update({
            "candidate_component": True,
            "component_candidate_sequence": sequence or None,
            "component_candidate_modification": modification or None,
            "component_candidate_confidence": "D_low_evidence",
            "component_candidate_feature_ids": sorted(
                feature_id for feature_id in feature_ids
                if feature_id in low_evidence_by_feature
            ),
        })

    component_rows.sort(key=lambda row: float(row.get("ranking_score") or 0.0), reverse=True)
    identified_index = 0
    inferred_index = 0
    for row in component_rows:
        if row.get("identified_component"):
            identified_index += 1
            row["component_group_id"] = f"MS1COMP_{identified_index:04d}"
        elif row.get("inferred_component"):
            inferred_index += 1
            row["component_group_id"] = f"MS1INFER_{inferred_index:04d}"
        else:
            row["component_group_id"] = f"MS1SINGLE_{row.get('component_primary_feature_id')}"
    return component_rows


def build_modified_peptide_findings(
    psms: list[dict[str, object]],
    annotations: list[dict[str, object]],
    payload: dict[str, object] | None = None,
    peptide_form_comparisons: list[dict[str, object]] | None = None,
) -> list[dict[str, object]]:
    """Summarize every identified modified peptide, even without an unmodified pair.

    Isobaric site alternatives for the same peptide and precursor mass are kept
    together. PSM counts and precursor intensities are supporting acquisition
    evidence only; linked normalized MS1 areas remain the preferred direction
    evidence when a significant Feature is available.
    """
    feature_rows = {
        str(row.get("feature_group_id") or ""): row
        for row in (payload or {}).get("global_feature_groups") or []
    }
    significant_feature_ids = {
        feature_id for feature_id, row in feature_rows.items()
        if row.get("difference_type") in {"presence_absence", "area_changed", "moderate_difference"}
    }
    paired_candidate_ids = {
        str(candidate_id)
        for comparison in peptide_form_comparisons or []
        for candidate_id in comparison.get("modified_candidate_ids") or [comparison.get("modified_candidate_id")]
        if candidate_id
    }
    unmodified_base_ids = {
        str(psm.get("base_peptide_id") or "")
        for psm in psms
        if str(psm.get("modification_text") or "Unmodified") == "Unmodified"
    }
    annotation_features: dict[str, set[str]] = {}
    for annotation in annotations:
        feature_id = str(annotation.get("feature_or_candidate_id") or "")
        if feature_id in significant_feature_ids:
            annotation_features.setdefault(str(annotation.get("candidate_id") or ""), set()).add(feature_id)

    grouped: dict[tuple[str, float], list[dict[str, object]]] = {}
    for psm in psms:
        if psm.get("is_decoy") or str(psm.get("modification_text") or "Unmodified") == "Unmodified":
            continue
        grouped_mass = psm.get("global_delta_cluster_mass")
        if grouped_mass is not None:
            grouped_mass = float(psm.get("backbone_neutral_mass") or 0.0) + float(grouped_mass)
        else:
            grouped_mass = float(psm.get("neutral_mass") or 0.0)
        key = (str(psm.get("base_peptide_id") or ""), round(grouped_mass, 5))
        grouped.setdefault(key, []).append(psm)

    confidence_order = {"B_high_confidence_inferred": 3, "C_tentative": 2, "D_low_evidence": 1}
    findings: list[dict[str, object]] = []
    for (base_peptide_id, neutral_mass), evidence in grouped.items():
        best = max(
            evidence,
            key=lambda row: (
                confidence_order[psm_confidence(row)],
                float(row.get("score") or 0.0),
                -float(row["q_value"] if row.get("q_value") is not None else 1.0),
            ),
        )
        candidate_ids = sorted({str(row.get("candidate_id") or "") for row in evidence})
        alternatives = sorted({str(row.get("modification_text") or "") for row in evidence})
        linked_feature_ids = {
            str(link.get("feature_group_id") or "")
            for row in evidence
            for link in row.get("feature_links") or []
            if str(link.get("feature_group_id") or "") in significant_feature_ids
        }
        for candidate_id in candidate_ids:
            linked_feature_ids.update(annotation_features.get(candidate_id) or set())
        linked_features = sorted(
            (feature_rows[feature_id] for feature_id in linked_feature_ids),
            key=lambda row: float(row.get("ranking_score") or 0.0),
            reverse=True,
        )
        leading_feature = linked_features[0] if linked_features else {}
        sample_evidence: dict[str, dict[str, object]] = {}
        for sample_id in sorted({str(row.get("sample_id") or "") for row in evidence}):
            rows = [row for row in evidence if str(row.get("sample_id") or "") == sample_id]
            sample_best = max(rows, key=lambda row: float(row.get("score") or 0.0))
            sample_evidence[sample_id] = {
                "psm_count": len(rows),
                "best_score": sample_best.get("score"),
                "best_q_value": min(float(row["q_value"] if row.get("q_value") is not None else 1.0) for row in rows),
                "best_scan_id": sample_best.get("scan_id"),
                "median_rt": statistics.median(float(row.get("rt") or 0.0) for row in rows),
                "charge_states": sorted({int(row.get("precursor_charge") or 0) for row in rows if int(row.get("precursor_charge") or 0) > 0}),
                "summed_precursor_intensity": sum(float(row.get("precursor_intensity") or 0.0) for row in rows),
            }
        is_paired = any(candidate_id in paired_candidate_ids for candidate_id in candidate_ids)
        has_unmodified_ms2 = base_peptide_id in unmodified_base_ids
        if is_paired:
            pairing_status = "paired_with_unmodified_form"
        elif has_unmodified_ms2:
            pairing_status = "unpaired_but_unmodified_ms2_seen"
        else:
            pairing_status = "modified_identified_without_unmodified_ms2"
        modifications = best.get("modifications") or []
        findings.append({
            "finding_id": "",
            "base_peptide_id": base_peptide_id,
            "candidate_ids": candidate_ids,
            "chain": best.get("chain"),
            "start": best.get("start"),
            "end": best.get("end"),
            "sequence": best.get("sequence"),
            "modification": format_modification_alternatives(
                alternatives,
                str(best.get("sequence") or ""),
            ),
            "modification_alternatives": alternatives,
            "modification_mass_delta": (
                sum(float(item[2]) for item in modifications)
                if modifications
                else float(best.get("mass_delta") or 0.0)
            ),
            "neutral_mass": neutral_mass,
            "proteolysis": best.get("proteolysis", "fully_tryptic"),
            "confidence": psm_confidence(best),
            "pairing_status": pairing_status,
            "evidence_scope": "differential_feature_linked" if linked_feature_ids else "ms2_identified_only",
            "linked_feature_ids": sorted(linked_feature_ids),
            "feature_ranking": leading_feature.get("ranking_score"),
            "difference_type": leading_feature.get("difference_type"),
            "max_fold_change": leading_feature.get("max_fold_change"),
            "higher_abundance_sample": leading_feature.get("higher_abundance_sample"),
            "lower_abundance_sample": leading_feature.get("lower_abundance_sample"),
            "normalized_area_by_sample": leading_feature.get("normalized_area_by_sample") or {},
            "sample_evidence": sample_evidence,
            "best_psm": {
                key: best.get(key) for key in (
                    "sample_id", "scan_id", "rt", "precursor_mz", "precursor_charge",
                    "chain", "start", "end", "sequence", "modification_text", "proteolysis",
                    "score", "q_value", "matched_ion_count", "fragment_coverage",
                    "explained_intensity", "spectrum_peaks", "mass_delta",
                    "localized_mass_offset_site", "global_delta_cluster_id",
                    "global_delta_cluster_mass", "global_delta_cluster_q_value",
                    "global_delta_cluster_support", "open_modification_confidence",
                )
            },
        })
    findings.sort(
        key=lambda row: (
            row["evidence_scope"] == "differential_feature_linked",
            confidence_order.get(str(row.get("confidence") or ""), 0),
            float(row.get("feature_ranking") or 0.0),
            max((float(item.get("best_score") or 0.0) for item in dict(row.get("sample_evidence") or {}).values()), default=0.0),
        ),
        reverse=True,
    )
    for rank, finding in enumerate(findings, start=1):
        finding["rank"] = rank
        finding["finding_id"] = f"MOD_{rank:04d}"
    return findings


def _feature_link_payload(row: dict[str, object], link_type: str, isotope_offset: int = 0) -> dict[str, object]:
    return {
        "feature_group_id": row.get("feature_group_id"),
        "feature_mz": row.get("representative_mz"),
        "feature_sample_mz": row.get("representative_mz"),
        "feature_rt": row.get("representative_rt"),
        "feature_sample_rt": row.get("representative_rt"),
        "feature_similarity": row.get("similarity_score"),
        "feature_difference": row.get("difference_score"),
        "feature_ranking": row.get("ranking_score"),
        "feature_difference_type": row.get("difference_type"),
        "feature_area_by_sample": row.get("area_by_sample") or {},
        "feature_normalized_area_by_sample": row.get("normalized_area_by_sample") or {},
        "feature_link_type": link_type,
        "feature_isotope_offset": isotope_offset,
        "feature_mz_error_da": 0.0,
        "feature_mz_error_ppm": 0.0,
        "feature_link_distance": 0.0,
    }


def _scan_isotope_relation(
    feature_mz: float,
    precursor_mz: float,
    charge: int | None,
    tolerance_da: float,
    min_isotope_offset: int = -1,
    max_isotope_offset: int = 4,
) -> tuple[int, float] | None:
    if not charge or charge <= 0:
        return (0, feature_mz - precursor_mz) if abs(feature_mz - precursor_mz) <= tolerance_da else None
    relations = [
        (offset, feature_mz - (precursor_mz + offset * ISOTOPE_MASS_DIFF / charge))
        for offset in range(min_isotope_offset, max_isotope_offset + 1)
    ]
    offset, error = min(relations, key=lambda item: abs(item[1]))
    return (offset, error) if abs(error) <= tolerance_da else None


def _top_spectrum_peaks(scan: LCMSSpectrumScan, limit: int = 60) -> list[dict[str, object]]:
    pairs = sorted(zip(scan.mz_array, scan.intensity_array), key=lambda pair: pair[0])
    top_indices = sorted(range(len(pairs)), key=lambda index: pairs[index][1], reverse=True)[:max(0, limit)]
    return [
        {"mz": pairs[index][0], "intensity": pairs[index][1], "label": ""}
        for index in sorted(top_indices, key=lambda index: pairs[index][0])
    ]


def _merge_consensus_ms2(
    scans: list[LCMSSpectrumScan],
    fragment_tolerance_ppm: float,
    max_points_per_scan: int = 300,
    max_consensus_points: int = 500,
    minimum_scan_support: int | None = None,
) -> tuple[list[float], list[float]]:
    points: list[tuple[float, float, int]] = []
    for scan_index, scan in enumerate(scans):
        if not scan.mz_array:
            continue
        selected = sorted(
            range(len(scan.mz_array)),
            key=scan.intensity_array.__getitem__,
            reverse=True,
        )[:max_points_per_scan]
        scale = max((float(scan.intensity_array[index]) for index in selected), default=0.0)
        if scale <= 0:
            continue
        points.extend(
            (float(scan.mz_array[index]), float(scan.intensity_array[index]) / scale, scan_index)
            for index in selected
        )
    clusters: list[dict[str, object]] = []
    for mz, intensity, scan_index in sorted(points):
        cluster = clusters[-1] if clusters else None
        center = float(cluster["weighted_mz"]) / max(float(cluster["weight"]), 1e-12) if cluster else 0.0
        tolerance = max(mz * fragment_tolerance_ppm / 1_000_000.0, 0.01)
        if cluster is None or abs(mz - center) > tolerance:
            cluster = {"weighted_mz": 0.0, "weight": 0.0, "intensity": 0.0, "scans": set()}
            clusters.append(cluster)
        cluster["weighted_mz"] = float(cluster["weighted_mz"]) + mz * intensity
        cluster["weight"] = float(cluster["weight"]) + intensity
        cluster["intensity"] = float(cluster["intensity"]) + intensity
        scan_support = cluster["scans"]
        assert isinstance(scan_support, set)
        scan_support.add(scan_index)
    minimum_support = (
        max(1, int(minimum_scan_support))
        if minimum_scan_support is not None
        else (2 if len(scans) >= 2 else 1)
    )
    supported = [cluster for cluster in clusters if len(cluster["scans"]) >= minimum_support]
    selected_clusters = sorted(
        supported,
        key=lambda cluster: float(cluster["intensity"]),
        reverse=True,
    )[:max_consensus_points]
    consensus = sorted(
        (
            float(cluster["weighted_mz"]) / max(float(cluster["weight"]), 1e-12),
            float(cluster["intensity"]) * 1_000_000.0,
        )
        for cluster in selected_clusters
    )
    return [mz for mz, _ in consensus], [intensity for _, intensity in consensus]


def search_selected_feature_consensus_scans(
    payload: dict[str, object],
    scans: list[LCMSSpectrumScan],
    candidates: list[PeptideCandidate],
    exclude_feature_ids: set[str] | None = None,
    precursor_tolerance_ppm: float = 20.0,
    fragment_tolerance_ppm: float = 20.0,
    rt_tolerance_min: float = 0.5,
    isolation_padding_da: float = 0.02,
    max_scans_per_sample: int = 8,
    fdr_threshold: float = 0.01,
) -> list[dict[str, object]]:
    """Merge repeated selected-precursor MS2 scans for unresolved MS1 Features."""
    excluded = exclude_feature_ids or set()
    significant_types = {"presence_absence", "area_changed", "moderate_difference"}
    features = [
        row for row in payload.get("global_feature_groups") or []
        if row.get("difference_type") in significant_types and str(row.get("feature_group_id")) not in excluded
    ]
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
    for scan in scans:
        if scan.ms_level == 2 and scan.precursor_mz is not None and scan.mz_array:
            scans_by_sample.setdefault(scan.sample_id, []).append(scan)
    synthetic_scans: list[LCMSSpectrumScan] = []
    metadata: dict[str, tuple[dict[str, object], list[LCMSSpectrumScan]]] = {}
    for feature in features:
        feature_id = str(feature.get("feature_group_id") or "")
        feature_mz = float(feature.get("representative_mz") or 0.0)
        feature_rt = float(feature.get("representative_rt") or 0.0)
        local_shifts = dict(feature.get("rt_correction_by_sample") or {})
        for sample_id, sample_scans in scans_by_sample.items():
            selected: list[tuple[float, float, LCMSSpectrumScan]] = []
            for scan in sample_scans:
                aligned_rt = scan.rt + float(shifts.get(sample_id) or 0.0) + float(local_shifts.get(sample_id) or 0.0)
                rt_error = aligned_rt - feature_rt
                if abs(rt_error) > rt_tolerance_min:
                    continue
                precursor_mz = float(scan.precursor_mz or 0.0)
                lower = float(scan.isolation_window_lower_offset or 0.8)
                upper = float(scan.isolation_window_upper_offset or 0.8)
                if not precursor_mz - lower - isolation_padding_da <= feature_mz <= precursor_mz + upper + isolation_padding_da:
                    continue
                tolerance = max(feature_mz * precursor_tolerance_ppm / 1_000_000.0, isolation_padding_da)
                if _scan_isotope_relation(feature_mz, precursor_mz, scan.precursor_charge, tolerance) is None:
                    continue
                selected.append((abs(rt_error), -float(scan.precursor_intensity or 0.0), scan))
            source_scans = [item[2] for item in sorted(selected)[:max_scans_per_sample]]
            if len(source_scans) < 2:
                continue
            mz_array, intensity_array = _merge_consensus_ms2(source_scans, fragment_tolerance_ppm)
            if not mz_array:
                continue
            charges = [int(scan.precursor_charge) for scan in source_scans if scan.precursor_charge and scan.precursor_charge > 0]
            charge = max(set(charges), key=charges.count) if charges else None
            precursor_values = [float(scan.precursor_mz or 0.0) for scan in source_scans]
            precursor_mz = statistics.median(precursor_values)
            synthetic_id = f"consensus_feature={feature_id}|sample={sample_id}"
            synthetic = LCMSSpectrumScan(
                scan_id=synthetic_id,
                raw_file_id="consensus",
                sample_id=sample_id,
                rt=statistics.median(scan.rt for scan in source_scans),
                ms_level=2,
                mz_array=mz_array,
                intensity_array=intensity_array,
                tic=sum(intensity_array),
                base_peak_mz=mz_array[max(range(len(mz_array)), key=intensity_array.__getitem__)],
                base_peak_intensity=max(intensity_array),
                precursor_mz=precursor_mz,
                precursor_charge=charge,
                precursor_intensity=max(float(scan.precursor_intensity or 0.0) for scan in source_scans),
                activation_method=source_scans[0].activation_method,
                collision_energy=source_scans[0].collision_energy,
            )
            synthetic_scans.append(synthetic)
            metadata[synthetic_id] = (feature, source_scans)
    if not synthetic_scans:
        return []
    searched = assign_q_values(search_scans(
        synthetic_scans,
        candidates,
        precursor_tolerance_ppm=precursor_tolerance_ppm,
        fragment_tolerance_ppm=fragment_tolerance_ppm,
        min_score=20.0,
    ))
    accepted: list[dict[str, object]] = []
    for psm in searched:
        q_value = float(psm["q_value"] if psm.get("q_value") is not None else 1.0)
        if psm.get("is_decoy") or q_value > fdr_threshold:
            continue
        if float(psm.get("score") or 0.0) < 50.0 or int(psm.get("matched_ion_count") or 0) < 6:
            continue
        feature, source_scans = metadata[str(psm["scan_id"])]
        link = _feature_link_payload(feature, "selected_precursor_consensus_search")
        psm.update({
            "search_origin": "selected_precursor_consensus_search",
            "consensus_source_scan_ids": [scan.scan_id for scan in source_scans],
            "consensus_source_scan_count": len(source_scans),
            "feature_links": [link],
            **{key: value for key, value in link.items() if key != "feature_link_distance"},
        })
        accepted.append(psm)
    return accepted


def search_component_consensus_scans(
    payload: dict[str, object],
    scans: list[LCMSSpectrumScan],
    candidates: list[PeptideCandidate],
    component_groups: list[dict[str, object]],
    exclude_feature_ids: set[str] | None = None,
    component_mass_tolerance_ppm: float = 20.0,
    fragment_tolerance_ppm: float = 20.0,
    rt_tolerance_min: float = 0.5,
    max_scans_per_member_sample: int = 4,
    fdr_threshold: float = 0.01,
    min_score: float = 45.0,
    min_matched_ions: int = 5,
    min_fragment_coverage: float = 0.15,
) -> list[dict[str, object]]:
    """Jointly search MS2 evidence across charge/isotope members of one MS1 component.

    The component neutral mass limits candidate sequences before fragment matching.
    Each acquisition scan is assigned to only its closest component member, while
    b/y evidence is combined at the ion level across repeated isotope and precursor
    charge-state scans. Isolation-window-only scans are deliberately excluded.
    """
    excluded = {str(value) for value in (exclude_feature_ids or set())}
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
    for scan in scans:
        if scan.ms_level == 2 and scan.precursor_mz is not None and scan.mz_array:
            scans_by_sample.setdefault(str(scan.sample_id), []).append(scan)
    for sample_scans in scans_by_sample.values():
        sample_scans.sort(key=lambda scan: float(scan.rt))
    for sample_scans in scans_by_sample.values():
        sample_scans.sort(key=lambda scan: float(scan.rt))

    ordered_candidates = sorted(candidates, key=lambda candidate: candidate.neutral_mass)
    candidate_masses = [candidate.neutral_mass for candidate in ordered_candidates]
    winners: list[dict[str, object]] = []

    def component_candidate_evidence(
        source_records: list[dict[str, object]],
        candidate: PeptideCandidate,
        component_mass: float,
    ) -> dict[str, object]:
        mass_error_ppm = (
            (component_mass - candidate.neutral_mass)
            / max(candidate.neutral_mass, 1e-12)
            * 1_000_000.0
        )
        position_set: set[tuple[str, int]] = set()
        label_matches: dict[str, dict[str, object]] = {}
        informative_records: list[dict[str, object]] = []
        explained_values: list[float] = []
        for record in source_records:
            scan = record["scan"]
            assert isinstance(scan, LCMSSpectrumScan)
            evidence = match_fragments(
                scan,
                candidate,
                mass_error_ppm,
                component_mass_tolerance_ppm,
                fragment_tolerance_ppm,
            )
            scan_positions = {
                (str(match["series"]), int(match["ordinal"]))
                for match in evidence["matched_fragments"]
            }
            position_set.update(scan_positions)
            for match in evidence["matched_fragments"]:
                label = str(match["label"])
                previous = label_matches.get(label)
                if previous is None or float(match["intensity"]) > float(previous["intensity"]):
                    label_matches[label] = dict(match)
            if len(scan_positions) >= 2:
                informative_records.append(record)
                explained_values.append(float(evidence["explained_intensity"]))

        denominator = max(1, 2 * (len(candidate.sequence) - 1))
        coverage = len(position_set) / denominator
        continuity = max(
            _longest_run({ordinal for series, ordinal in position_set if series == "b"}),
            _longest_run({ordinal for series, ordinal in position_set if series == "y"}),
        ) / max(1, len(candidate.sequence) - 1)
        scan_support = len(informative_records)
        member_support_ids = {
            str(record["member_id"])
            for record in informative_records
        }
        charge_support = {
            int(record["member_charge"])
            for record in informative_records
            if int(record.get("member_charge") or 0) > 0
        }
        explained = (
            statistics.mean(sorted(explained_values, reverse=True)[:4])
            if explained_values else 0.0
        )
        mass_score = max(
            0.0,
            1.0
            - abs(mass_error_ppm)
            / max(component_mass_tolerance_ppm, 1e-9),
        )
        score = (
            25.0 * mass_score
            + 30.0 * min(1.0, coverage)
            + 15.0 * min(1.0, explained)
            + 10.0 * min(1.0, continuity)
            + 10.0 * min(1.0, scan_support / 2.0)
            + 5.0 * min(1.0, len(member_support_ids) / 2.0)
            + 5.0 * min(1.0, len(charge_support) / 2.0)
        )
        return {
            "score": score,
            "matched_ion_count": len(label_matches),
            "fragment_coverage": coverage,
            "explained_intensity": explained,
            "ion_continuity": continuity,
            "matched_fragments": list(label_matches.values()),
            "component_scan_support_count": scan_support,
            "component_member_support_count": len(member_support_ids),
            "component_charge_support_count": len(charge_support),
            "component_source_member_ids": sorted(member_support_ids),
            "component_source_charge_states": sorted(charge_support),
        }

    for component in component_groups:
        if not component.get("inferred_component"):
            continue
        component_mass = float(component.get("component_neutral_mass") or 0.0)
        if component_mass <= 0.0:
            continue
        members = [
            member for member in component.get("members") or []
            if str(member.get("feature_group_id") or "") not in excluded
        ]
        if len(members) < 2:
            continue
        tolerance = (
            component_mass
            * max(component_mass_tolerance_ppm, 0.1)
            / 1_000_000.0
        )
        left = bisect_left(candidate_masses, component_mass - tolerance)
        right = bisect_right(candidate_masses, component_mass + tolerance)
        local_candidates = ordered_candidates[left:right]
        if not local_candidates:
            continue

        for sample_id, sample_scans in scans_by_sample.items():
            best_assignment_by_scan: dict[str, tuple[tuple[float, ...], dict[str, object]]] = {}
            for member in members:
                member_id = str(member.get("feature_group_id") or "")
                member_mz = float(
                    member.get("true_peak_mz")
                    or member.get("representative_mz")
                    or 0.0
                )
                member_rt = float(member.get("representative_rt") or 0.0)
                member_charge = int(member.get("component_charge") or 0)
                local_shift = float(
                    dict(member.get("rt_correction_by_sample") or {}).get(sample_id)
                    or 0.0
                )
                for scan in sample_scans:
                    if (
                        member_charge > 0
                        and scan.precursor_charge
                        and int(scan.precursor_charge) != member_charge
                    ):
                        continue
                    aligned_rt = (
                        float(scan.rt)
                        + float(shifts.get(sample_id) or 0.0)
                        + local_shift
                    )
                    rt_error = aligned_rt - member_rt
                    if abs(rt_error) > max(0.0, rt_tolerance_min):
                        continue
                    precursor_mz = float(scan.precursor_mz or 0.0)
                    mz_tolerance = max(
                        member_mz
                        * max(component_mass_tolerance_ppm, 0.1)
                        / 1_000_000.0,
                        0.02,
                    )
                    relation = _scan_isotope_relation(
                        member_mz,
                        precursor_mz,
                        scan.precursor_charge,
                        mz_tolerance,
                    )
                    if relation is None:
                        continue
                    isotope_offset, mz_error = relation
                    key = (
                        float(abs(isotope_offset)),
                        abs(float(mz_error)) / max(mz_tolerance, 1e-12),
                        abs(rt_error) / max(rt_tolerance_min, 1e-12),
                        -float(scan.precursor_intensity or 0.0),
                    )
                    record = {
                        "scan": scan,
                        "member_id": member_id,
                        "member_charge": member_charge or int(scan.precursor_charge or 0),
                        "member_isotope_offset": int(
                            member.get("component_isotope_offset") or 0
                        ),
                    }
                    previous = best_assignment_by_scan.get(str(scan.scan_id))
                    if previous is None or key < previous[0]:
                        best_assignment_by_scan[str(scan.scan_id)] = (key, record)

            records_by_member: dict[str, list[tuple[tuple[float, ...], dict[str, object]]]] = {}
            for key, record in best_assignment_by_scan.values():
                records_by_member.setdefault(str(record["member_id"]), []).append((key, record))
            source_records = [
                record
                for member_records in records_by_member.values()
                for _, record in sorted(member_records, key=lambda item: item[0])[
                    :max(1, int(max_scans_per_member_sample))
                ]
            ]
            if len(source_records) < 2 or len({str(record["member_id"]) for record in source_records}) < 2:
                continue

            candidate_rows: list[dict[str, object]] = []
            for candidate in local_candidates:
                joint = component_candidate_evidence(
                    source_records,
                    candidate,
                    component_mass,
                )
                if float(joint["score"]) < 20.0:
                    continue
                representative_charge = int(
                    component.get("component_representative_charge")
                    or next(
                        (
                            record["member_charge"] for record in source_records
                            if int(record.get("member_charge") or 0) > 0
                        ),
                        2,
                    )
                )
                source_scans = [
                    record["scan"] for record in source_records
                    if isinstance(record.get("scan"), LCMSSpectrumScan)
                ]
                mz_array, intensity_array = _merge_consensus_ms2(
                    source_scans,
                    fragment_tolerance_ppm,
                    minimum_scan_support=1,
                )
                synthetic = LCMSSpectrumScan(
                    scan_id=(
                        f"component={component.get('component_group_id')}"
                        f"|sample={sample_id}"
                    ),
                    raw_file_id="component_consensus",
                    sample_id=sample_id,
                    rt=statistics.median(float(scan.rt) for scan in source_scans),
                    ms_level=2,
                    mz_array=mz_array,
                    intensity_array=intensity_array,
                    tic=sum(intensity_array),
                    base_peak_mz=(
                        mz_array[max(range(len(mz_array)), key=intensity_array.__getitem__)]
                        if mz_array else None
                    ),
                    base_peak_intensity=max(intensity_array, default=0.0),
                    precursor_mz=component_mass / representative_charge + PROTON,
                    precursor_charge=representative_charge,
                    precursor_intensity=max(
                        (float(scan.precursor_intensity or 0.0) for scan in source_scans),
                        default=0.0,
                    ),
                    activation_method=source_scans[0].activation_method,
                    collision_energy=source_scans[0].collision_energy,
                )
                display_evidence = match_fragments(
                    synthetic,
                    candidate,
                    (
                        (component_mass - candidate.neutral_mass)
                        / max(candidate.neutral_mass, 1e-12)
                        * 1_000_000.0
                    ),
                    component_mass_tolerance_ppm,
                    fragment_tolerance_ppm,
                )
                display_pairs = sorted(
                    zip(synthetic.mz_array, synthetic.intensity_array),
                    key=lambda pair: pair[0],
                )
                labels = {
                    int(match["observed_index"]): str(match["label"])
                    for match in display_evidence["matched_fragments"]
                }
                top_indices = sorted(
                    range(len(display_pairs)),
                    key=lambda index: display_pairs[index][1],
                    reverse=True,
                )[:60]
                spectrum_peaks = [
                    {
                        "mz": display_pairs[index][0],
                        "intensity": display_pairs[index][1],
                        "label": labels.get(index, ""),
                    }
                    for index in sorted(top_indices, key=lambda index: display_pairs[index][0])
                ]
                row = {
                    "sample_id": sample_id,
                    "scan_id": synthetic.scan_id,
                    "rt": synthetic.rt,
                    "precursor_mz": synthetic.precursor_mz,
                    "precursor_charge": representative_charge,
                    "precursor_intensity": synthetic.precursor_intensity,
                    "precursor_error_ppm": (
                        (component_mass - candidate.neutral_mass)
                        / max(candidate.neutral_mass, 1e-12)
                        * 1_000_000.0
                    ),
                    "activation_method": synthetic.activation_method,
                    "collision_energy": synthetic.collision_energy,
                    **asdict(candidate),
                    "modification_text": candidate.modification_text,
                    **joint,
                    "spectrum_peaks": spectrum_peaks,
                    "search_origin": "component_charge_isotope_consensus_search",
                    "component_group_id": component.get("component_group_id"),
                    "component_neutral_mass": component_mass,
                    "component_confidence": component.get("component_confidence"),
                    "component_source_scan_ids": [
                        str(scan.scan_id) for scan in source_scans
                    ],
                    "component_source_scan_count": len(source_scans),
                }
                candidate_rows.append(row)
            if candidate_rows:
                winners.append(max(candidate_rows, key=lambda row: float(row["score"])))

    searched = assign_q_values(winners)
    accepted: list[dict[str, object]] = []
    for psm in searched:
        q_value = float(psm["q_value"] if psm.get("q_value") is not None else 1.0)
        if psm.get("is_decoy") or q_value > fdr_threshold:
            continue
        if (
            float(psm.get("score") or 0.0) < min_score
            or int(psm.get("matched_ion_count") or 0) < min_matched_ions
            or float(psm.get("fragment_coverage") or 0.0) < min_fragment_coverage
            or int(psm.get("component_scan_support_count") or 0) < 2
            or int(psm.get("component_member_support_count") or 0) < 2
        ):
            continue
        component = next(
            (
                row for row in component_groups
                if str(row.get("component_group_id") or "")
                == str(psm.get("component_group_id") or "")
            ),
            None,
        )
        if component is None:
            continue
        links: list[dict[str, object]] = []
        for member in component.get("members") or []:
            member_id = str(member.get("feature_group_id") or "")
            if not member_id or member_id in excluded:
                continue
            isotope_offset = int(member.get("component_isotope_offset") or 0)
            charge = int(member.get("component_charge") or 0)
            link = _feature_link_payload(
                member,
                "component_charge_isotope_consensus_search",
                isotope_offset,
            )
            link["feature_charge"] = charge or None
            if charge > 0:
                expected_mz = (
                    (float(psm["neutral_mass"]) + isotope_offset * ISOTOPE_MASS_DIFF)
                    / charge
                    + PROTON
                )
                observed_mz = float(
                    member.get("true_peak_mz")
                    or member.get("representative_mz")
                    or 0.0
                )
                link["feature_mz_error_da"] = observed_mz - expected_mz
                link["feature_mz_error_ppm"] = (
                    (observed_mz - expected_mz)
                    / max(expected_mz, 1e-12)
                    * 1_000_000.0
                )
            links.append(link)
        if len(links) < 2:
            continue
        psm["feature_links"] = links
        psm.update({
            key: value
            for key, value in links[0].items()
            if key != "feature_link_distance"
        })
        accepted.append(psm)
    return accepted


def _fragment_charge(label: str) -> int:
    return int(label.rsplit("^", 1)[1]) if "^" in label else 1


def _best_observed_peak(
    observed_mz: list[float],
    observed_intensity: list[float],
    target_mz: float,
    tolerance_ppm: float,
) -> dict[str, object] | None:
    tolerance = max(target_mz * tolerance_ppm / 1_000_000.0, 0.005)
    left = bisect_left(observed_mz, target_mz - tolerance)
    right = bisect_right(observed_mz, target_mz + tolerance)
    if left >= right:
        return None
    index = max(range(left, right), key=observed_intensity.__getitem__)
    return {
        "observed_index": index,
        "observed_mz": observed_mz[index],
        "intensity": observed_intensity[index],
        "error_ppm": (
            (observed_mz[index] - target_mz)
            / max(target_mz, 1e-12)
            * 1_000_000.0
        ),
    }


def _sequence_tag_metrics(
    matches: list[dict[str, object]],
    sequence_length: int,
    explained_intensity: float,
) -> dict[str, object]:
    positions = {
        (str(match["series"]), int(match["ordinal"]))
        for match in matches
    }
    b_positions = {ordinal for series, ordinal in positions if series == "b"}
    y_positions = {ordinal for series, ordinal in positions if series == "y"}
    b_run = _longest_run(b_positions)
    y_run = _longest_run(y_positions)
    longest_run = max(b_run, y_run)
    sequence_tag_length = max(0, longest_run - 1)
    complementary_pairs = sum(
        ("b", ordinal) in positions
        and ("y", sequence_length - ordinal) in positions
        for ordinal in range(1, sequence_length)
    )
    denominator = max(1, 2 * (sequence_length - 1))
    coverage = len(positions) / denominator
    series_count = int(bool(b_positions)) + int(bool(y_positions))
    structural_score = (
        20.0 * min(1.0, coverage / 0.25)
        + 25.0 * min(1.0, sequence_tag_length / 3.0)
        + 15.0 * min(1.0, complementary_pairs / 2.0)
        + 10.0 * min(1.0, series_count / 2.0)
        + 10.0 * min(1.0, max(0.0, explained_intensity))
        + 10.0 * min(1.0, len(positions) / 8.0)
    )
    return {
        "score": structural_score,
        "matched_ion_count": len(positions),
        "fragment_coverage": coverage,
        "explained_intensity": explained_intensity,
        "ion_continuity": longest_run / max(1, sequence_length - 1),
        "sequence_tag_length": sequence_tag_length,
        "b_ion_longest_run": b_run,
        "y_ion_longest_run": y_run,
        "complementary_ion_pair_count": complementary_pairs,
        "fragment_series_count": series_count,
    }


def _match_mass_offset_fragments(
    candidate: PeptideCandidate,
    mass_delta: float,
    observed_mz: list[float],
    observed_intensity: list[float],
    fragment_tolerance_ppm: float,
    fixed_site_index: int | None = None,
) -> dict[str, object]:
    """Match regular and delta-shifted b/y ions and choose one coherent site."""
    options: list[dict[str, object]] = []
    for label, theoretical_mz, series, ordinal in theoretical_fragments(candidate):
        charge = _fragment_charge(label)
        regular = _best_observed_peak(
            observed_mz,
            observed_intensity,
            theoretical_mz,
            fragment_tolerance_ppm,
        )
        shifted_mz = theoretical_mz + mass_delta / charge
        shifted = (
            _best_observed_peak(
                observed_mz,
                observed_intensity,
                shifted_mz,
                fragment_tolerance_ppm,
            )
            if shifted_mz > 0.0 else None
        )
        options.append({
            "label": label,
            "series": series,
            "ordinal": ordinal,
            "charge": charge,
            "regular_mz": theoretical_mz,
            "shifted_mz": shifted_mz,
            "regular": regular,
            "shifted": shifted,
        })

    def site_evidence(site_index: int) -> dict[str, object]:
        possible: list[dict[str, object]] = []
        for option in options:
            series = str(option["series"])
            ordinal = int(option["ordinal"])
            contains_site = (
                abs(mass_delta) > 1e-9
                and (
                    (series == "b" and ordinal > site_index)
                    or (
                        series == "y"
                        and ordinal >= len(candidate.sequence) - site_index
                    )
                )
            )
            peak = option["shifted"] if contains_site else option["regular"]
            if not isinstance(peak, dict):
                continue
            possible.append({
                **peak,
                "label": (
                    f"{option['label']}+Δ"
                    if contains_site else str(option["label"])
                ),
                "series": series,
                "ordinal": ordinal,
                "fragment_charge": option["charge"],
                "theoretical_mz": (
                    option["shifted_mz"]
                    if contains_site else option["regular_mz"]
                ),
                "mass_offset_shifted": contains_site,
            })
        chosen: list[dict[str, object]] = []
        used_indices: set[int] = set()
        used_positions: set[tuple[str, int]] = set()
        for match in sorted(
            possible,
            key=lambda row: float(row["intensity"]),
            reverse=True,
        ):
            index = int(match["observed_index"])
            position = (str(match["series"]), int(match["ordinal"]))
            if index in used_indices or position in used_positions:
                continue
            used_indices.add(index)
            used_positions.add(position)
            chosen.append(match)
        total_intensity = sum(observed_intensity)
        explained = (
            sum(float(match["intensity"]) for match in chosen)
            / total_intensity
            if total_intensity else 0.0
        )
        metrics = _sequence_tag_metrics(
            chosen,
            len(candidate.sequence),
            explained,
        )
        return {
            **metrics,
            "matched_fragments": chosen,
            "mass_offset_site_index": site_index,
            "regular_fragment_count": sum(
                not bool(match["mass_offset_shifted"]) for match in chosen
            ),
            "shifted_fragment_count": sum(
                bool(match["mass_offset_shifted"]) for match in chosen
            ),
        }

    sites = (
        [fixed_site_index]
        if fixed_site_index is not None
        else list(range(len(candidate.sequence)))
    )
    evaluated = [site_evidence(site) for site in sites]
    if not evaluated:
        return {
            **_sequence_tag_metrics([], len(candidate.sequence), 0.0),
            "matched_fragments": [],
            "mass_offset_site_index": None,
            "localized_mass_offset_site": None,
            "mass_offset_localization_margin": 0.0,
            "regular_fragment_count": 0,
            "shifted_fragment_count": 0,
        }
    evaluated.sort(key=lambda row: float(row["score"]), reverse=True)
    best = evaluated[0]
    margin = (
        float(best["score"]) - float(evaluated[1]["score"])
        if len(evaluated) > 1 else float(best["score"])
    )
    localized = (
        int(best["mass_offset_site_index"]) + 1
        if margin >= 2.0
        and int(best["regular_fragment_count"]) > 0
        and int(best["shifted_fragment_count"]) > 0
        else None
    )
    return {
        **best,
        "localized_mass_offset_site": localized,
        "mass_offset_localization_margin": margin,
    }


def _screen_mass_offset_fragments(
    candidate: PeptideCandidate,
    mass_delta: float,
    observed_mz: list[float],
    observed_intensity: list[float],
    fragment_tolerance_ppm: float,
) -> dict[str, object]:
    """Fast optimistic screen before coherent site localization."""
    possible: list[dict[str, object]] = []
    for label, theoretical_mz, series, ordinal in theoretical_fragments(
        candidate
    ):
        charge = _fragment_charge(label)
        choices: list[tuple[bool, float, dict[str, object]]] = []
        regular = _best_observed_peak(
            observed_mz,
            observed_intensity,
            theoretical_mz,
            fragment_tolerance_ppm,
        )
        if regular is not None:
            choices.append((False, theoretical_mz, regular))
        shifted_mz = theoretical_mz + mass_delta / charge
        if shifted_mz > 0.0 and abs(mass_delta) > 1e-9:
            shifted = _best_observed_peak(
                observed_mz,
                observed_intensity,
                shifted_mz,
                fragment_tolerance_ppm,
            )
            if shifted is not None:
                choices.append((True, shifted_mz, shifted))
        if not choices:
            continue
        is_shifted, target_mz, peak = max(
            choices,
            key=lambda item: float(item[2]["intensity"]),
        )
        possible.append({
            **peak,
            "label": f"{label}+Δ" if is_shifted else label,
            "series": series,
            "ordinal": ordinal,
            "fragment_charge": charge,
            "theoretical_mz": target_mz,
            "mass_offset_shifted": is_shifted,
        })
    chosen: list[dict[str, object]] = []
    used_indices: set[int] = set()
    used_positions: set[tuple[str, int]] = set()
    for match in sorted(
        possible,
        key=lambda row: float(row["intensity"]),
        reverse=True,
    ):
        index = int(match["observed_index"])
        position = (str(match["series"]), int(match["ordinal"]))
        if index in used_indices or position in used_positions:
            continue
        used_indices.add(index)
        used_positions.add(position)
        chosen.append(match)
    total_intensity = sum(observed_intensity)
    explained = (
        sum(float(match["intensity"]) for match in chosen)
        / total_intensity
        if total_intensity else 0.0
    )
    return {
        **_sequence_tag_metrics(
            chosen,
            len(candidate.sequence),
            explained,
        ),
        "matched_fragments": chosen,
    }


def _mass_offset_description(
    candidate: PeptideCandidate,
    mass_delta: float,
    localized_site: int | None,
) -> tuple[str, list[str]]:
    truncation_labels = {
        "n_terminal_truncation_candidate": "N-terminal truncation candidate",
        "c_terminal_truncation_candidate": "C-terminal truncation candidate",
    }
    truncation_text = truncation_labels.get(candidate.proteolysis)
    if truncation_text and abs(mass_delta) < 0.5:
        return (f"{truncation_text}; no additional ΔMass", [])
    compatible: list[str] = []
    for name, residues, delta, location in COMMON_MODIFICATIONS:
        if abs(mass_delta - delta) > 0.05:
            continue
        positions = (
            [localized_site - 1]
            if localized_site is not None else list(range(len(candidate.sequence)))
        )
        for position in positions:
            aa = candidate.sequence[position]
            valid = (
                (location == "residue" and aa in residues)
                or (location == "n_term" and position == 0 and aa in residues)
                or (
                    location == "n_glycan_sequon"
                    and aa in residues
                    and position + 2 < len(candidate.sequence)
                    and candidate.sequence[position + 1] != "P"
                    and candidate.sequence[position + 2] in "ST"
                )
            )
            if valid:
                compatible.append(name)
                break
    compatible = sorted(set(compatible))
    site_text = (
        f"near {candidate.sequence[localized_site - 1]}{localized_site}"
        if localized_site is not None else "site unresolved"
    )
    compatibility = (
        f"; compatible with {' / '.join(compatible)}"
        if compatible else ""
    )
    prefix = "Open" if compatible else "Unknown"
    truncation_prefix = f"{truncation_text}; " if truncation_text else ""
    return (
        f"{truncation_prefix}{prefix} ΔMass {mass_delta:+.4f} Da"
        f"{compatibility}; {site_text}",
        compatible,
    )


def _unique_compatible_mass_offset_site(
    candidate: PeptideCandidate,
    mass_delta: float,
) -> int | None:
    """Return a unique chemically compatible 1-based site, if one exists."""
    positions: set[int] = set()
    for _, residues, delta, location in COMMON_MODIFICATIONS:
        if abs(mass_delta - delta) > 0.05:
            continue
        for position, aa in enumerate(candidate.sequence):
            valid = (
                (location == "residue" and aa in residues)
                or (location == "n_term" and position == 0 and aa in residues)
                or (
                    location == "n_glycan_sequon"
                    and aa in residues
                    and position + 2 < len(candidate.sequence)
                    and candidate.sequence[position + 1] != "P"
                    and candidate.sequence[position + 2] in "ST"
                )
            )
            if valid:
                positions.add(position + 1)
    return next(iter(positions)) if len(positions) == 1 else None


def _mass_offset_display_label(
    candidate: PeptideCandidate,
    mass_delta: float,
    localized_site: int | None,
    compatible_offsets: list[str],
) -> str:
    """Return a concise UI label without hiding an unresolved open mass."""
    if compatible_offsets:
        names = " / ".join(sorted(set(str(value) for value in compatible_offsets)))
        if localized_site is not None and 1 <= localized_site <= len(candidate.sequence):
            return f"{names}@{candidate.sequence[localized_site - 1]}{localized_site}"
        return f"{names} (site unresolved)"
    return f"Open ΔMass {mass_delta:+.4f} Da"


def _joint_mass_offset_evidence(
    source_records: list[dict[str, object]],
    candidate: PeptideCandidate,
    mass_delta: float,
    site_index: int,
    fragment_tolerance_ppm: float,
) -> dict[str, object]:
    label_matches: dict[tuple[str, int], dict[str, object]] = {}
    informative_records: list[dict[str, object]] = []
    explained_values: list[float] = []
    for record in source_records:
        scan = record["scan"]
        assert isinstance(scan, LCMSSpectrumScan)
        pairs = sorted(
            zip(scan.mz_array, scan.intensity_array),
            key=lambda pair: pair[0],
        )
        evidence = _match_mass_offset_fragments(
            candidate,
            mass_delta,
            [float(pair[0]) for pair in pairs],
            [max(0.0, float(pair[1])) for pair in pairs],
            fragment_tolerance_ppm,
            fixed_site_index=site_index,
        )
        if int(evidence["matched_ion_count"]) < 2:
            continue
        informative_records.append(record)
        explained_values.append(float(evidence["explained_intensity"]))
        for match in evidence["matched_fragments"]:
            position = (str(match["series"]), int(match["ordinal"]))
            previous = label_matches.get(position)
            if (
                previous is None
                or float(match["intensity"]) > float(previous["intensity"])
            ):
                label_matches[position] = dict(match)
    matches = list(label_matches.values())
    explained = (
        statistics.mean(sorted(explained_values, reverse=True)[:4])
        if explained_values else 0.0
    )
    metrics = _sequence_tag_metrics(
        matches,
        len(candidate.sequence),
        explained,
    )
    member_ids = {
        str(record["member_id"]) for record in informative_records
    }
    charge_states = {
        int(record["member_charge"])
        for record in informative_records
        if int(record.get("member_charge") or 0) > 0
    }
    support_score = (
        5.0 * min(1.0, len(informative_records) / 2.0)
        + 5.0 * min(1.0, len(member_ids) / 2.0)
    )
    return {
        **metrics,
        "score": float(metrics["score"]) + support_score,
        "matched_fragments": matches,
        "component_scan_support_count": len(informative_records),
        "component_member_support_count": len(member_ids),
        "component_charge_support_count": len(charge_states),
        "component_source_member_ids": sorted(member_ids),
        "component_source_charge_states": sorted(charge_states),
    }


def search_component_sequence_tag_scans(
    payload: dict[str, object],
    scans: list[LCMSSpectrumScan],
    backbone_candidates: list[PeptideCandidate],
    component_groups: list[dict[str, object]],
    exclude_feature_ids: set[str] | None = None,
    fragment_tolerance_ppm: float = 20.0,
    rt_tolerance_min: float = 0.5,
    precursor_relation_ppm: float = 20.0,
    minimum_mass_delta_da: float = -250.0,
    maximum_mass_delta_da: float = 2500.0,
    minimum_absolute_delta_da: float = 0.5,
    max_scans_per_member_sample: int = 4,
    max_screened_candidates: int = 12,
    fdr_threshold: float = 0.01,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Infer a known-sequence backbone using an unknown mass offset and b/y tags.

    Only selected precursor/isotope scans assigned to at least two component
    members are used. Accepted backbone candidates and weaker sequence-region
    candidates are returned separately.
    """
    excluded = {str(value) for value in (exclude_feature_ids or set())}
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
    for scan in scans:
        if scan.ms_level == 2 and scan.precursor_mz is not None and scan.mz_array:
            scans_by_sample.setdefault(str(scan.sample_id), []).append(scan)
    candidates = [
        candidate for candidate in backbone_candidates
        if len(candidate.sequence) >= 6
    ]
    winners: list[dict[str, object]] = []

    for component in component_groups:
        if not component.get("inferred_component"):
            continue
        component_mass = float(component.get("component_neutral_mass") or 0.0)
        if component_mass <= 0.0:
            continue
        members = [
            member for member in component.get("members") or []
            if str(member.get("feature_group_id") or "") not in excluded
        ]
        if len(members) < 2:
            continue
        local_candidates = [
            (candidate, component_mass - candidate.neutral_mass)
            for candidate in candidates
            if minimum_mass_delta_da
            <= component_mass - candidate.neutral_mass
            <= maximum_mass_delta_da
            and (
                abs(component_mass - candidate.neutral_mass)
                >= minimum_absolute_delta_da
                or candidate.proteolysis in {
                    "n_terminal_truncation_candidate",
                    "c_terminal_truncation_candidate",
                }
            )
        ]
        if not local_candidates:
            continue

        for sample_id, sample_scans in scans_by_sample.items():
            best_assignment_by_scan: dict[
                str,
                tuple[tuple[float, ...], dict[str, object]],
            ] = {}
            for member in members:
                member_id = str(member.get("feature_group_id") or "")
                member_mz = float(
                    member.get("true_peak_mz")
                    or member.get("representative_mz")
                    or 0.0
                )
                member_rt = float(member.get("representative_rt") or 0.0)
                member_charge = int(member.get("component_charge") or 0)
                local_shift = float(
                    dict(member.get("rt_correction_by_sample") or {}).get(
                        sample_id
                    )
                    or 0.0
                )
                for scan in sample_scans:
                    if (
                        member_charge > 0
                        and scan.precursor_charge
                        and int(scan.precursor_charge) != member_charge
                    ):
                        continue
                    aligned_rt = (
                        float(scan.rt)
                        + float(shifts.get(sample_id) or 0.0)
                        + local_shift
                    )
                    rt_error = aligned_rt - member_rt
                    if abs(rt_error) > max(0.0, rt_tolerance_min):
                        continue
                    mz_tolerance = max(
                        member_mz
                        * max(precursor_relation_ppm, 0.1)
                        / 1_000_000.0,
                        0.02,
                    )
                    relation = _scan_isotope_relation(
                        member_mz,
                        float(scan.precursor_mz or 0.0),
                        scan.precursor_charge,
                        mz_tolerance,
                    )
                    if relation is None:
                        continue
                    isotope_offset, mz_error = relation
                    key = (
                        float(abs(isotope_offset)),
                        abs(float(mz_error)) / max(mz_tolerance, 1e-12),
                        abs(rt_error) / max(rt_tolerance_min, 1e-12),
                        -float(scan.precursor_intensity or 0.0),
                    )
                    record = {
                        "scan": scan,
                        "member_id": member_id,
                        "member_charge": (
                            member_charge or int(scan.precursor_charge or 0)
                        ),
                    }
                    previous = best_assignment_by_scan.get(str(scan.scan_id))
                    if previous is None or key < previous[0]:
                        best_assignment_by_scan[str(scan.scan_id)] = (
                            key,
                            record,
                        )
            records_by_member: dict[
                str,
                list[tuple[tuple[float, ...], dict[str, object]]],
            ] = {}
            for key, record in best_assignment_by_scan.values():
                records_by_member.setdefault(
                    str(record["member_id"]),
                    [],
                ).append((key, record))
            source_records = [
                record
                for member_records in records_by_member.values()
                for _, record in sorted(
                    member_records,
                    key=lambda item: item[0],
                )[:max(1, int(max_scans_per_member_sample))]
            ]
            if (
                len(source_records) < 2
                or len({
                    str(record["member_id"]) for record in source_records
                }) < 2
            ):
                continue
            source_scans = [
                record["scan"] for record in source_records
                if isinstance(record.get("scan"), LCMSSpectrumScan)
            ]
            consensus_mz, consensus_intensity = _merge_consensus_ms2(
                source_scans,
                fragment_tolerance_ppm,
                minimum_scan_support=1,
            )
            if not consensus_mz:
                continue

            optimistic: list[
                tuple[float, PeptideCandidate, float, dict[str, object]]
            ] = []
            for candidate, mass_delta in local_candidates:
                evidence = _screen_mass_offset_fragments(
                    candidate,
                    mass_delta,
                    consensus_mz,
                    consensus_intensity,
                    fragment_tolerance_ppm,
                )
                if (
                    int(evidence["matched_ion_count"]) < 3
                    or (
                        int(evidence["sequence_tag_length"]) < 1
                        and int(evidence["complementary_ion_pair_count"]) < 1
                    )
                ):
                    continue
                optimistic.append((
                    float(evidence["score"]),
                    candidate,
                    mass_delta,
                    evidence,
                ))
            screened: list[
                tuple[float, PeptideCandidate, float, dict[str, object]]
            ] = []
            for _, candidate, mass_delta, __ in sorted(
                optimistic,
                key=lambda item: item[0],
                reverse=True,
            )[:max(4, int(max_screened_candidates) * 4)]:
                evidence = _match_mass_offset_fragments(
                    candidate,
                    mass_delta,
                    consensus_mz,
                    consensus_intensity,
                    fragment_tolerance_ppm,
                )
                if (
                    int(evidence["matched_ion_count"]) < 3
                    or (
                        int(evidence["sequence_tag_length"]) < 1
                        and int(evidence["complementary_ion_pair_count"]) < 1
                    )
                ):
                    continue
                screened.append((
                    float(evidence["score"]),
                    candidate,
                    mass_delta,
                    evidence,
                ))
            evaluated: list[dict[str, object]] = []
            for _, candidate, mass_delta, screen in sorted(
                screened,
                key=lambda item: item[0],
                reverse=True,
            )[:max(1, int(max_screened_candidates))]:
                site_index = int(screen["mass_offset_site_index"])
                joint = _joint_mass_offset_evidence(
                    source_records,
                    candidate,
                    mass_delta,
                    site_index,
                    fragment_tolerance_ppm,
                )
                if int(joint["matched_ion_count"]) < 3:
                    continue
                localized_site = (
                    int(screen["localized_mass_offset_site"])
                    if screen.get("localized_mass_offset_site") is not None
                    else None
                )
                modification_text, compatible_offsets = (
                    _mass_offset_description(
                        candidate,
                        mass_delta,
                        localized_site,
                    )
                )
                display_modification_text = _mass_offset_display_label(
                    candidate,
                    mass_delta,
                    localized_site,
                    compatible_offsets,
                )
                candidate_model_prior = (
                    5.0
                    if (
                        candidate.proteolysis in {
                            "n_terminal_truncation_candidate",
                            "c_terminal_truncation_candidate",
                        }
                        and abs(mass_delta) < 0.5
                    )
                    else (2.0 if compatible_offsets else 0.0)
                )
                joint["score"] = (
                    float(joint["score"]) + candidate_model_prior
                )
                display = _match_mass_offset_fragments(
                    candidate,
                    mass_delta,
                    consensus_mz,
                    consensus_intensity,
                    fragment_tolerance_ppm,
                    fixed_site_index=site_index,
                )
                labels = {
                    int(match["observed_index"]): str(match["label"])
                    for match in display["matched_fragments"]
                }
                top_indices = sorted(
                    range(len(consensus_mz)),
                    key=consensus_intensity.__getitem__,
                    reverse=True,
                )[:60]
                spectrum_peaks = [
                    {
                        "mz": consensus_mz[index],
                        "intensity": consensus_intensity[index],
                        "label": labels.get(index, ""),
                    }
                    for index in sorted(
                        top_indices,
                        key=consensus_mz.__getitem__,
                    )
                ]
                representative_charge = int(
                    component.get("component_representative_charge")
                    or next(
                        (
                            record["member_charge"]
                            for record in source_records
                            if int(record.get("member_charge") or 0) > 0
                        ),
                        2,
                    )
                )
                row = {
                    **asdict(candidate),
                    "candidate_id": (
                        f"{candidate.candidate_id}:OPEN_DELTA="
                        f"{mass_delta:+.4f}"
                    ),
                    "neutral_mass": component_mass,
                    "backbone_neutral_mass": candidate.neutral_mass,
                    "modification_text": modification_text,
                    "mass_offset_description": modification_text,
                    "mass_delta": mass_delta,
                    "compatible_mass_offsets": compatible_offsets,
                    "candidate_model_prior": candidate_model_prior,
                    "localized_mass_offset_site": localized_site,
                    "mass_offset_site_index": site_index,
                    "mass_offset_localization_margin": screen[
                        "mass_offset_localization_margin"
                    ],
                    "sample_id": sample_id,
                    "scan_id": (
                        f"component_tag={component.get('component_group_id')}"
                        f"|sample={sample_id}"
                    ),
                    "rt": statistics.median(
                        float(scan.rt) for scan in source_scans
                    ),
                    "precursor_mz": (
                        component_mass / representative_charge + PROTON
                    ),
                    "precursor_charge": representative_charge,
                    "precursor_intensity": max(
                        (
                            float(scan.precursor_intensity or 0.0)
                            for scan in source_scans
                        ),
                        default=0.0,
                    ),
                    "precursor_error_ppm": 0.0,
                    "activation_method": source_scans[0].activation_method,
                    "collision_energy": source_scans[0].collision_energy,
                    **joint,
                    "spectrum_peaks": spectrum_peaks,
                    "search_origin": (
                        "component_mass_offset_sequence_tag_search"
                    ),
                    "component_group_id": component.get("component_group_id"),
                    "component_neutral_mass": component_mass,
                    "component_confidence": component.get(
                        "component_confidence"
                    ),
                    "component_source_scan_ids": [
                        str(scan.scan_id) for scan in source_scans
                    ],
                    "component_source_scan_count": len(source_scans),
                }
                evaluated.append(row)
            if not evaluated:
                continue
            evaluated.sort(
                key=lambda row: float(row["score"]),
                reverse=True,
            )
            winner = evaluated[0]
            competing_scores = [
                float(row["score"]) for row in evaluated[1:]
                if (
                    row.get("sequence") != winner.get("sequence")
                    or row.get("chain") != winner.get("chain")
                    or row.get("start") != winner.get("start")
                )
            ]
            winner["candidate_score_margin"] = (
                float(winner["score"]) - max(competing_scores)
                if competing_scores else float(winner["score"])
            )
            winner["alternative_sequence_candidates"] = [
                {
                    "sequence": row.get("sequence"),
                    "chain": row.get("chain"),
                    "start": row.get("start"),
                    "end": row.get("end"),
                    "mass_delta": row.get("mass_delta"),
                    "score": row.get("score"),
                    "is_decoy": row.get("is_decoy"),
                }
                for row in evaluated[1:4]
            ]
            links: list[dict[str, object]] = []
            for member in members:
                member_id = str(member.get("feature_group_id") or "")
                if not member_id:
                    continue
                isotope_offset = int(
                    member.get("component_isotope_offset") or 0
                )
                charge = int(member.get("component_charge") or 0)
                link = _feature_link_payload(
                    member,
                    "component_mass_offset_sequence_tag_search",
                    isotope_offset,
                )
                link["feature_charge"] = charge or None
                if charge > 0:
                    expected_mz = (
                        (
                            component_mass
                            + isotope_offset * ISOTOPE_MASS_DIFF
                        )
                        / charge
                        + PROTON
                    )
                    observed_mz = float(
                        member.get("true_peak_mz")
                        or member.get("representative_mz")
                        or 0.0
                    )
                    link["feature_mz_error_da"] = (
                        observed_mz - expected_mz
                    )
                    link["feature_mz_error_ppm"] = (
                        (observed_mz - expected_mz)
                        / max(expected_mz, 1e-12)
                        * 1_000_000.0
                    )
                links.append(link)
            winner["feature_links"] = links
            if links:
                winner.update({
                    key: value for key, value in links[0].items()
                    if key != "feature_link_distance"
                })
            winners.append(winner)

    searched = assign_q_values(winners)
    accepted: list[dict[str, object]] = []
    region_candidates: list[dict[str, object]] = []
    for row in searched:
        if row.get("is_decoy"):
            continue
        q_value = float(
            row["q_value"] if row.get("q_value") is not None else 1.0
        )
        score = float(row.get("score") or 0.0)
        ions = int(row.get("matched_ion_count") or 0)
        tag_length = int(row.get("sequence_tag_length") or 0)
        complementary = int(
            row.get("complementary_ion_pair_count") or 0
        )
        b_run = int(row.get("b_ion_longest_run") or 0)
        y_run = int(row.get("y_ion_longest_run") or 0)
        scan_support = int(row.get("component_scan_support_count") or 0)
        member_support = int(
            row.get("component_member_support_count") or 0
        )
        margin = float(row.get("candidate_score_margin") or 0.0)
        sequence_length = len(str(row.get("sequence") or ""))
        mass_delta = abs(float(row.get("mass_delta") or 0.0))
        exact_terminal_truncation = (
            row.get("proteolysis") in {
                "n_terminal_truncation_candidate",
                "c_terminal_truncation_candidate",
            }
            and mass_delta < 0.5
        )
        specific_enough_mass_model = (
            sequence_length >= 8
            or exact_terminal_truncation
            or bool(row.get("compatible_mass_offsets"))
        )
        backbone_supported = (
            q_value <= fdr_threshold
            and score >= 45.0
            and ions >= 5
            and float(row.get("fragment_coverage") or 0.0) >= 0.06
            and tag_length >= 2
            and (
                tag_length >= 3
                or complementary >= 1
                or (b_run >= 3 and y_run >= 3)
            )
            and scan_support >= 2
            and member_support >= 2
            and margin >= 3.0
            and specific_enough_mass_model
        )
        if backbone_supported:
            row["sequence_inference_level"] = (
                "truncation_sequence_supported"
                if row.get("proteolysis") in {
                    "n_terminal_truncation_candidate",
                    "c_terminal_truncation_candidate",
                }
                and abs(float(row.get("mass_delta") or 0.0)) < 0.5
                else "backbone_sequence_supported"
            )
            accepted.append(row)
            continue
        sequence_region_supported = (
            q_value <= max(0.05, fdr_threshold)
            and score >= 32.0
            and ions >= 4
            and (tag_length >= 2 or complementary >= 1)
            and scan_support >= 2
            and member_support >= 2
            and margin >= 1.0
        )
        if sequence_region_supported:
            row["sequence_inference_level"] = (
                "sequence_region_candidate"
            )
            region_candidates.append(row)
    return accepted, region_candidates


def search_feature_guided_scans(
    payload: dict[str, object],
    scans: list[LCMSSpectrumScan],
    candidates: list[PeptideCandidate],
    exclude_feature_ids: set[str] | None = None,
    precursor_tolerance_ppm: float = 20.0,
    fragment_tolerance_ppm: float = 20.0,
    rt_tolerance_min: float = 0.5,
    isolation_padding_da: float = 0.02,
    max_scans_per_sample: int = 2,
    fdr_threshold: float = 0.01,
) -> list[dict[str, object]]:
    """Re-search unresolved MS1 features in co-isolating MS2 scans.

    These matches are deliberately marked tentative because the selected ion can
    differ from the MS1 feature inside the same isolation window.
    """
    excluded = exclude_feature_ids or set()
    significant_types = {"presence_absence", "area_changed", "moderate_difference"}
    features = [
        row for row in payload.get("global_feature_groups") or []
        if row.get("difference_type") in significant_types and str(row.get("feature_group_id")) not in excluded
    ]
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
    for scan in scans:
        if scan.ms_level == 2 and scan.precursor_mz is not None and scan.mz_array:
            scans_by_sample.setdefault(scan.sample_id, []).append(scan)
    synthetic_scans: list[LCMSSpectrumScan] = []
    metadata: dict[str, tuple[dict[str, object], LCMSSpectrumScan]] = {}
    counter = 0
    for feature in features:
        feature_mz = float(feature.get("representative_mz") or 0.0)
        feature_rt = float(feature.get("representative_rt") or 0.0)
        local_shifts = dict(feature.get("rt_correction_by_sample") or {})
        for sample_id, sample_scans in scans_by_sample.items():
            possible: list[tuple[float, float, LCMSSpectrumScan]] = []
            center_rt = feature_rt - float(shifts.get(sample_id) or 0.0) - float(local_shifts.get(sample_id) or 0.0)
            candidate_scans = [
                scan for scan in sample_scans
                if abs(float(scan.rt) - center_rt) <= max(0.0, rt_tolerance_min)
            ]
            for scan in candidate_scans:
                aligned_rt = scan.rt + float(shifts.get(sample_id) or 0.0) + float(local_shifts.get(sample_id) or 0.0)
                rt_error = aligned_rt - feature_rt
                if abs(rt_error) > rt_tolerance_min:
                    continue
                precursor_mz = float(scan.precursor_mz or 0.0)
                lower = float(scan.isolation_window_lower_offset or 0.8)
                upper = float(scan.isolation_window_upper_offset or 0.8)
                if not precursor_mz - lower - isolation_padding_da <= feature_mz <= precursor_mz + upper + isolation_padding_da:
                    continue
                tolerance = max(feature_mz * precursor_tolerance_ppm / 1_000_000.0, isolation_padding_da)
                if _scan_isotope_relation(feature_mz, precursor_mz, scan.precursor_charge, tolerance) is not None:
                    continue
                window = upper if feature_mz >= precursor_mz else lower
                distance = abs(rt_error) / max(rt_tolerance_min, 1e-9) + abs(feature_mz - precursor_mz) / max(window, 1e-9)
                possible.append((distance, -float(scan.precursor_intensity or 0.0), scan))
            for _, __, scan in sorted(possible)[:max(0, max_scans_per_sample)]:
                synthetic_id = f"feature_guided_scan={counter}"
                counter += 1
                synthetic_scans.append(replace(
                    scan,
                    scan_id=synthetic_id,
                    precursor_mz=feature_mz,
                    precursor_charge=None,
                ))
                metadata[synthetic_id] = (feature, scan)
    if not synthetic_scans:
        return []
    searched = assign_q_values(search_scans(
        synthetic_scans,
        candidates,
        precursor_tolerance_ppm=precursor_tolerance_ppm,
        fragment_tolerance_ppm=fragment_tolerance_ppm,
        min_score=15.0,
    ))
    accepted: list[dict[str, object]] = []
    for psm in searched:
        q_value = float(psm["q_value"] if psm.get("q_value") is not None else 1.0)
        if psm.get("is_decoy") or q_value > fdr_threshold:
            continue
        if float(psm.get("score") or 0.0) < 40.0 or int(psm.get("matched_ion_count") or 0) < 4:
            continue
        feature, original_scan = metadata[str(psm["scan_id"])]
        link = _feature_link_payload(feature, "isolation_window_targeted_search")
        psm.update({
            "targeted_search_scan_id": psm["scan_id"],
            "scan_id": original_scan.scan_id,
            "rt": original_scan.rt,
            "acquisition_precursor_mz": original_scan.precursor_mz,
            "acquisition_precursor_charge": original_scan.precursor_charge,
            "isolation_window_lower_offset": original_scan.isolation_window_lower_offset,
            "isolation_window_upper_offset": original_scan.isolation_window_upper_offset,
            "search_origin": "feature_guided_isolation_window",
            "feature_links": [link],
            **{key: value for key, value in link.items() if key != "feature_link_distance"},
        })
        accepted.append(psm)
    return accepted


def _n_glycan_sequon_sites(sequence: str) -> list[int]:
    return [
        position
        for position, aa in enumerate(sequence)
        if aa == "N"
        and position + 2 < len(sequence)
        and sequence[position + 1] != "P"
        and sequence[position + 2] in "ST"
    ]


def _targeted_n_glycan_candidates(
    backbone_candidates: list[PeptideCandidate],
) -> list[PeptideCandidate]:
    candidates: list[PeptideCandidate] = []
    seen: set[tuple[str, int, str]] = set()
    for backbone in backbone_candidates:
        if backbone.is_decoy or backbone.modifications:
            continue
        for site_index in _n_glycan_sequon_sites(backbone.sequence):
            for glycan in TARGETED_N_GLYCANS:
                name = str(glycan["name"])
                key = (backbone.base_peptide_id, site_index, name)
                if key in seen:
                    continue
                seen.add(key)
                target = _candidate(
                    backbone.chain,
                    backbone.start,
                    backbone.end,
                    backbone.sequence,
                    ((site_index, name, float(glycan["mass"])),),
                    carbamidomethyl_cys=backbone.carbamidomethyl_cys,
                    proteolysis=backbone.proteolysis,
                    base_peptide_id=backbone.base_peptide_id,
                )
                candidates.extend((target, _decoy(target)))
    return sorted(candidates, key=lambda candidate: candidate.neutral_mass)


def _targeted_glycan_diagnostic_evidence(
    observed_mz: list[float],
    observed_intensity: list[float],
    fragment_tolerance_ppm: float,
) -> dict[str, object]:
    diagnostic_ions = (
        ("HexNAc-126", 126.0550),
        ("HexNAc-138", 138.0550),
        ("HexNAc-144", 144.0655),
        ("HexNAc-168", 168.0655),
        ("HexNAc-186", 186.0761),
        ("HexNAc", 204.086649),
        ("NeuAc-H2O", 274.092127),
        ("NeuAc", 292.102692),
        ("HexHexNAc", 366.139472),
    )
    matches: list[dict[str, object]] = []
    used: set[int] = set()
    for label, target_mz in diagnostic_ions:
        match = _best_observed_peak(
            observed_mz,
            observed_intensity,
            target_mz,
            fragment_tolerance_ppm,
        )
        if match is None or int(match["observed_index"]) in used:
            continue
        used.add(int(match["observed_index"]))
        matches.append({
            **match,
            "label": label,
            "theoretical_mz": target_mz,
            "series": "glycan_diagnostic",
            "ordinal": len(matches) + 1,
        })
    total_intensity = sum(observed_intensity)
    diagnostic_intensity_fraction = (
        sum(float(match["intensity"]) for match in matches)
        / total_intensity
        if total_intensity else 0.0
    )
    labels = {str(match["label"]) for match in matches}
    return {
        "glycan_diagnostic_ion_count": len(matches),
        "glycan_diagnostic_intensity_fraction": diagnostic_intensity_fraction,
        "glycan_has_hexnac_204": "HexNAc" in labels,
        "glycan_has_sialic_diagnostic": bool(
            labels.intersection({"NeuAc-H2O", "NeuAc"})
        ),
        "glycan_diagnostic_matches": matches,
    }


def _targeted_glycan_core_y_evidence(
    backbone: PeptideCandidate,
    composition: dict[str, int],
    observed_mz: list[float],
    observed_intensity: list[float],
    precursor_charge: int,
    fragment_tolerance_ppm: float,
) -> dict[str, object]:
    core_stubs = [
        ("Y0", {}),
        ("Y1", {"HexNAc": 1}),
        ("Y2", {"HexNAc": 2}),
        ("Y3", {"HexNAc": 2, "Hex": 1}),
        ("Y4", {"HexNAc": 2, "Hex": 2}),
        ("Y5", {"HexNAc": 2, "Hex": 3}),
    ]
    if int(composition.get("Fuc") or 0) > 0:
        core_stubs.extend((
            ("Y1F", {"HexNAc": 1, "Fuc": 1}),
            ("Y2F", {"HexNAc": 2, "Fuc": 1}),
            ("Y3F", {"HexNAc": 2, "Hex": 1, "Fuc": 1}),
        ))
    matches: list[dict[str, object]] = []
    used: set[int] = set()
    for label, stub in core_stubs:
        if any(
            int(stub_count) > int(composition.get(name) or 0)
            for name, stub_count in stub.items()
        ):
            continue
        stub_mass = _glycan_composition_mass(stub)
        for charge in range(1, max(1, min(3, precursor_charge)) + 1):
            target_mz = (
                backbone.neutral_mass + stub_mass + charge * PROTON
            ) / charge
            match = _best_observed_peak(
                observed_mz,
                observed_intensity,
                target_mz,
                fragment_tolerance_ppm,
            )
            if match is None or int(match["observed_index"]) in used:
                continue
            used.add(int(match["observed_index"]))
            matches.append({
                **match,
                "label": f"{label}^{charge}" if charge > 1 else label,
                "theoretical_mz": target_mz,
                "series": "glycan_core_y",
                "ordinal": len(matches) + 1,
                "fragment_charge": charge,
            })
    labels = {str(match["label"]).split("^", 1)[0] for match in matches}
    return {
        "glycan_core_y_ion_count": len(matches),
        "glycan_core_y_type_count": len(labels),
        "glycan_has_y1_support": bool(labels.intersection({"Y1", "Y1F"})),
        "glycan_core_y_matches": matches,
    }


def _targeted_glycan_backbone_evidence(
    backbone: PeptideCandidate,
    glycan_site_index: int,
    observed_mz: list[float],
    observed_intensity: list[float],
    fragment_tolerance_ppm: float,
) -> dict[str, object]:
    possible: list[dict[str, object]] = []
    for label, theoretical_mz, series, ordinal in theoretical_fragments(
        backbone,
        max_charge=2,
    ):
        charge = _fragment_charge(label)
        regular = _best_observed_peak(
            observed_mz,
            observed_intensity,
            theoretical_mz,
            fragment_tolerance_ppm,
        )
        if regular is not None:
            possible.append({
                **regular,
                "label": label,
                "series": series,
                "ordinal": ordinal,
                "fragment_charge": charge,
                "theoretical_mz": theoretical_mz,
                "hexnac_retained": False,
            })
        contains_site = (
            (series == "b" and ordinal > glycan_site_index)
            or (
                series == "y"
                and ordinal >= len(backbone.sequence) - glycan_site_index
            )
        )
        if not contains_site:
            continue
        shifted_mz = theoretical_mz + GLYCAN_RESIDUE_MASS["HexNAc"] / charge
        shifted = _best_observed_peak(
            observed_mz,
            observed_intensity,
            shifted_mz,
            fragment_tolerance_ppm,
        )
        if shifted is not None:
            possible.append({
                **shifted,
                "label": f"{label}+HexNAc",
                "series": series,
                "ordinal": ordinal,
                "fragment_charge": charge,
                "theoretical_mz": shifted_mz,
                "hexnac_retained": True,
            })
    chosen: list[dict[str, object]] = []
    used_indices: set[int] = set()
    used_positions: set[tuple[str, int]] = set()
    for match in sorted(
        possible,
        key=lambda row: (
            bool(row["hexnac_retained"]),
            float(row["intensity"]),
        ),
        reverse=True,
    ):
        observed_index = int(match["observed_index"])
        position = (str(match["series"]), int(match["ordinal"]))
        if observed_index in used_indices or position in used_positions:
            continue
        used_indices.add(observed_index)
        used_positions.add(position)
        chosen.append(match)
    total_intensity = sum(observed_intensity)
    explained_intensity = (
        sum(float(match["intensity"]) for match in chosen) / total_intensity
        if total_intensity else 0.0
    )
    metrics = _sequence_tag_metrics(
        chosen,
        len(backbone.sequence),
        explained_intensity,
    )
    return {
        **metrics,
        "matched_fragments": chosen,
        "glycan_hexnac_fragment_count": sum(
            bool(match["hexnac_retained"]) for match in chosen
        ),
    }


def search_feature_glycopeptide_scans(
    payload: dict[str, object],
    scans: list[LCMSSpectrumScan],
    backbone_candidates: list[PeptideCandidate],
    exclude_feature_ids: set[str] | None = None,
    rt_tolerance_min: float = 0.5,
    precursor_tolerance_ppm: float = 20.0,
    fragment_tolerance_ppm: float = 20.0,
    max_scans_per_sample: int = 2,
    fdr_threshold: float = 0.01,
) -> list[dict[str, object]]:
    """Identify N-glycopeptides only for differential Features with linked MS2."""
    excluded = {str(value) for value in (exclude_feature_ids or set())}
    significant_types = {"presence_absence", "area_changed", "moderate_difference"}
    features = [
        row for row in payload.get("global_feature_groups") or []
        if row.get("difference_type") in significant_types
        and str(row.get("feature_group_id") or "") not in excluded
    ]
    glycopeptides = _targeted_n_glycan_candidates(backbone_candidates)
    if not features or not glycopeptides:
        return []
    glycopeptide_masses = [candidate.neutral_mass for candidate in glycopeptides]
    glycan_by_name = {
        str(glycan["name"]): glycan for glycan in TARGETED_N_GLYCANS
    }
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
    for scan in scans:
        if scan.ms_level == 2 and scan.precursor_mz is not None and scan.mz_array:
            scans_by_sample.setdefault(str(scan.sample_id), []).append(scan)
    for sample_scans in scans_by_sample.values():
        sample_scans.sort(key=lambda scan: float(scan.rt))

    diagnostic_cache: dict[str, tuple[list[float], list[float], dict[str, object]]] = {}
    winners: list[dict[str, object]] = []
    for feature in features:
        feature_id = str(feature.get("feature_group_id") or "")
        feature_mz = float(
            feature.get("true_peak_mz")
            or feature.get("representative_mz")
            or 0.0
        )
        feature_rt = float(feature.get("representative_rt") or 0.0)
        local_shifts = dict(feature.get("rt_correction_by_sample") or {})
        selected_scans: list[tuple[LCMSSpectrumScan, int]] = []
        for sample_id, sample_scans in scans_by_sample.items():
            center_rt = (
                feature_rt
                - float(shifts.get(sample_id) or 0.0)
                - float(local_shifts.get(sample_id) or 0.0)
            )
            possible: list[tuple[tuple[float, float], LCMSSpectrumScan, int]] = []
            for scan in sample_scans:
                if abs(float(scan.rt) - center_rt) > max(0.0, rt_tolerance_min):
                    continue
                charge = int(scan.precursor_charge or 0)
                if charge <= 0:
                    continue
                aligned_rt = (
                    float(scan.rt)
                    + float(shifts.get(sample_id) or 0.0)
                    + float(local_shifts.get(sample_id) or 0.0)
                )
                relation = _scan_isotope_relation(
                    feature_mz,
                    float(scan.precursor_mz or 0.0),
                    charge,
                    max(
                        feature_mz * max(0.1, precursor_tolerance_ppm)
                        / 1_000_000.0,
                        0.02,
                    ),
                )
                if relation is None:
                    continue
                isotope_offset, mz_error = relation
                possible.append((
                    (abs(aligned_rt - feature_rt), abs(float(mz_error))),
                    scan,
                    isotope_offset,
                ))
            selected_scans.extend(
                (scan, isotope_offset)
                for _, scan, isotope_offset in sorted(
                    possible,
                    key=lambda item: item[0],
                )[:max(1, max_scans_per_sample)]
            )

        for scan, isotope_offset in selected_scans:
            cached = diagnostic_cache.get(str(scan.scan_id))
            if cached is None:
                observed_pairs = sorted(
                    zip(scan.mz_array, scan.intensity_array),
                    key=lambda pair: pair[0],
                )
                observed_mz = [float(pair[0]) for pair in observed_pairs]
                observed_intensity = [
                    max(0.0, float(pair[1])) for pair in observed_pairs
                ]
                diagnostic = _targeted_glycan_diagnostic_evidence(
                    observed_mz,
                    observed_intensity,
                    fragment_tolerance_ppm,
                )
                cached = (observed_mz, observed_intensity, diagnostic)
                diagnostic_cache[str(scan.scan_id)] = cached
            observed_mz, observed_intensity, diagnostic = cached
            if (
                not diagnostic["glycan_has_hexnac_204"]
                or int(diagnostic["glycan_diagnostic_ion_count"]) < 2
            ):
                continue
            charge = int(scan.precursor_charge or 0)
            component_mass = (
                (feature_mz - PROTON) * charge
                - isotope_offset * ISOTOPE_MASS_DIFF
            )
            tolerance = (
                component_mass * max(0.1, precursor_tolerance_ppm)
                / 1_000_000.0
            )
            left = bisect_left(glycopeptide_masses, component_mass - tolerance)
            right = bisect_right(glycopeptide_masses, component_mass + tolerance)
            evaluated: list[dict[str, object]] = []
            for candidate in glycopeptides[left:right]:
                site_index, glycan_name, glycan_mass = candidate.modifications[0]
                glycan = glycan_by_name[glycan_name]
                composition = dict(glycan["composition"])
                if (
                    int(composition.get("NeuAc") or 0) > 0
                    and not diagnostic["glycan_has_sialic_diagnostic"]
                ):
                    continue
                backbone = _candidate(
                    candidate.chain,
                    candidate.start,
                    candidate.end,
                    candidate.sequence,
                    carbamidomethyl_cys=candidate.carbamidomethyl_cys,
                    proteolysis=candidate.proteolysis,
                    is_decoy=candidate.is_decoy,
                    base_peptide_id=candidate.base_peptide_id,
                )
                core_y = _targeted_glycan_core_y_evidence(
                    backbone,
                    composition,
                    observed_mz,
                    observed_intensity,
                    charge,
                    fragment_tolerance_ppm,
                )
                if int(core_y["glycan_core_y_ion_count"]) < 2:
                    continue
                backbone_evidence = _targeted_glycan_backbone_evidence(
                    backbone,
                    site_index,
                    observed_mz,
                    observed_intensity,
                    fragment_tolerance_ppm,
                )
                mass_error_ppm = (
                    (component_mass - candidate.neutral_mass)
                    / max(candidate.neutral_mass, 1e-12)
                    * 1_000_000.0
                )
                mass_score = max(
                    0.0,
                    1.0
                    - abs(mass_error_ppm)
                    / max(precursor_tolerance_ppm, 1e-9),
                )
                diagnostic_score = min(
                    18.0,
                    3.0 * int(diagnostic["glycan_diagnostic_ion_count"])
                    + 20.0
                    * float(diagnostic["glycan_diagnostic_intensity_fraction"]),
                )
                core_y_score = min(
                    25.0,
                    5.0 * int(core_y["glycan_core_y_ion_count"]),
                )
                hexnac_score = min(
                    15.0,
                    5.0
                    * int(backbone_evidence["glycan_hexnac_fragment_count"]),
                )
                score = (
                    15.0 * mass_score
                    + diagnostic_score
                    + core_y_score
                    + hexnac_score
                    + 0.55 * float(backbone_evidence["score"])
                )
                label_matches = {
                    int(match["observed_index"]): str(match["label"])
                    for match in (
                        list(diagnostic["glycan_diagnostic_matches"])
                        + list(core_y["glycan_core_y_matches"])
                        + list(backbone_evidence["matched_fragments"])
                    )
                }
                top_indices = sorted(
                    range(len(observed_mz)),
                    key=observed_intensity.__getitem__,
                    reverse=True,
                )[:80]
                spectrum_peaks = [
                    {
                        "mz": observed_mz[index],
                        "intensity": observed_intensity[index],
                        "label": label_matches.get(index, ""),
                    }
                    for index in sorted(top_indices, key=observed_mz.__getitem__)
                ]
                link = _feature_link_payload(
                    feature,
                    "feature_targeted_glycopeptide_search",
                    isotope_offset,
                )
                row = {
                    **asdict(candidate),
                    "modification_text": candidate.modification_text,
                    "sample_id": scan.sample_id,
                    "scan_id": scan.scan_id,
                    "rt": scan.rt,
                    "precursor_mz": scan.precursor_mz,
                    "precursor_charge": charge,
                    "precursor_intensity": scan.precursor_intensity or 0.0,
                    "precursor_error_ppm": mass_error_ppm,
                    "activation_method": scan.activation_method,
                    "collision_energy": scan.collision_energy,
                    "spectrum_peaks": spectrum_peaks,
                    **backbone_evidence,
                    **diagnostic,
                    **core_y,
                    "score": score,
                    "mass_delta": float(glycan_mass),
                    "localized_mass_offset_site": site_index + 1,
                    "glycan_name": glycan_name,
                    "glycan_composition": composition,
                    "glycan_site": f"{candidate.sequence[site_index]}{site_index + 1}",
                    "search_origin": "feature_targeted_glycopeptide_search",
                    "feature_links": [link],
                    **{key: value for key, value in link.items() if key != "feature_link_distance"},
                }
                evaluated.append(row)
            if not evaluated:
                continue
            evaluated.sort(key=lambda row: float(row["score"]), reverse=True)
            winner = evaluated[0]
            competitor_scores = [
                float(row["score"])
                for row in evaluated[1:]
                if (
                    row.get("candidate_id") != winner.get("candidate_id")
                    or row.get("is_decoy") != winner.get("is_decoy")
                )
            ]
            winner["candidate_score_margin"] = (
                float(winner["score"]) - max(competitor_scores)
                if competitor_scores else float(winner["score"])
            )
            winners.append(winner)

    searched = assign_q_values(winners)
    accepted: list[dict[str, object]] = []
    for row in searched:
        if row.get("is_decoy"):
            continue
        if float(row.get("q_value") if row.get("q_value") is not None else 1.0) > fdr_threshold:
            continue
        if float(row.get("score") or 0.0) < 48.0:
            continue
        if int(row.get("glycan_diagnostic_ion_count") or 0) < 2:
            continue
        if int(row.get("glycan_core_y_ion_count") or 0) < 2:
            continue
        if int(row.get("matched_ion_count") or 0) < 3:
            continue
        if not (
            int(row.get("sequence_tag_length") or 0) >= 1
            or int(row.get("complementary_ion_pair_count") or 0) >= 1
            or int(row.get("glycan_hexnac_fragment_count") or 0) >= 1
        ):
            continue
        if float(row.get("candidate_score_margin") or 0.0) < 2.0:
            continue
        row["feature_link_type"] = "feature_targeted_glycopeptide_search"
        accepted.append(row)
    return accepted


def search_feature_open_mass_scans(
    payload: dict[str, object],
    scans: list[LCMSSpectrumScan],
    backbone_candidates: list[PeptideCandidate],
    exclude_feature_ids: set[str] | None = None,
    rt_tolerance_min: float = 0.5,
    precursor_tolerance_ppm: float = 20.0,
    fragment_tolerance_ppm: float = 20.0,
    minimum_mass_delta_da: float = -250.0,
    maximum_mass_delta_da: float = 2500.0,
    minimum_absolute_delta_da: float = 0.5,
    max_scans_per_sample: int = 2,
    max_screened_candidates: int = 12,
    fdr_threshold: float = 0.01,
) -> list[dict[str, object]]:
    """Search selected precursor/isotope scans with one unrestricted mass delta.

    This is deliberately feature-centric: only a Feature with an actually
    selected precursor or selected isotope MS2 scan is considered.  Features
    covered only by an isolation window and Features without any MS2 are not
    sent through this search.
    """
    excluded = {str(value) for value in (exclude_feature_ids or set())}
    significant_types = {"presence_absence", "area_changed", "moderate_difference"}
    features = [
        row for row in payload.get("global_feature_groups") or []
        if row.get("difference_type") in significant_types
        and str(row.get("feature_group_id") or "") not in excluded
    ]
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
    for scan in scans:
        if scan.ms_level == 2 and scan.precursor_mz is not None and scan.mz_array:
            scans_by_sample.setdefault(str(scan.sample_id), []).append(scan)

    for sample_scans in scans_by_sample.values():
        sample_scans.sort(key=lambda scan: float(scan.rt))

    candidates = [candidate for candidate in backbone_candidates if len(candidate.sequence) >= 6]
    ordered = sorted(candidates, key=lambda candidate: candidate.neutral_mass)
    winners: list[dict[str, object]] = []

    for feature in features:
        feature_id = str(feature.get("feature_group_id") or "")
        feature_mz = float(
            feature.get("true_peak_mz")
            or feature.get("representative_mz")
            or 0.0
        )
        feature_rt = float(feature.get("representative_rt") or 0.0)
        local_shifts = dict(feature.get("rt_correction_by_sample") or {})
        selected_by_sample: dict[str, list[tuple[tuple[float, float], LCMSSpectrumScan, int]]] = {}
        for sample_id, sample_scans in scans_by_sample.items():
            possible: list[tuple[tuple[float, float], LCMSSpectrumScan, int]] = []
            center_rt = feature_rt - float(shifts.get(sample_id) or 0.0) - float(local_shifts.get(sample_id) or 0.0)
            candidate_scans = [
                scan for scan in sample_scans
                if abs(float(scan.rt) - center_rt) <= max(0.0, rt_tolerance_min)
            ]
            for scan in candidate_scans:
                aligned_rt = (
                    float(scan.rt)
                    + float(shifts.get(sample_id) or 0.0)
                    + float(local_shifts.get(sample_id) or 0.0)
                )
                rt_error = aligned_rt - feature_rt
                if abs(rt_error) > max(0.0, rt_tolerance_min):
                    continue
                relation_tolerance = max(
                    feature_mz * max(precursor_tolerance_ppm, 0.1) / 1_000_000.0,
                    0.02,
                )
                relation = _scan_isotope_relation(
                    feature_mz,
                    float(scan.precursor_mz or 0.0),
                    scan.precursor_charge,
                    relation_tolerance,
                )
                if relation is None:
                    continue
                isotope_offset, mz_error = relation
                possible.append((
                    (
                        abs(float(isotope_offset)),
                        abs(float(rt_error)) + abs(float(mz_error)),
                    ),
                    scan,
                    isotope_offset,
                ))
            if possible:
                selected_by_sample[sample_id] = sorted(possible, key=lambda item: item[0])[:max(1, max_scans_per_sample)]

        for sample_id, selected in selected_by_sample.items():
            for _, scan, isotope_offset in selected:
                charge = int(scan.precursor_charge or 0)
                if charge <= 0:
                    continue
                component_mass = (feature_mz - PROTON) * charge - isotope_offset * ISOTOPE_MASS_DIFF
                if component_mass <= 0.0:
                    continue
                local_candidates = [
                    candidate for candidate in ordered
                    if minimum_mass_delta_da
                    <= component_mass - candidate.neutral_mass
                    <= maximum_mass_delta_da
                    and (
                        abs(component_mass - candidate.neutral_mass) >= minimum_absolute_delta_da
                        or candidate.proteolysis in {
                            "n_terminal_truncation_candidate",
                            "c_terminal_truncation_candidate",
                        }
                    )
                ]
                # Open search is a rescue path.  Limit expensive fragment
                # matching to the closest precursor-mass candidates first.
                local_candidates.sort(
                    key=lambda candidate: abs(component_mass - candidate.neutral_mass)
                )
                local_candidates = local_candidates[: max(48, max_screened_candidates * 2)]
                optimistic: list[tuple[float, PeptideCandidate, float]] = []
                observed_pairs = sorted(
                    zip(scan.mz_array, scan.intensity_array), key=lambda pair: pair[0]
                )
                observed_mz = [float(pair[0]) for pair in observed_pairs]
                observed_intensity = [max(0.0, float(pair[1])) for pair in observed_pairs]
                for candidate in local_candidates:
                    delta = component_mass - candidate.neutral_mass
                    screen = _screen_mass_offset_fragments(
                        candidate,
                        delta,
                        observed_mz,
                        observed_intensity,
                        fragment_tolerance_ppm,
                    )
                    if int(screen.get("matched_ion_count") or 0) < 3:
                        continue
                    if not (
                        int(screen.get("sequence_tag_length") or 0) >= 1
                        or int(screen.get("complementary_ion_pair_count") or 0) >= 1
                    ):
                        continue
                    optimistic.append((float(screen.get("score") or 0.0), candidate, delta))
                if not optimistic:
                    continue
                evaluated: list[dict[str, object]] = []
                for _, candidate, delta in sorted(optimistic, key=lambda item: item[0], reverse=True)[:max(1, max_screened_candidates)]:
                    evidence = _match_mass_offset_fragments(
                        candidate,
                        delta,
                        observed_mz,
                        observed_intensity,
                        fragment_tolerance_ppm,
                    )
                    if int(evidence.get("matched_ion_count") or 0) < 4:
                        continue
                    localized_site = evidence.get("localized_mass_offset_site")
                    modification_description, compatible_offsets = (
                        _mass_offset_description(
                            candidate,
                            delta,
                            localized_site,
                        )
                    )
                    display_modification_text = _mass_offset_display_label(
                        candidate,
                        delta,
                        localized_site,
                        compatible_offsets,
                    )
                    labels = {
                        int(match["observed_index"]): str(match["label"])
                        for match in evidence.get("matched_fragments") or []
                    }
                    top_indices = sorted(
                        range(len(observed_mz)),
                        key=observed_intensity.__getitem__,
                        reverse=True,
                    )[:60]
                    spectrum_peaks = [
                        {
                            "mz": observed_mz[index],
                            "intensity": observed_intensity[index],
                            "label": labels.get(index, ""),
                        }
                        for index in sorted(
                            top_indices,
                            key=observed_mz.__getitem__,
                        )
                    ]
                    row = {
                        **asdict(candidate),
                        "candidate_id": f"{candidate.candidate_id}:OPEN_DELTA={delta:+.4f}",
                        "neutral_mass": component_mass,
                        "backbone_neutral_mass": candidate.neutral_mass,
                        "modification_text": display_modification_text,
                        "mass_offset_description": modification_description,
                        "compatible_mass_offsets": compatible_offsets,
                        "mass_delta": delta,
                        "localized_mass_offset_site": localized_site,
                        "sample_id": sample_id,
                        "scan_id": scan.scan_id,
                        "rt": scan.rt,
                        "precursor_mz": scan.precursor_mz,
                        "precursor_charge": charge,
                        "precursor_intensity": scan.precursor_intensity or 0.0,
                        "precursor_error_ppm": 0.0,
                        "activation_method": scan.activation_method,
                        "collision_energy": scan.collision_energy,
                        "spectrum_peaks": spectrum_peaks,
                        **evidence,
                        "search_origin": "feature_open_mass_search",
                        "feature_links": [
                            _feature_link_payload(feature, "feature_open_mass_search", isotope_offset)
                        ],
                    }
                    evaluated.append(row)
                if not evaluated:
                    continue
                evaluated.sort(key=lambda row: float(row.get("score") or 0.0), reverse=True)
                winner = evaluated[0]
                competitors = [
                    float(row.get("score") or 0.0)
                    for row in evaluated[1:]
                    if (
                        row.get("sequence") != winner.get("sequence")
                        or row.get("chain") != winner.get("chain")
                        or row.get("start") != winner.get("start")
                    )
                ]
                winner["candidate_score_margin"] = (
                    float(winner.get("score") or 0.0) - max(competitors)
                    if competitors else float(winner.get("score") or 0.0)
                )
                winner["feature_link_type"] = "feature_open_mass_search"
                winner["feature_isotope_offset"] = isotope_offset
                winner.update(winner["feature_links"][0])
                winners.append(winner)

    searched = assign_q_values(winners)
    accepted: list[dict[str, object]] = []
    for row in searched:
        if row.get("is_decoy"):
            continue
        q_value = float(row.get("q_value") if row.get("q_value") is not None else 1.0)
        if q_value > fdr_threshold:
            continue
        if float(row.get("score") or 0.0) < 42.0:
            continue
        if int(row.get("matched_ion_count") or 0) < 5:
            continue
        if int(row.get("sequence_tag_length") or 0) < 2 and int(row.get("complementary_ion_pair_count") or 0) < 1:
            continue
        if float(row.get("candidate_score_margin") or 0.0) < 1.5:
            continue
        accepted.append(row)
    return accepted


def _cluster_open_mass_results(
    rows: list[dict[str, object]],
    tolerance_da: float = 0.02,
) -> list[dict[str, object]]:
    """Cluster recurrent open-search mass shifts and estimate cluster q-values."""
    clusters: list[dict[str, object]] = []
    for row in sorted(rows, key=lambda item: float(item.get("mass_delta") or 0.0)):
        delta = float(row.get("mass_delta") or 0.0)
        cluster = clusters[-1] if clusters else None
        if cluster is None or abs(delta - float(cluster["mass_delta"])) > tolerance_da:
            cluster = {"mass_delta": delta, "rows": []}
            clusters.append(cluster)
        cluster_rows = cluster["rows"]
        assert isinstance(cluster_rows, list)
        cluster_rows.append(row)
        cluster["mass_delta"] = statistics.mean(
            float(item.get("mass_delta") or 0.0) for item in cluster_rows
        )

    summaries: list[dict[str, object]] = []
    for cluster in clusters:
        cluster_rows = list(cluster["rows"])
        targets = [row for row in cluster_rows if not row.get("is_decoy")]
        decoys = [row for row in cluster_rows if row.get("is_decoy")]
        scored_rows = targets or decoys
        interpretations = sorted({
            str(value)
            for row in targets
            for value in row.get("compatible_mass_offsets") or []
            if value
        })
        best_score = max(float(row.get("score") or 0.0) for row in scored_rows)
        summaries.append({
            "mass_delta": float(cluster["mass_delta"]),
            "mass_delta_min": min(float(row.get("mass_delta") or 0.0) for row in cluster_rows),
            "mass_delta_max": max(float(row.get("mass_delta") or 0.0) for row in cluster_rows),
            "target_psm_count": len(targets),
            "decoy_psm_count": len(decoys),
            "sample_count": len({str(row.get("sample_id") or "") for row in targets}),
            "samples": sorted({str(row.get("sample_id") or "") for row in targets}),
            "sequence_count": len({str(row.get("sequence") or "") for row in targets}),
            "sequences": sorted({str(row.get("sequence") or "") for row in targets})[:12],
            "localized_site_count": sum(row.get("localized_mass_offset_site") is not None for row in targets),
            "candidate_interpretations": interpretations,
            "best_score": best_score,
            "cluster_score": best_score + 4.0 * math.log1p(max(len(targets), len(decoys))),
            "rows": cluster_rows,
        })
    summaries.sort(key=lambda item: float(item["cluster_score"]), reverse=True)
    running_targets = 0
    running_decoys = 0
    fdrs: list[float] = []
    for summary in summaries:
        running_targets += int(summary["target_psm_count"])
        running_decoys += int(summary["decoy_psm_count"])
        fdrs.append(running_decoys / max(1, running_targets))
    q_value = 1.0
    for index in range(len(summaries) - 1, -1, -1):
        q_value = min(q_value, fdrs[index])
        summaries[index]["cluster_q_value"] = q_value
    for rank, summary in enumerate(summaries, start=1):
        cluster_id = f"OPEN_DELTA_{rank:04d}"
        summary["cluster_id"] = cluster_id
        for row in summary.pop("rows"):
            row["global_delta_cluster_id"] = cluster_id
            row["global_delta_cluster_mass"] = summary["mass_delta"]
            row["global_delta_cluster_q_value"] = summary["cluster_q_value"]
            row["global_delta_cluster_support"] = summary["target_psm_count"]
    return [
        summary for summary in summaries
        if int(summary.get("target_psm_count") or 0) > 0
    ]


def search_global_open_modification_scans(
    scans: list[LCMSSpectrumScan],
    backbone_candidates: list[PeptideCandidate],
    exclude_scan_ids: set[str] | None = None,
    fragment_tolerance_ppm: float = 20.0,
    minimum_mass_delta_da: float = -250.0,
    maximum_mass_delta_da: float = 2500.0,
    minimum_absolute_delta_da: float = 0.5,
    delta_cluster_tolerance_da: float = 0.02,
    max_discovery_scans: int = 15000,
    max_screened_candidates: int = 12,
    max_delta_clusters: int = 40,
    fdr_threshold: float = 0.01,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Run a two-pass unrestricted modification search over acquired MS2.

    Pass one discovers recurrent mass shifts on representative unresolved
    precursor scans using a fragment index.  Pass two re-searches every
    unresolved acquired MS2 scan against only the discovered delta clusters.
    Features without an acquired MS2 scan never enter this workflow.
    """
    excluded = {str(value) for value in (exclude_scan_ids or set())}
    valid_scans = [
        scan for scan in scans
        if scan.ms_level == 2
        and scan.precursor_mz is not None
        and int(scan.precursor_charge or 0) > 0
        and bool(scan.mz_array)
        and str(scan.scan_id) not in excluded
    ]
    candidates = [candidate for candidate in backbone_candidates if len(candidate.sequence) >= 6]
    if not valid_scans or not candidates:
        return [], [], []

    # A coarse fragment index supplies sequence-backbone candidates without
    # constraining precursor delta.  Exact ppm matching is still performed by
    # the regular open-mass matcher after this inexpensive screen.
    fragment_bin_width = 0.02
    fragment_index: dict[int, set[int]] = {}
    for candidate_index, candidate in enumerate(candidates):
        for _, fragment_mz, series, __ in theoretical_fragments(candidate):
            if series not in {"b", "y"} or fragment_mz <= 0.0:
                continue
            fragment_index.setdefault(int(round(fragment_mz / fragment_bin_width)), set()).add(candidate_index)

    def scan_quality(scan: LCMSSpectrumScan) -> float:
        strongest = sorted((max(0.0, float(value)) for value in scan.intensity_array), reverse=True)[:20]
        return math.log1p(float(scan.precursor_intensity or 0.0)) + math.log1p(sum(strongest))

    representative_by_bin: dict[tuple[str, int, int, int], LCMSSpectrumScan] = {}
    for scan in valid_scans:
        charge = int(scan.precursor_charge or 0)
        neutral_mass = (float(scan.precursor_mz or 0.0) - PROTON) * charge
        key = (
            str(scan.sample_id),
            charge,
            int(round(neutral_mass / 0.05)),
            int(round(float(scan.rt) / 0.20)),
        )
        previous = representative_by_bin.get(key)
        if previous is None or scan_quality(scan) > scan_quality(previous):
            representative_by_bin[key] = scan
    discovery_scans = sorted(
        representative_by_bin.values(),
        key=scan_quality,
        reverse=True,
    )[:max(1, int(max_discovery_scans))]

    def candidate_indices_for_scan(scan: LCMSSpectrumScan) -> list[int]:
        charge = int(scan.precursor_charge or 0)
        component_mass = (float(scan.precursor_mz or 0.0) - PROTON) * charge
        selected_peaks = sorted(
            zip(scan.mz_array, scan.intensity_array),
            key=lambda pair: float(pair[1]),
            reverse=True,
        )[:60]
        hits: dict[int, set[int]] = {}
        for observed_mz, __ in selected_peaks:
            center = int(round(float(observed_mz) / fragment_bin_width))
            for fragment_bin in range(center - 2, center + 3):
                for candidate_index in fragment_index.get(fragment_bin, set()):
                    candidate = candidates[candidate_index]
                    delta = component_mass - candidate.neutral_mass
                    if not minimum_mass_delta_da <= delta <= maximum_mass_delta_da:
                        continue
                    if (
                        abs(delta) < minimum_absolute_delta_da
                        and candidate.proteolysis not in {
                            "n_terminal_truncation_candidate",
                            "c_terminal_truncation_candidate",
                        }
                    ):
                        continue
                    hits.setdefault(candidate_index, set()).add(fragment_bin)
        ranked = sorted(
            hits,
            key=lambda index: (
                len(hits[index]),
                -abs(component_mass - candidates[index].neutral_mass),
            ),
            reverse=True,
        )
        return ranked[:max(24, int(max_screened_candidates) * 4)]

    def evaluate_scan(
        scan: LCMSSpectrumScan,
        candidate_indices: list[int],
        cluster: dict[str, object] | None = None,
    ) -> dict[str, object] | None:
        charge = int(scan.precursor_charge or 0)
        component_mass = (float(scan.precursor_mz or 0.0) - PROTON) * charge
        observed_pairs = sorted(zip(scan.mz_array, scan.intensity_array), key=lambda pair: pair[0])
        observed_mz = [float(pair[0]) for pair in observed_pairs]
        observed_intensity = [max(0.0, float(pair[1])) for pair in observed_pairs]
        optimistic: list[tuple[float, int, float]] = []
        for candidate_index in candidate_indices:
            candidate = candidates[candidate_index]
            delta = component_mass - candidate.neutral_mass
            if not minimum_mass_delta_da <= delta <= maximum_mass_delta_da:
                continue
            screen = _screen_mass_offset_fragments(
                candidate,
                delta,
                observed_mz,
                observed_intensity,
                fragment_tolerance_ppm,
            )
            if int(screen.get("matched_ion_count") or 0) < 3:
                continue
            if int(screen.get("sequence_tag_length") or 0) < 1 and int(screen.get("complementary_ion_pair_count") or 0) < 1:
                continue
            optimistic.append((float(screen.get("score") or 0.0), candidate_index, delta))
        evaluated: list[dict[str, object]] = []
        for _, candidate_index, delta in sorted(optimistic, reverse=True)[:max(1, int(max_screened_candidates))]:
            candidate = candidates[candidate_index]
            evidence = _match_mass_offset_fragments(
                candidate,
                delta,
                observed_mz,
                observed_intensity,
                fragment_tolerance_ppm,
            )
            if int(evidence.get("matched_ion_count") or 0) < 4:
                continue
            localized_site = evidence.get("localized_mass_offset_site")
            if localized_site is None:
                unique_site = _unique_compatible_mass_offset_site(
                    candidate,
                    delta,
                )
                if (
                    unique_site is not None
                    and int(evidence.get("regular_fragment_count") or 0) > 0
                    and int(evidence.get("shifted_fragment_count") or 0) > 0
                ):
                    localized_site = unique_site
                    evidence["localized_mass_offset_site"] = unique_site
                    evidence["mass_offset_localization_basis"] = (
                        "unique_compatible_residue_with_fragment_boundary_support"
                    )
            modification_description, compatible_offsets = _mass_offset_description(
                candidate,
                delta,
                localized_site,
            )
            display_modification_text = _mass_offset_display_label(
                candidate,
                delta,
                localized_site,
                compatible_offsets,
            )
            model_prior = 2.0 if compatible_offsets else 0.0
            evidence["score"] = float(evidence.get("score") or 0.0) + model_prior
            labels = {
                int(match["observed_index"]): str(match["label"])
                for match in evidence.get("matched_fragments") or []
            }
            top_indices = sorted(
                range(len(observed_mz)),
                key=observed_intensity.__getitem__,
                reverse=True,
            )[:60]
            row = {
                **asdict(candidate),
                "candidate_id": f"{candidate.candidate_id}:GLOBAL_OPEN_DELTA={delta:+.4f}",
                "neutral_mass": component_mass,
                "backbone_neutral_mass": candidate.neutral_mass,
                "modification_text": display_modification_text,
                "mass_offset_description": modification_description,
                "compatible_mass_offsets": compatible_offsets,
                "candidate_model_prior": model_prior,
                "mass_delta": delta,
                "localized_mass_offset_site": localized_site,
                "sample_id": scan.sample_id,
                "scan_id": scan.scan_id,
                "rt": scan.rt,
                "precursor_mz": scan.precursor_mz,
                "precursor_charge": charge,
                "precursor_intensity": scan.precursor_intensity or 0.0,
                "precursor_error_ppm": 0.0,
                "activation_method": scan.activation_method,
                "collision_energy": scan.collision_energy,
                "spectrum_peaks": [
                    {
                        "mz": observed_mz[index],
                        "intensity": observed_intensity[index],
                        "label": labels.get(index, ""),
                    }
                    for index in sorted(top_indices, key=observed_mz.__getitem__)
                ],
                **evidence,
                "search_origin": "global_open_modification_search",
            }
            if cluster is not None:
                row["discovery_delta_cluster_mass"] = cluster.get("mass_delta")
            evaluated.append(row)
        if not evaluated:
            return None
        evaluated.sort(key=lambda row: float(row.get("score") or 0.0), reverse=True)
        winner = evaluated[0]
        competing_scores = [
            float(row.get("score") or 0.0)
            for row in evaluated[1:]
            if (
                row.get("sequence") != winner.get("sequence")
                or row.get("chain") != winner.get("chain")
                or row.get("start") != winner.get("start")
            )
        ]
        winner["candidate_score_margin"] = (
            float(winner.get("score") or 0.0) - max(competing_scores)
            if competing_scores else float(winner.get("score") or 0.0)
        )
        return winner

    discovery_rows = [
        row for row in (
            evaluate_scan(scan, candidate_indices_for_scan(scan))
            for scan in discovery_scans
        )
        if row is not None
    ]
    discovery_rows = assign_q_values(discovery_rows)
    discovery_clusters = _cluster_open_mass_results(
        discovery_rows,
        delta_cluster_tolerance_da,
    )
    selected_clusters = [
        cluster for cluster in discovery_clusters
        if (
            int(cluster.get("target_psm_count") or 0) >= 2
            or bool(cluster.get("candidate_interpretations"))
        )
        and float(cluster.get("best_score") or 0.0) >= 34.0
        and float(
            cluster["cluster_q_value"]
            if cluster.get("cluster_q_value") is not None
            else 1.0
        ) <= 0.20
    ][:max(1, int(max_delta_clusters))]
    if not selected_clusters:
        return [], [], discovery_clusters

    ordered_mass_candidates = sorted(enumerate(candidates), key=lambda item: item[1].neutral_mass)
    candidate_masses = [candidate.neutral_mass for _, candidate in ordered_mass_candidates]

    def restricted_candidates(scan: LCMSSpectrumScan, cluster: dict[str, object]) -> list[int]:
        charge = int(scan.precursor_charge or 0)
        component_mass = (float(scan.precursor_mz or 0.0) - PROTON) * charge
        expected_backbone_mass = component_mass - float(cluster.get("mass_delta") or 0.0)
        tolerance = max(delta_cluster_tolerance_da, expected_backbone_mass * 20.0 / 1_000_000.0)
        left = bisect_left(candidate_masses, expected_backbone_mass - tolerance)
        right = bisect_right(candidate_masses, expected_backbone_mass + tolerance)
        return [ordered_mass_candidates[index][0] for index in range(left, right)]

    second_pass_rows: list[dict[str, object]] = []
    for scan in valid_scans:
        winners = [
            row for row in (
                evaluate_scan(scan, restricted_candidates(scan, cluster), cluster)
                for cluster in selected_clusters
            )
            if row is not None
        ]
        if winners:
            second_pass_rows.append(max(winners, key=lambda row: float(row.get("score") or 0.0)))
    searched = assign_q_values(second_pass_rows)
    final_clusters = _cluster_open_mass_results(searched, delta_cluster_tolerance_da)
    cluster_by_id = {str(cluster["cluster_id"]): cluster for cluster in final_clusters}
    accepted: list[dict[str, object]] = []
    candidates_out: list[dict[str, object]] = []
    for row in searched:
        if row.get("is_decoy"):
            continue
        cluster = cluster_by_id.get(str(row.get("global_delta_cluster_id") or "")) or {}
        q_value = float(row.get("q_value") if row.get("q_value") is not None else 1.0)
        cluster_q = float(
            row["global_delta_cluster_q_value"]
            if row.get("global_delta_cluster_q_value") is not None
            else 1.0
        )
        score = float(row.get("score") or 0.0)
        ions = int(row.get("matched_ion_count") or 0)
        tag_length = int(row.get("sequence_tag_length") or 0)
        complementary = int(row.get("complementary_ion_pair_count") or 0)
        margin = float(row.get("candidate_score_margin") or 0.0)
        localized = row.get("localized_mass_offset_site") is not None
        strong = (
            q_value <= fdr_threshold
            and cluster_q <= max(0.05, fdr_threshold)
            and score >= 45.0
            and ions >= 5
            and (tag_length >= 2 or complementary >= 1)
            and margin >= 2.0
            and localized
            and int(cluster.get("target_psm_count") or 0) >= 2
        )
        tentative = (
            q_value <= max(0.05, fdr_threshold)
            and score >= 34.0
            and ions >= 4
            and (tag_length >= 1 or complementary >= 1)
            and margin >= 1.0
            and int(cluster.get("target_psm_count") or 0) >= 2
        )
        if strong:
            row["open_modification_confidence"] = "B_open_modification_supported"
            row["sequence_inference_level"] = "global_open_modification_supported"
            accepted.append(row)
        elif tentative:
            row["open_modification_confidence"] = "C_open_mass_candidate"
            row["sequence_inference_level"] = "global_open_mass_candidate"
            candidates_out.append(row)

    accepted_cluster_ids = {str(row.get("global_delta_cluster_id") or "") for row in accepted}
    candidate_cluster_ids = {str(row.get("global_delta_cluster_id") or "") for row in candidates_out}
    for cluster in final_clusters:
        cluster_id = str(cluster.get("cluster_id") or "")
        if cluster_id in accepted_cluster_ids:
            cluster["classification"] = "supported_open_modification"
        elif cluster_id in candidate_cluster_ids:
            cluster["classification"] = "open_mass_candidate"
        else:
            cluster["classification"] = "unresolved_mass_shift_cluster"
    return accepted, candidates_out, final_clusters


def build_global_known_modification_summary(
    supported_psms: list[dict[str, object]],
    candidate_psms: list[dict[str, object]],
    clusters: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Summarize chemically interpretable open-search results before unknown shifts."""
    grouped: dict[str, dict[str, object]] = {}

    def bucket(name: str) -> dict[str, object]:
        return grouped.setdefault(name, {
            "supported": [],
            "candidates": [],
            "clusters": [],
        })

    for cluster in clusters:
        for name in cluster.get("candidate_interpretations") or []:
            bucket(str(name))["clusters"].append(cluster)
    for level, rows in (("supported", supported_psms), ("candidates", candidate_psms)):
        for psm in rows:
            for name in psm.get("compatible_mass_offsets") or []:
                bucket(str(name))[level].append(psm)

    summary: list[dict[str, object]] = []
    for name, evidence in grouped.items():
        supported = list(evidence["supported"])
        candidates = list(evidence["candidates"])
        related_clusters = list(evidence["clusters"])
        all_psms = supported + candidates
        sample_evidence: dict[str, dict[str, object]] = {}
        for sample_id in sorted({str(row.get("sample_id") or "") for row in all_psms if row.get("sample_id")}):
            sample_supported = [row for row in supported if str(row.get("sample_id") or "") == sample_id]
            sample_candidates = [row for row in candidates if str(row.get("sample_id") or "") == sample_id]
            sample_rows = sample_supported + sample_candidates
            sample_evidence[sample_id] = {
                "supported_psm_count": len(sample_supported),
                "candidate_psm_count": len(sample_candidates),
                "sequence_count": len({str(row.get("sequence") or "") for row in sample_rows}),
            }

        localized_sites: set[str] = set()
        for row in all_psms:
            local_site = row.get("localized_mass_offset_site")
            if local_site is None:
                continue
            local_site = int(local_site)
            sequence = str(row.get("sequence") or "")
            residue = sequence[local_site - 1] if 1 <= local_site <= len(sequence) else "?"
            absolute_site = int(row.get("start") or 1) + local_site - 1
            localized_sites.add(f"{row.get('chain') or '?'}:{residue}{absolute_site}")

        delta_values = [float(row.get("mass_delta") or 0.0) for row in all_psms]
        if not delta_values:
            delta_values = [float(cluster.get("mass_delta") or 0.0) for cluster in related_clusters]
        best_psm = max(all_psms, key=lambda row: float(row.get("score") or 0.0), default=None)
        if supported:
            classification = "supported_known_modification"
        elif candidates:
            classification = "known_modification_candidate"
        else:
            classification = "unresolved_known_mass_shift"
        summary.append({
            "known_modification_id": "",
            "modification": name,
            "mass_delta": statistics.median(delta_values) if delta_values else None,
            "classification": classification,
            "supported_psm_count": len(supported),
            "candidate_psm_count": len(candidates),
            "target_winner_spectrum_count": sum(int(cluster.get("target_psm_count") or 0) for cluster in related_clusters),
            "sample_count": len(sample_evidence),
            "sample_evidence": sample_evidence,
            "sequence_count": len({str(row.get("sequence") or "") for row in all_psms}),
            "sequences": sorted({str(row.get("sequence") or "") for row in all_psms})[:20],
            "localized_site_count": len(localized_sites),
            "localized_sites": sorted(localized_sites)[:30],
            "cluster_ids": sorted({str(cluster.get("cluster_id") or "") for cluster in related_clusters}),
            "best_score": max((float(row.get("score") or 0.0) for row in all_psms), default=None),
            "best_psm_q_value": min((float(row["q_value"] if row.get("q_value") is not None else 1.0) for row in all_psms), default=None),
            "best_cluster_q_value": min((float(cluster["cluster_q_value"] if cluster.get("cluster_q_value") is not None else 1.0) for cluster in related_clusters), default=None),
            "best_psm": {
                key: best_psm.get(key) for key in (
                    "sample_id", "scan_id", "rt", "precursor_mz", "precursor_charge",
                    "chain", "start", "end", "sequence", "modification_text",
                    "mass_delta", "localized_mass_offset_site", "score", "q_value",
                    "global_delta_cluster_id", "spectrum_peaks",
                )
            } if best_psm is not None else None,
        })
    classification_order = {
        "supported_known_modification": 3,
        "known_modification_candidate": 2,
        "unresolved_known_mass_shift": 1,
    }
    summary.sort(key=lambda row: (
        classification_order.get(str(row.get("classification") or ""), 0),
        int(row.get("supported_psm_count") or 0),
        int(row.get("candidate_psm_count") or 0),
        int(row.get("target_winner_spectrum_count") or 0),
    ), reverse=True)
    for rank, row in enumerate(summary, start=1):
        row["rank"] = rank
        row["known_modification_id"] = f"KNOWN_MOD_{rank:03d}"
    return summary


def build_feature_ms2_evidence(
    payload: dict[str, object],
    scans: list[LCMSSpectrumScan],
    psms: list[dict[str, object]],
    annotations: list[dict[str, object]],
    candidates: list[PeptideCandidate] | None = None,
    sequence_region_candidates: list[dict[str, object]] | None = None,
    rt_tolerance_min: float = 0.5,
    mz_tolerance_ppm: float = 20.0,
    isolation_padding_da: float = 0.02,
) -> tuple[list[dict[str, object]], dict[str, int]]:
    features = sorted(
        list(payload.get("global_feature_groups") or []),
        key=lambda row: float(row.get("ranking_score") or 0.0),
        reverse=True,
    )
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    scans_by_sample: dict[str, list[LCMSSpectrumScan]] = {}
    for scan in scans:
        if scan.ms_level == 2 and scan.precursor_mz is not None:
            scans_by_sample.setdefault(scan.sample_id, []).append(scan)
    psms_by_feature: dict[str, list[dict[str, object]]] = {}
    for psm in psms:
        for link in psm.get("feature_links") or []:
            feature_id = str(link.get("feature_group_id") or "")
            if feature_id:
                psms_by_feature.setdefault(feature_id, []).append(psm)
    annotations_by_feature: dict[str, list[dict[str, object]]] = {}
    for annotation in annotations:
        feature_id = str(annotation.get("feature_or_candidate_id") or "")
        if feature_id:
            annotations_by_feature.setdefault(feature_id, []).append(annotation)
    region_candidates_by_feature: dict[str, list[dict[str, object]]] = {}
    for candidate in sequence_region_candidates or []:
        for link in candidate.get("feature_links") or []:
            feature_id = str(link.get("feature_group_id") or "")
            if feature_id:
                region_candidates_by_feature.setdefault(
                    feature_id,
                    [],
                ).append(candidate)

    ordered_targets = sorted(
        (candidate for candidate in candidates or [] if not candidate.is_decoy),
        key=lambda candidate: candidate.neutral_mass,
    )
    target_masses = [candidate.neutral_mass for candidate in ordered_targets]

    def mass_hypotheses(feature_mz: float, charge: int | None, isotope_offset: int) -> list[dict[str, object]]:
        matches: list[dict[str, object]] = []
        charges = [charge] if charge and charge > 0 else list(range(1, 8))
        for candidate_charge in charges:
            neutral_mass = (feature_mz - PROTON) * candidate_charge - isotope_offset * ISOTOPE_MASS_DIFF
            tolerance = neutral_mass * mz_tolerance_ppm / 1_000_000.0
            left = bisect_left(target_masses, neutral_mass - tolerance)
            right = bisect_right(target_masses, neutral_mass + tolerance)
            for candidate in ordered_targets[left:right]:
                error_ppm = (neutral_mass - candidate.neutral_mass) / candidate.neutral_mass * 1_000_000.0
                matches.append({
                    "candidate_id": candidate.candidate_id,
                    "sequence": candidate.sequence,
                    "modification": candidate.modification_text,
                    "proteolysis": candidate.proteolysis,
                    "charge": candidate_charge,
                    "neutral_mass": candidate.neutral_mass,
                    "mass_error_ppm": error_ppm,
                })
        return sorted(matches, key=lambda row: abs(float(row["mass_error_ppm"])))[:5]

    status_counts: dict[str, int] = {}
    evidence_rows: list[dict[str, object]] = []
    for rank, feature in enumerate(features, start=1):
        feature_id = str(feature.get("feature_group_id"))
        feature_mz = float(feature.get("representative_mz") or 0.0)
        feature_rt = float(feature.get("representative_rt") or 0.0)
        local_shifts = dict(feature.get("rt_correction_by_sample") or {})
        coverage_by_sample: dict[str, dict[str, object]] = {}
        best_coverage_scan: tuple[tuple[float, float, float], LCMSSpectrumScan, str, int] | None = None
        for sample_id, sample_scans in scans_by_sample.items():
            acquired = selected = isotope = 0
            best_sample_scan: tuple[tuple[float, float, float], LCMSSpectrumScan, str, int] | None = None
            for scan in sample_scans:
                aligned_rt = scan.rt + float(shifts.get(sample_id) or 0.0) + float(local_shifts.get(sample_id) or 0.0)
                rt_error = aligned_rt - feature_rt
                if abs(rt_error) > rt_tolerance_min:
                    continue
                precursor_mz = float(scan.precursor_mz or 0.0)
                lower = float(scan.isolation_window_lower_offset or 0.8)
                upper = float(scan.isolation_window_upper_offset or 0.8)
                if not precursor_mz - lower - isolation_padding_da <= feature_mz <= precursor_mz + upper + isolation_padding_da:
                    continue
                acquired += 1
                tolerance = max(feature_mz * mz_tolerance_ppm / 1_000_000.0, isolation_padding_da)
                relation = _scan_isotope_relation(feature_mz, precursor_mz, scan.precursor_charge, tolerance)
                relation_type = "isolation_window_only"
                isotope_offset = 0
                relation_priority = 2.0
                if relation is not None:
                    isotope_offset = relation[0]
                    if isotope_offset == 0:
                        selected += 1
                        relation_type = "selected_precursor"
                        relation_priority = 0.0
                    else:
                        isotope += 1
                        relation_type = "selected_isotope_envelope"
                        relation_priority = 1.0
                key = (relation_priority, abs(rt_error), abs(feature_mz - precursor_mz))
                candidate = (key, scan, relation_type, isotope_offset)
                if best_sample_scan is None or key < best_sample_scan[0]:
                    best_sample_scan = candidate
                if best_coverage_scan is None or key < best_coverage_scan[0]:
                    best_coverage_scan = candidate
            coverage_by_sample[sample_id] = {
                "isolation_window_scan_count": acquired,
                "selected_precursor_scan_count": selected,
                "selected_isotope_scan_count": isotope,
                "best_scan_id": best_sample_scan[1].scan_id if best_sample_scan else None,
                "best_relation": best_sample_scan[2] if best_sample_scan else "none",
            }

        feature_annotations = annotations_by_feature.get(feature_id) or []
        best_annotation = max(
            feature_annotations,
            key=lambda row: max((float(item.get("score") or 0.0) for item in dict(row.get("sample_evidence") or {}).values()), default=0.0),
            default=None,
        )
        feature_psms = psms_by_feature.get(feature_id) or []
        best_psm = max(feature_psms, key=lambda row: float(row.get("score") or 0.0), default=None)
        feature_region_candidates = sorted(
            region_candidates_by_feature.get(feature_id) or [],
            key=lambda row: float(row.get("score") or 0.0),
            reverse=True,
        )
        best_region_candidate = (
            feature_region_candidates[0]
            if feature_region_candidates else None
        )
        annotation_confidence = str(
            best_annotation.get("confidence") or ""
        ) if best_annotation is not None else ""
        if best_annotation is not None and annotation_confidence == "D_low_evidence":
            status = "low_evidence_sequence_candidate"
            hypothesis = (
                "MS2 存在序列匹配线索，但综合得分未达到当前暂定鉴定阈值；"
                "该序列仅作为低证据候选显示，不替代 Unknown 成分名称。"
            )
        elif best_annotation is not None:
            link_type = str(best_annotation.get("feature_link_type") or "direct_precursor")
            if link_type == "isotope_envelope":
                status = "identified_isotope_envelope"
            elif link_type == "charge_state_envelope":
                status = "identified_charge_state_envelope"
            elif link_type == "isolation_window_targeted_search":
                status = "tentative_feature_guided_identification"
            elif link_type == "selected_precursor_consensus_search":
                status = "tentative_feature_consensus_identification"
            elif link_type == "component_charge_isotope_consensus_search":
                status = "tentative_component_consensus_identification"
            elif link_type == "component_mass_offset_sequence_tag_search":
                status = (
                    "tentative_truncation_sequence_support"
                    if best_annotation.get("sequence_inference_level")
                    == "truncation_sequence_supported"
                    else "tentative_backbone_sequence_support"
                )
            elif link_type == "feature_open_mass_search":
                status = "tentative_feature_open_mass_identification"
            elif link_type == "feature_targeted_glycopeptide_search":
                status = "tentative_feature_glycopeptide_identification"
            elif link_type == "global_open_modification_search":
                status = "tentative_global_open_modification_identification"
            else:
                status = "identified_direct_precursor"
            hypothesis = str(best_annotation.get("interpretation") or "")
        elif best_region_candidate is not None:
            status = "tentative_sequence_region_candidate"
            hypothesis = (
                "组件中性质量和局部连续 b/y 标签支持该已知序列区段，"
                "但候选唯一性或碎片覆盖尚不足以升级为骨架序列定性。"
            )
        else:
            selected_count = sum(int(item["selected_precursor_scan_count"]) + int(item["selected_isotope_scan_count"]) for item in coverage_by_sample.values())
            acquired_count = sum(int(item["isolation_window_scan_count"]) for item in coverage_by_sample.values())
            if selected_count:
                status = "selected_precursor_unidentified"
                hypothesis = "已采集该前体/同位素包络的 MS2，但没有序列或修饰候选通过当前碎片证据阈值。"
            elif acquired_count:
                status = "coisolated_ms2_unresolved"
                hypothesis = "MS2 隔离窗覆盖该 MS1 Feature，但可能是共隔离/嵌合谱；当前不能可靠归属序列或修饰。"
            else:
                status = "no_ms2_acquired"
                hypothesis = "DDA 未采集覆盖该 Feature 的 MS2，无法仅凭 MS1 可靠推断具体修饰。"
        status_counts[status] = status_counts.get(status, 0) + 1
        unidentified_reason = None
        candidate_mass_hypotheses: list[dict[str, object]] = []
        selected_precursor_charge = None
        if best_annotation is not None and annotation_confidence == "D_low_evidence":
            unidentified_reason = "low_evidence_sequence_candidate_below_identification_threshold"
            selected_precursor_charge = (
                best_psm.get("precursor_charge") if best_psm is not None else None
            )
        elif best_annotation is None and best_region_candidate is not None:
            unidentified_reason = (
                "partial_sequence_tag_without_unique_identification"
            )
            selected_precursor_charge = best_region_candidate.get(
                "precursor_charge"
            )
        elif best_annotation is None:
            if status == "no_ms2_acquired":
                unidentified_reason = "no_ms2_scan"
            elif status == "coisolated_ms2_unresolved":
                unidentified_reason = "coisolated_not_selected"
            elif best_coverage_scan is not None:
                _, diagnostic_scan, diagnostic_relation, diagnostic_isotope_offset = best_coverage_scan
                selected_precursor_charge = diagnostic_scan.precursor_charge
                candidate_mass_hypotheses = mass_hypotheses(
                    feature_mz,
                    diagnostic_scan.precursor_charge,
                    diagnostic_isotope_offset,
                )
                if candidate_mass_hypotheses:
                    unidentified_reason = "candidate_mass_matched_but_fragments_insufficient"
                elif (
                    diagnostic_scan.precursor_charge == 1 and feature_mz < 650.0
                ) or (
                    diagnostic_scan.precursor_charge is None and feature_mz < 400.0
                ):
                    unidentified_reason = "likely_nonpeptide_or_small_component"
                elif diagnostic_relation in {"selected_precursor", "selected_isotope_envelope"}:
                    unidentified_reason = "outside_current_sequence_or_modification_space"
                else:
                    unidentified_reason = "coisolated_not_selected"
        if status == "low_evidence_sequence_candidate":
            peptide_likelihood = "medium"
            peptide_likelihood_reason = (
                "sequence_match_below_current_identification_threshold"
            )
        elif status in {
            "tentative_component_consensus_identification",
            "tentative_backbone_sequence_support",
            "tentative_truncation_sequence_support",
        }:
            peptide_likelihood = "medium"
            peptide_likelihood_reason = (
                "component_neutral_mass_and_joint_sequence_fragment_support"
            )
        elif status == "tentative_sequence_region_candidate":
            peptide_likelihood = "medium"
            peptide_likelihood_reason = (
                "partial_sequence_tag_and_component_mass_support"
            )
        elif best_annotation is not None:
            peptide_likelihood = "high"
            peptide_likelihood_reason = "sequence_informative_ms2"
        elif unidentified_reason == "candidate_mass_matched_but_fragments_insufficient":
            peptide_likelihood = "medium"
            peptide_likelihood_reason = "product_peptide_mass_candidate_without_sufficient_fragments"
        elif unidentified_reason == "likely_nonpeptide_or_small_component":
            peptide_likelihood = "low"
            peptide_likelihood_reason = "small_singly_charged_or_charge_unknown_without_product_mass_match"
        else:
            peptide_likelihood = "undetermined"
            peptide_likelihood_reason = "target_ms2_evidence_unavailable_or_search_space_mismatch"
        coverage_scan = None
        if best_coverage_scan is not None:
            _, scan, relation_type, isotope_offset = best_coverage_scan
            coverage_scan = {
                "sample_id": scan.sample_id,
                "scan_id": scan.scan_id,
                "rt": scan.rt,
                "precursor_mz": scan.precursor_mz,
                "precursor_charge": scan.precursor_charge,
                "relation_type": relation_type,
                "isotope_offset": isotope_offset,
                "activation_method": scan.activation_method,
                "collision_energy": scan.collision_energy,
                "spectrum_peaks": _top_spectrum_peaks(scan),
            }
        compact_best_psm = None
        if best_psm is not None:
            compact_best_psm = {
                key: best_psm.get(key)
                for key in (
                    "sample_id", "scan_id", "rt", "precursor_mz", "precursor_charge", "chain", "start", "end",
                    "sequence", "modification_text", "score", "q_value", "matched_ion_count", "fragment_coverage",
                    "explained_intensity", "feature_link_type", "feature_isotope_offset", "search_origin", "spectrum_peaks",
                    "feature_charge", "component_group_id", "component_neutral_mass",
                    "component_source_scan_count", "component_source_scan_ids",
                    "component_scan_support_count", "component_member_support_count",
                    "component_charge_support_count", "component_source_member_ids",
                    "component_source_charge_states",
                    "mass_delta", "backbone_neutral_mass",
                    "localized_mass_offset_site", "sequence_tag_length",
                    "b_ion_longest_run", "y_ion_longest_run",
                    "complementary_ion_pair_count", "candidate_score_margin",
                    "sequence_inference_level", "alternative_sequence_candidates",
                    "glycan_name", "glycan_composition", "glycan_site",
                    "glycan_diagnostic_ion_count",
                    "glycan_diagnostic_intensity_fraction",
                    "glycan_core_y_ion_count", "glycan_core_y_type_count",
                    "glycan_has_y1_support", "glycan_hexnac_fragment_count",
                )
            }
        compact_region_candidate = None
        if best_region_candidate is not None:
            compact_region_candidate = {
                key: best_region_candidate.get(key)
                for key in (
                    "sample_id", "scan_id", "rt", "precursor_mz",
                    "precursor_charge", "chain", "start", "end", "sequence",
                    "proteolysis", "modification_text", "score", "q_value",
                    "matched_ion_count", "fragment_coverage",
                    "explained_intensity", "feature_link_type",
                    "search_origin", "spectrum_peaks", "component_group_id",
                    "component_neutral_mass", "mass_delta",
                    "backbone_neutral_mass", "localized_mass_offset_site",
                    "sequence_tag_length", "b_ion_longest_run",
                    "y_ion_longest_run", "complementary_ion_pair_count",
                    "candidate_score_margin", "sequence_inference_level",
                    "alternative_sequence_candidates",
                )
            }
        evidence_rows.append({
            "rank": rank,
            "feature_group_id": feature_id,
            "parent_tic_peak_id": feature.get("parent_tic_peak_id"),
            "representative_rt": feature_rt,
            "representative_mz": feature_mz,
            "true_peak_mz": feature.get("true_peak_mz", feature_mz),
            "envelope_representative_mz": feature.get(
                "envelope_representative_mz",
                feature_mz,
            ),
            "difference_type": feature.get("difference_type"),
            "max_fold_change": feature.get("max_fold_change"),
            "ranking_score": feature.get("ranking_score"),
            "higher_abundance_sample": feature.get("higher_abundance_sample"),
            "lower_abundance_sample": feature.get("lower_abundance_sample"),
            "normalized_area_by_sample": feature.get("normalized_area_by_sample") or {},
            "ms2_status": status,
            "sequence": (
                best_annotation.get("sequence")
                if best_annotation else (
                    best_region_candidate.get("sequence")
                    if best_region_candidate else None
                )
            ),
            "modification": (
                best_annotation.get("modification")
                if best_annotation else (
                    best_region_candidate.get("modification_text")
                    if best_region_candidate else None
                )
            ),
            "confidence": (
                best_annotation.get("confidence")
                if best_annotation else (
                    "D_sequence_region_candidate"
                    if best_region_candidate else None
                )
            ),
            "feature_link_type": (
                best_annotation.get("feature_link_type")
                if best_annotation else (
                    "component_mass_offset_sequence_tag_candidate"
                    if best_region_candidate else None
                )
            ),
            "feature_isotope_offset": best_annotation.get("feature_isotope_offset") if best_annotation else None,
            "feature_charge": best_annotation.get("feature_charge") if best_annotation else None,
            "hypothesis": hypothesis,
            "unidentified_reason": unidentified_reason,
            "candidate_mass_hypotheses": candidate_mass_hypotheses,
            "sequence_region_candidates": [
                {
                    "sequence": row.get("sequence"),
                    "chain": row.get("chain"),
                    "start": row.get("start"),
                    "end": row.get("end"),
                    "proteolysis": row.get("proteolysis"),
                    "mass_delta": row.get("mass_delta"),
                    "score": row.get("score"),
                    "sequence_tag_length": row.get(
                        "sequence_tag_length"
                    ),
                    "complementary_ion_pair_count": row.get(
                        "complementary_ion_pair_count"
                    ),
                    "candidate_score_margin": row.get(
                        "candidate_score_margin"
                    ),
                }
                for row in feature_region_candidates[:3]
            ],
            "selected_precursor_charge": selected_precursor_charge,
            "peptide_likelihood": peptide_likelihood,
            "peptide_likelihood_reason": peptide_likelihood_reason,
            "coverage_by_sample": coverage_by_sample,
            "best_psm": compact_best_psm,
            "candidate_psm": compact_region_candidate,
            "coverage_scan": coverage_scan,
        })
    return evidence_rows, status_counts


def _candidate_counterpart(
    modified: PeptideCandidate,
    candidates: list[PeptideCandidate],
) -> PeptideCandidate | None:
    if any(name == "C-terminal Lys clipping" for _, name, _ in modified.modifications):
        return next((
            candidate for candidate in candidates
            if not candidate.is_decoy
            and not candidate.modifications
            and candidate.chain == modified.chain
            and candidate.start == modified.start
            and candidate.sequence == modified.sequence + "K"
        ), None)
    unmodified = [
        candidate for candidate in candidates
        if not candidate.is_decoy and not candidate.modifications and candidate.base_peptide_id == modified.base_peptide_id
    ]
    if unmodified:
        return unmodified[0]
    return None


def _candidate_from_psm(psm: dict[str, object]) -> PeptideCandidate | None:
    """Rebuild a searchable form candidate from accepted PSM metadata.

    Targeted glycopeptide and sequence-inference PSMs may use candidates that
    are not present in the bounded standard-search catalog passed to the
    quantitation layer.  Their accepted PSM payload still carries all fields
    needed for targeted MS1 extraction.
    """
    if psm.get("is_decoy"):
        return None
    sequence = str(psm.get("sequence") or "").strip()
    chain = str(psm.get("chain") or "").strip()
    candidate_id = str(psm.get("candidate_id") or "").strip()
    if not sequence or not chain or not candidate_id:
        return None
    try:
        start = int(psm.get("start") or 0)
        end = int(psm.get("end") or 0)
        neutral_mass = float(psm.get("neutral_mass") or 0.0)
    except (TypeError, ValueError):
        return None
    if start <= 0 or end < start or neutral_mass <= 0:
        return None
    modifications: list[tuple[int, str, float]] = []
    for item in psm.get("modifications") or []:
        if not isinstance(item, (list, tuple)) or len(item) < 3:
            continue
        try:
            modifications.append((int(item[0]), str(item[1]), float(item[2])))
        except (TypeError, ValueError):
            continue
    return PeptideCandidate(
        candidate_id=candidate_id,
        base_peptide_id=str(psm.get("base_peptide_id") or f"{chain}:{start}-{end}:{sequence}"),
        chain=chain,
        start=start,
        end=end,
        sequence=sequence,
        neutral_mass=neutral_mass,
        modifications=tuple(modifications),
        carbamidomethyl_cys=bool(psm.get("carbamidomethyl_cys", True)),
        proteolysis=str(psm.get("proteolysis") or "fully_tryptic"),
        is_decoy=False,
    )


def _psm_matches_form(psm: dict[str, object], candidate: PeptideCandidate) -> bool:
    if str(psm.get("candidate_id") or "") == candidate.candidate_id:
        return True
    if not any(name == "C-terminal Lys clipping" for _, name, _ in candidate.modifications):
        return False
    try:
        neutral_mass = float(psm.get("neutral_mass") or 0.0)
    except (TypeError, ValueError):
        return False
    return (
        str(psm.get("chain") or "") == candidate.chain
        and int(psm.get("start") or 0) == candidate.start
        and int(psm.get("end") or 0) == candidate.end
        and str(psm.get("sequence") or "") == candidate.sequence
        and abs(neutral_mass - candidate.neutral_mass) <= 0.01
    )


def _aligned_psm_rt(psm: dict[str, object], shifts: dict[str, object]) -> float:
    if psm.get("feature_rt") is not None:
        return float(psm["feature_rt"])
    return float(psm.get("rt") or 0.0) + float(shifts.get(str(psm.get("sample_id") or "")) or 0.0)


def _envelope_peaks(
    spectra: list[dict[str, object]],
    sample_id: str,
    neutral_mass: float,
    charge: int,
    mz_tolerance_ppm: float,
) -> list[dict[str, float]]:
    targets = [
        (neutral_mass + isotope_offset * ISOTOPE_MASS_DIFF) / charge + PROTON
        for isotope_offset in range(3)
    ]
    rt_values: list[float] = []
    intensity_values: list[float] = []
    for scan in spectra:
        mz_values = scan.get("mz") or []
        intensities = scan.get("intensity") or []
        total = 0.0
        for target in targets:
            tolerance = max(target * mz_tolerance_ppm / 1_000_000.0, 0.01)
            left = bisect_left(mz_values, target - tolerance)
            right = bisect_right(mz_values, target + tolerance)
            total += sum(float(value) for value in intensities[left:right])
        rt_values.append(float(scan.get("aligned_rt") or scan.get("rt") or 0.0))
        intensity_values.append(total)
    trace = XICTrace(
        sample_id=sample_id,
        raw_file_id=sample_id,
        target_mz=targets[0],
        mz_tolerance=max(targets[0] * mz_tolerance_ppm / 1_000_000.0, 0.01),
        rt_array=rt_values,
        intensity_array=intensity_values,
        max_intensity=max(intensity_values, default=0.0),
    )
    return [
        {
            "charge": float(charge),
            "target_mz": targets[0],
            "rt_start": feature.rt_start,
            "rt_apex": feature.rt_apex,
            "rt_end": feature.rt_end,
            "area": feature.area,
            "height": feature.height,
            "signal_to_noise": feature.signal_to_noise,
        }
        for feature in detect_xic_peaks(
            trace,
            min_peak_height=500.0,
            min_peak_area=0.0,
            min_snr=3.0,
            max_peak_width=3.0,
        )
    ]


def _consensus_peak_rt(
    peaks: dict[str, dict[int, list[dict[str, float]]]],
    seed_rt: float | None,
    tolerance_min: float,
) -> float | None:
    if seed_rt is not None:
        return seed_rt
    clusters: list[dict[str, object]] = []
    for sample_id, by_charge in peaks.items():
        for charge, rows in by_charge.items():
            for row in rows[:5]:
                rt = float(row["rt_apex"])
                cluster = min(
                    (item for item in clusters if abs(rt - float(item["rt"])) <= tolerance_min),
                    key=lambda item: abs(rt - float(item["rt"])),
                    default=None,
                )
                if cluster is None:
                    cluster = {"rt": rt, "weighted_rt": 0.0, "weight": 0.0, "members": []}
                    clusters.append(cluster)
                weight = math.log1p(float(row["area"]))
                cluster["weighted_rt"] = float(cluster["weighted_rt"]) + rt * weight
                cluster["weight"] = float(cluster["weight"]) + weight
                members = cluster["members"]
                assert isinstance(members, list)
                members.append((sample_id, charge, row))
                cluster["rt"] = float(cluster["weighted_rt"]) / max(float(cluster["weight"]), 1e-12)
    if not clusters:
        return None
    best = max(
        clusters,
        key=lambda item: (
            len({sample for sample, _, __ in item["members"]}),
            len({charge for _, charge, __ in item["members"]}),
            float(item["weight"]),
        ),
    )
    return float(best["rt"])


def _summarize_form(
    neutral_mass: float,
    psm_rows: list[dict[str, object]],
    peaks: dict[str, dict[int, list[dict[str, float]]]],
    sample_ids: list[str],
    normalization: dict[str, float],
    seed_rt: float | None,
    peak_tolerance_min: float,
) -> dict[str, object]:
    consensus_rt = _consensus_peak_rt(peaks, seed_rt, peak_tolerance_min)
    raw_by_charge: dict[str, dict[str, float]] = {}
    normalized_by_charge: dict[str, dict[str, float]] = {}
    rt_by_sample: dict[str, float | None] = {}
    selected_by_sample: dict[str, list[dict[str, float]]] = {}
    for sample_id in sample_ids:
        selected: list[dict[str, float]] = []
        for charge, rows in peaks.get(sample_id, {}).items():
            if consensus_rt is None:
                continue
            nearest = min(rows, key=lambda row: abs(float(row["rt_apex"]) - consensus_rt), default=None)
            if nearest is not None and abs(float(nearest["rt_apex"]) - consensus_rt) <= peak_tolerance_min:
                selected.append(nearest)
        selected_by_sample[sample_id] = selected
        factor = float(normalization.get(sample_id) or 1.0)
        raw_by_charge[sample_id] = {str(int(row["charge"])): float(row["area"]) for row in selected}
        normalized_by_charge[sample_id] = {
            charge: area * factor for charge, area in raw_by_charge[sample_id].items()
        }
        total_area = sum(float(row["area"]) for row in selected)
        rt_by_sample[sample_id] = (
            sum(float(row["rt_apex"]) * float(row["area"]) for row in selected) / total_area
            if total_area > 0 else None
        )
    normalized_area = {
        sample_id: sum(normalized_by_charge[sample_id].values()) for sample_id in sample_ids
    }
    psm_samples = sorted({str(row.get("sample_id") or "") for row in psm_rows if row.get("sample_id")})
    max_charge_count = max((len(rows) for rows in normalized_by_charge.values()), default=0)
    if psm_samples:
        status = "ms2_confirmed"
    elif max_charge_count >= 2:
        status = "ms1_multi_charge_supported"
    elif any(normalized_area.values()):
        status = "ms1_candidate"
    else:
        status = "not_detected"
    best_psm = max(
        psm_rows,
        key=lambda row: (
            float(row.get("score") or 0.0),
            int(row.get("matched_ion_count") or 0),
            float(row.get("fragment_coverage") or 0.0),
        ),
        default=None,
    )
    return {
        "neutral_mass": neutral_mass,
        "consensus_rt": consensus_rt,
        "status": status,
        "ms2_confirmed_samples": psm_samples,
        "normalized_area_by_sample": normalized_area,
        "normalized_area_by_charge_by_sample": normalized_by_charge,
        "raw_area_by_charge_by_sample": raw_by_charge,
        "charge_states_by_sample": {
            sample_id: sorted(int(charge) for charge in normalized_by_charge[sample_id])
            for sample_id in sample_ids
        },
        "rt_apex_by_sample": rt_by_sample,
        "selected_peaks_by_sample": selected_by_sample,
        "best_psm": best_psm,
    }


def _fold_direction(test: float, reference: float, threshold: float = 1.2) -> tuple[float | None, str]:
    if test <= 0 or reference <= 0:
        return None, "unavailable"
    fold = test / reference
    if fold >= threshold:
        return fold, "increased"
    if fold <= 1.0 / threshold:
        return fold, "decreased"
    return fold, "stable"


def _charge_trend(
    form: dict[str, object],
    reference_sample: str,
    test_sample: str,
) -> dict[str, object]:
    by_sample = form.get("normalized_area_by_charge_by_sample") or {}
    reference = dict(by_sample.get(reference_sample) or {})
    test = dict(by_sample.get(test_sample) or {})
    total_fold, total_direction = _fold_direction(
        float(dict(form["normalized_area_by_sample"]).get(test_sample) or 0.0),
        float(dict(form["normalized_area_by_sample"]).get(reference_sample) or 0.0),
    )
    charge_directions: dict[str, str] = {}
    for charge in sorted(set(reference) & set(test), key=int):
        _, direction = _fold_direction(float(test[charge]), float(reference[charge]))
        charge_directions[charge] = direction
    comparable = [direction for direction in charge_directions.values() if direction != "unavailable"]
    agreeing = sum(direction == total_direction for direction in comparable)
    if len(comparable) >= 2:
        status = "consistent" if agreeing / len(comparable) >= 0.75 else "mixed"
    elif len(comparable) == 1:
        status = "single_charge"
    else:
        status = "unavailable"
    return {
        "status": status,
        "total_fold_test_over_reference": total_fold,
        "total_direction": total_direction,
        "comparable_charge_count": len(comparable),
        "agreeing_charge_count": agreeing,
        "direction_by_charge": charge_directions,
    }


def build_peptide_form_comparisons(
    payload: dict[str, object],
    psms: list[dict[str, object]],
    annotations: list[dict[str, object]],
    candidates: list[PeptideCandidate],
    spectra_loader: Callable[[str], list[dict[str, object]]],
    mz_tolerance_ppm: float = 20.0,
    max_charge: int = 7,
    rt_cluster_tolerance_min: float = 0.5,
) -> list[dict[str, object]]:
    """Pair differential modified peptides with their unmodified parent form.

    Charge states and the first three isotope peaks are quantified as one peptide
    form so they support one conclusion instead of becoming duplicate Features.
    """
    candidate_by_id = {candidate.candidate_id: candidate for candidate in candidates if not candidate.is_decoy}
    for psm in psms:
        rebuilt = _candidate_from_psm(psm)
        if rebuilt is not None:
            candidate_by_id.setdefault(rebuilt.candidate_id, rebuilt)
    available_candidates = list(candidate_by_id.values())
    shifts = dict((payload.get("alignment") or {}).get("rt_shift_by_sample") or {})
    sample_ids = [str(sample) for sample in payload.get("sample_ids") or []]
    reference_sample = str(payload.get("reference_sample") or (sample_ids[0] if sample_ids else ""))
    normalization = dict((payload.get("feature_area_normalization") or {}).get("factor_by_sample") or {})
    all_feature_ids = {
        str(row.get("feature_group_id") or "")
        for row in payload.get("global_feature_groups") or []
        if row.get("feature_group_id")
    }
    significant_feature_ids = {
        str(row.get("feature_group_id") or "")
        for row in payload.get("global_feature_groups") or []
        if row.get("difference_type") in {"presence_absence", "area_changed", "moderate_difference"}
    }

    annotations_by_candidate: dict[str, list[dict[str, object]]] = {}
    for annotation in annotations:
        annotations_by_candidate.setdefault(str(annotation.get("candidate_id") or ""), []).append(annotation)
    triggered_family_candidate_ids: set[str] = set()
    for psm in psms:
        feature_ids = {
            str(link.get("feature_group_id") or "")
            for link in psm.get("feature_links") or []
        }
        if not (feature_ids & significant_feature_ids):
            continue
        candidate = candidate_by_id.get(str(psm.get("candidate_id") or ""))
        if candidate is None:
            continue
        counterpart = _candidate_counterpart(candidate, available_candidates) if candidate.modifications else None
        triggered_family_candidate_ids.add((counterpart or candidate).candidate_id)
    for annotation in annotations:
        if str(annotation.get("feature_or_candidate_id") or "") not in significant_feature_ids:
            continue
        candidate = candidate_by_id.get(str(annotation.get("candidate_id") or ""))
        if candidate is None:
            continue
        counterpart = _candidate_counterpart(candidate, available_candidates) if candidate.modifications else None
        triggered_family_candidate_ids.add((counterpart or candidate).candidate_id)

    groups: list[dict[str, object]] = []
    for annotation in sorted(annotations, key=lambda row: float(row.get("feature_rt") or 0.0)):
        if str(annotation.get("feature_or_candidate_id") or "") not in significant_feature_ids:
            continue
        if str(annotation.get("modification") or "Unmodified") == "Unmodified":
            continue
        modified = candidate_by_id.get(str(annotation.get("candidate_id") or ""))
        if modified is None or not modified.modifications:
            continue
        unmodified = _candidate_counterpart(modified, available_candidates)
        if unmodified is None:
            continue
        rt = float(annotation.get("feature_rt") or 0.0)
        group = next((
            item for item in groups
            if item["unmodified_candidate"].candidate_id == unmodified.candidate_id
            and abs(item["modified_candidate"].neutral_mass - modified.neutral_mass) <= 1e-6
            and abs(rt - float(item["modified_seed_rt"])) <= rt_cluster_tolerance_min
        ), None)
        if group is None:
            group = {
                "modified_candidate": modified,
                "modified_candidates": [modified],
                "unmodified_candidate": unmodified,
                "modified_seed_rt": rt,
                "annotations": [],
                "supporting_annotations": [],
            }
            groups.append(group)
        elif all(candidate.candidate_id != modified.candidate_id for candidate in group["modified_candidates"]):
            group["modified_candidates"].append(modified)
        group["annotations"].append(annotation)
        group["supporting_annotations"].append(annotation)
        group["modified_seed_rt"] = statistics.median(
            float(row.get("feature_rt") or rt) for row in group["annotations"]
        )

    # Expand every triggered family with all accepted modified forms.  This is
    # essential for events such as heavy-chain C-terminal Lys processing: the
    # retained-K peptide can be differential while the clipped counterpart is
    # a common Feature.  A significant member therefore triggers targeted
    # quantitation of every MS2-supported counterpart in the same family.
    for modified in available_candidates:
        if modified.is_decoy or not modified.modifications:
            continue
        unmodified = _candidate_counterpart(modified, available_candidates)
        if unmodified is None:
            continue
        modified_psms = [psm for psm in psms if _psm_matches_form(psm, modified)]
        direct_modified_annotations = annotations_by_candidate.get(modified.candidate_id, [])
        if not modified_psms and not direct_modified_annotations:
            continue
        unmodified_psms = [psm for psm in psms if _psm_matches_form(psm, unmodified)]
        related_candidate_ids = {
            str(psm.get("candidate_id") or "")
            for psm in modified_psms + unmodified_psms
        } | {modified.candidate_id, unmodified.candidate_id}
        related_annotations = [
            row
            for candidate_id in related_candidate_ids
            for row in annotations_by_candidate.get(candidate_id, [])
        ]
        related_feature_ids = {
            str(link.get("feature_group_id") or "")
            for psm in modified_psms + unmodified_psms
            for link in psm.get("feature_links") or []
        } | {
            str(row.get("feature_or_candidate_id") or "")
            for row in related_annotations
        }
        if (
            unmodified.candidate_id not in triggered_family_candidate_ids
            and not (related_feature_ids & significant_feature_ids)
        ):
            continue
        modified_annotations = [
            row
            for candidate_id in {
                str(psm.get("candidate_id") or "") for psm in modified_psms
            } | {modified.candidate_id}
            for row in annotations_by_candidate.get(candidate_id, [])
        ]
        modified_rts = [_aligned_psm_rt(psm, shifts) for psm in modified_psms]
        annotation_rts = [
            float(row.get("feature_rt") or 0.0)
            for row in modified_annotations
            if float(row.get("feature_rt") or 0.0) > 0
        ]
        seed_rt = statistics.median(modified_rts or annotation_rts) if (modified_rts or annotation_rts) else 0.0
        group = next((
            item for item in groups
            if item["unmodified_candidate"].candidate_id == unmodified.candidate_id
            and abs(item["modified_candidate"].neutral_mass - modified.neutral_mass) <= 1e-6
            and (
                not seed_rt
                or not float(item.get("modified_seed_rt") or 0.0)
                or abs(seed_rt - float(item["modified_seed_rt"])) <= rt_cluster_tolerance_min
            )
        ), None)
        if group is None:
            group = {
                "modified_candidate": modified,
                "modified_candidates": [modified],
                "unmodified_candidate": unmodified,
                "modified_seed_rt": seed_rt,
                "annotations": list(modified_annotations),
                "supporting_annotations": list(related_annotations),
            }
            groups.append(group)
        else:
            if all(candidate.candidate_id != modified.candidate_id for candidate in group["modified_candidates"]):
                group["modified_candidates"].append(modified)
            group["annotations"].extend(
                row for row in modified_annotations if row not in group["annotations"]
            )
            group["supporting_annotations"].extend(
                row for row in related_annotations if row not in group["supporting_annotations"]
            )
            if modified_rts:
                group["modified_seed_rt"] = statistics.median(modified_rts)
    if not groups or not sample_ids:
        return []

    forms: dict[tuple[str, str], dict[str, object]] = {}
    for group in groups:
        for label in ("modified", "unmodified"):
            candidate = group[f"{label}_candidate"]
            key = (str(id(group)), label)
            form_psms = (
                [
                    psm
                    for psm in psms
                    if any(_psm_matches_form(psm, modified_candidate) for modified_candidate in group["modified_candidates"])
                ]
                if label == "modified" else [psm for psm in psms if _psm_matches_form(psm, candidate)]
            )
            observed_charges = {
                int(psm.get("precursor_charge") or 0) for psm in form_psms
                if int(psm.get("precursor_charge") or 0) > 0
            }
            charges = sorted({
                charge + offset
                for charge in observed_charges
                for offset in (-1, 0, 1)
                if 1 <= charge + offset <= max_charge
            }) if observed_charges else list(range(1, max_charge + 1))
            psm_rts = [_aligned_psm_rt(psm, shifts) for psm in form_psms]
            if label == "modified" and not psm_rts:
                seed_rt = float(group["modified_seed_rt"])
            else:
                seed_rt = statistics.median(psm_rts) if psm_rts else None
            forms[key] = {
                "candidate": candidate,
                "psms": form_psms,
                "charges": charges,
                "seed_rt": seed_rt,
                "peaks": {sample_id: {} for sample_id in sample_ids},
            }

    for sample_id in sample_ids:
        spectra = spectra_loader(sample_id)
        cache: dict[tuple[float, int], list[dict[str, float]]] = {}
        for form in forms.values():
            candidate = form["candidate"]
            assert isinstance(candidate, PeptideCandidate)
            for charge in form["charges"]:
                cache_key = (round(candidate.neutral_mass, 6), charge)
                if cache_key not in cache:
                    cache[cache_key] = _envelope_peaks(
                        spectra, sample_id, candidate.neutral_mass, charge, mz_tolerance_ppm
                    )
                form["peaks"][sample_id][charge] = cache[cache_key]
        del spectra

    comparisons: list[dict[str, object]] = []
    for index, group in enumerate(groups, start=1):
        modified = group["modified_candidate"]
        unmodified = group["unmodified_candidate"]
        modified_form = forms[(str(id(group)), "modified")]
        unmodified_form = forms[(str(id(group)), "unmodified")]
        modified_summary = _summarize_form(
            modified.neutral_mass,
            modified_form["psms"],
            modified_form["peaks"],
            sample_ids,
            normalization,
            modified_form["seed_rt"],
            rt_cluster_tolerance_min,
        )
        unmodified_summary = _summarize_form(
            unmodified.neutral_mass,
            unmodified_form["psms"],
            unmodified_form["peaks"],
            sample_ids,
            normalization,
            unmodified_form["seed_rt"],
            rt_cluster_tolerance_min,
        )
        modified_areas = dict(modified_summary["normalized_area_by_sample"])
        unmodified_areas = dict(unmodified_summary["normalized_area_by_sample"])
        ratios = {
            sample_id: modified_areas[sample_id] / unmodified_areas[sample_id]
            if unmodified_areas[sample_id] > 0 else None
            for sample_id in sample_ids
        }
        occupancies = {
            sample_id: modified_areas[sample_id] / (modified_areas[sample_id] + unmodified_areas[sample_id])
            if modified_areas[sample_id] > 0 and unmodified_areas[sample_id] > 0 else None
            for sample_id in sample_ids
        }
        pairwise = []
        for test_sample in sample_ids:
            if test_sample == reference_sample:
                continue
            reference_ratio = ratios.get(reference_sample)
            test_ratio = ratios.get(test_sample)
            ratio_fold, direction = _fold_direction(
                float(test_ratio or 0.0), float(reference_ratio or 0.0)
            )
            pairwise.append({
                "reference_sample": reference_sample,
                "test_sample": test_sample,
                "modified_unmodified_ratio_fold": ratio_fold,
                "modification_direction": direction,
                "occupancy_difference": (
                    float(occupancies[test_sample]) - float(occupancies[reference_sample])
                    if occupancies.get(test_sample) is not None and occupancies.get(reference_sample) is not None
                    else None
                ),
                "modified_charge_trend": _charge_trend(modified_summary, reference_sample, test_sample),
                "unmodified_charge_trend": _charge_trend(unmodified_summary, reference_sample, test_sample),
            })
        annotations_for_group = group["annotations"]
        supporting_annotations = group.get("supporting_annotations") or annotations_for_group
        confidence_order = {"B_high_confidence_inferred": 3, "C_tentative": 2, "D_low_evidence": 1}
        best_annotation = max(
            annotations_for_group or supporting_annotations,
            key=lambda row: confidence_order.get(str(row.get("confidence") or ""), 0),
            default={},
        )
        if not best_annotation:
            modified_psm_rows = list(modified_form.get("psms") or [])
            best_psm = max(
                modified_psm_rows,
                key=lambda row: (
                    confidence_order.get(psm_confidence(row), 0),
                    float(row.get("score") or 0.0),
                ),
                default={},
            )
            best_annotation = {
                "confidence": psm_confidence(best_psm) if best_psm else "D_low_evidence",
                "feature_ranking": max(
                    (
                        float(link.get("feature_ranking") or 0.0)
                        for psm in modified_psm_rows
                        for link in psm.get("feature_links") or []
                    ),
                    default=0.0,
                ),
            }
        modification_alternatives = sorted({
            candidate.modification_text for candidate in group["modified_candidates"]
        })
        modified_candidate_ids = {
            candidate.candidate_id for candidate in group["modified_candidates"]
        }
        modified_linked_feature_ids = ({
            str(row.get("feature_or_candidate_id") or "")
            for row in supporting_annotations
            if row.get("feature_or_candidate_id")
            and str(row.get("candidate_id") or "") in modified_candidate_ids
        } | {
            str(link.get("feature_group_id") or "")
            for psm in modified_form.get("psms") or []
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        }) & all_feature_ids
        unmodified_linked_feature_ids = ({
            str(row.get("feature_or_candidate_id") or "")
            for row in supporting_annotations
            if row.get("feature_or_candidate_id")
            and str(row.get("candidate_id") or "") == unmodified.candidate_id
        } | {
            str(link.get("feature_group_id") or "")
            for psm in unmodified_form.get("psms") or []
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        }) & all_feature_ids
        linked_feature_ids = ({
            str(row.get("feature_or_candidate_id") or "")
            for row in supporting_annotations
            if row.get("feature_or_candidate_id")
        } | {
            str(link.get("feature_group_id") or "")
            for form in (modified_form, unmodified_form)
            for psm in form.get("psms") or []
            for link in psm.get("feature_links") or []
            if link.get("feature_group_id")
        }) & all_feature_ids
        feature_rankings = [
            float(row.get("feature_ranking") or 0.0)
            for row in supporting_annotations
        ] + [
            float(link.get("feature_ranking") or 0.0)
            for form in (modified_form, unmodified_form)
            for psm in form.get("psms") or []
            for link in psm.get("feature_links") or []
        ]
        comparisons.append({
            "rank": index,
            "pair_id": f"PAIR_{index:04d}",
            "base_peptide_id": unmodified.base_peptide_id,
            "chain": unmodified.chain,
            "start": unmodified.start,
            "end": unmodified.end,
            "sequence": unmodified.sequence,
            "modification": format_modification_alternatives(
                modification_alternatives,
                str(unmodified.sequence or ""),
            ),
            "modification_alternatives": modification_alternatives,
            "modification_mass_delta": modified.neutral_mass - unmodified.neutral_mass,
            "modified_candidate_id": modified.candidate_id,
            "modified_candidate_ids": sorted(modified_candidate_ids),
            "unmodified_candidate_id": unmodified.candidate_id,
            "linked_feature_ids": sorted(linked_feature_ids),
            "modified_linked_feature_ids": sorted(modified_linked_feature_ids),
            "unmodified_linked_feature_ids": sorted(unmodified_linked_feature_ids),
            "modified_ms2_confidence": best_annotation.get("confidence"),
            "feature_ranking": max(feature_rankings, default=0.0),
            "reference_sample": reference_sample,
            "modified_form": modified_summary,
            "unmodified_form": unmodified_summary,
            "modified_unmodified_ratio_by_sample": ratios,
            "apparent_occupancy_by_sample": occupancies,
            "pairwise_comparisons": pairwise,
        })
    comparisons.sort(key=lambda row: float(row.get("feature_ranking") or 0.0), reverse=True)
    for rank, row in enumerate(comparisons, start=1):
        row["rank"] = rank
        row["pair_id"] = f"PAIR_{rank:04d}"
    return comparisons


def build_modification_level_quantitation(
    peptide_form_comparisons: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Aggregate paired peptide forms into site/proteoform-level proportions.

    The result is deliberately labelled as relative MS-response composition.
    B-level localized evidence can enter the formal relative summary; C-level
    evidence remains visible as tentative, while D-level/mass-only hypotheses
    are excluded from the denominator.
    """
    confidence_order = {
        "B_high_confidence_inferred": 3,
        "C_tentative": 2,
        "D_low_evidence": 1,
    }
    families: dict[str, dict[str, object]] = {}
    for comparison in peptide_form_comparisons:
        confidence = str(comparison.get("modified_ms2_confidence") or "")
        if confidence_order.get(confidence, 0) < confidence_order["C_tentative"]:
            continue
        modification = str(comparison.get("modification") or "").strip()
        is_terminal_lys = "C-terminal Lys clipping" in modification
        is_glycan = "glycan" in modification.lower()
        if is_terminal_lys:
            event_type = "c_terminal_lys_processing"
            event_name = "重链 C 端 Lys 剪切/保留"
            family_key = (
                f"terminal-lys:{comparison.get('chain')}:{comparison.get('start')}:"
                f"{comparison.get('unmodified_candidate_id')}"
            )
        else:
            event_type = "glycoform_distribution" if is_glycan else "site_modification_distribution"
            event_name = "糖型相对构成" if is_glycan else "位点/肽段修饰水平"
            family_key = f"peptide:{comparison.get('unmodified_candidate_id')}"
        family = families.setdefault(family_key, {
            "event_type": event_type,
            "event_name": event_name,
            "base_peptide_id": comparison.get("base_peptide_id"),
            "chain": comparison.get("chain"),
            "start": comparison.get("start"),
            "end": comparison.get("end"),
            "sequence": comparison.get("sequence"),
            "reference_sample": comparison.get("reference_sample"),
            "forms": {},
            "linked_feature_ids": set(),
            "feature_ranking": 0.0,
        })
        family["linked_feature_ids"].update(comparison.get("linked_feature_ids") or [])
        family["feature_ranking"] = max(
            float(family.get("feature_ranking") or 0.0),
            float(comparison.get("feature_ranking") or 0.0),
        )

        unmodified_summary = dict(comparison.get("unmodified_form") or {})
        unmodified_label = "C端 Lys 保留" if is_terminal_lys else "未修饰"
        unmodified_id = f"unmodified:{comparison.get('unmodified_candidate_id')}"
        unmodified_form = family["forms"].setdefault(unmodified_id, {
            "form_id": unmodified_id,
            "label": unmodified_label,
            "modification": "Unmodified",
            "is_unmodified": True,
            "ms2_confidence": "B_high_confidence_inferred" if unmodified_summary.get("status") == "ms2_confirmed" else "C_tentative",
            "status": unmodified_summary.get("status"),
            "normalized_area_by_sample": {},
            "charge_states_by_sample": unmodified_summary.get("charge_states_by_sample") or {},
            "consensus_rts": [],
            "linked_feature_ids": set(),
            "neutral_mass": unmodified_summary.get("neutral_mass"),
            "best_psm": unmodified_summary.get("best_psm"),
            "candidate_ids": [comparison.get("unmodified_candidate_id")],
        })
        unmodified_form["linked_feature_ids"].update(
            comparison.get("unmodified_linked_feature_ids") or []
        )
        if float(dict(unmodified_summary.get("best_psm") or {}).get("score") or 0.0) > float(dict(unmodified_form.get("best_psm") or {}).get("score") or 0.0):
            unmodified_form["best_psm"] = unmodified_summary.get("best_psm")
        for sample_id, area in dict(unmodified_summary.get("normalized_area_by_sample") or {}).items():
            # The same unmodified denominator is extracted once per paired form;
            # retain its largest estimate rather than counting it repeatedly.
            unmodified_form["normalized_area_by_sample"][sample_id] = max(
                float(unmodified_form["normalized_area_by_sample"].get(sample_id) or 0.0),
                float(area or 0.0),
            )
        if unmodified_summary.get("consensus_rt") is not None:
            rt = float(unmodified_summary["consensus_rt"])
            if rt not in unmodified_form["consensus_rts"]:
                unmodified_form["consensus_rts"].append(rt)

        modified_summary = dict(comparison.get("modified_form") or {})
        modified_label = "C端 Lys 已剪切" if is_terminal_lys else modification
        modified_id = (
            f"modified:{round(float(comparison.get('modification_mass_delta') or 0.0), 5)}:"
            f"{modified_label}"
        )
        modified_form = family["forms"].setdefault(modified_id, {
            "form_id": modified_id,
            "label": modified_label,
            "modification": modification,
            "is_unmodified": False,
            "ms2_confidence": confidence,
            "status": modified_summary.get("status"),
            "normalized_area_by_sample": {},
            "charge_states_by_sample": modified_summary.get("charge_states_by_sample") or {},
            "consensus_rts": [],
            "linked_feature_ids": set(),
            "neutral_mass": modified_summary.get("neutral_mass"),
            "best_psm": modified_summary.get("best_psm"),
            "candidate_ids": list(comparison.get("modified_candidate_ids") or []),
        })
        modified_form["linked_feature_ids"].update(
            comparison.get("modified_linked_feature_ids") or []
        )
        modified_form["candidate_ids"] = sorted({
            *[str(candidate_id) for candidate_id in modified_form.get("candidate_ids") or [] if candidate_id],
            *[str(candidate_id) for candidate_id in comparison.get("modified_candidate_ids") or [] if candidate_id],
        })
        if float(dict(modified_summary.get("best_psm") or {}).get("score") or 0.0) > float(dict(modified_form.get("best_psm") or {}).get("score") or 0.0):
            modified_form["best_psm"] = modified_summary.get("best_psm")
        if confidence_order.get(confidence, 0) > confidence_order.get(str(modified_form.get("ms2_confidence") or ""), 0):
            modified_form["ms2_confidence"] = confidence
        for sample_id, area in dict(modified_summary.get("normalized_area_by_sample") or {}).items():
            modified_form["normalized_area_by_sample"][sample_id] = (
                float(modified_form["normalized_area_by_sample"].get(sample_id) or 0.0)
                + float(area or 0.0)
            )
        if modified_summary.get("consensus_rt") is not None:
            rt = float(modified_summary["consensus_rt"])
            if rt not in modified_form["consensus_rts"]:
                modified_form["consensus_rts"].append(rt)

    rows: list[dict[str, object]] = []
    for family_key, family in families.items():
        forms = list(dict(family.pop("forms")).values())
        if family["event_type"] == "glycoform_distribution":
            denominator_forms = [form for form in forms if not form["is_unmodified"]]
            denominator_scope = "identified_glycoforms_only"
        else:
            denominator_forms = forms
            denominator_scope = "unmodified_plus_eligible_modified_forms"
        sample_ids = sorted({
            str(sample_id)
            for form in denominator_forms
            for sample_id in dict(form.get("normalized_area_by_sample") or {})
        })
        quantifiable = len(denominator_forms) >= 2 and bool(sample_ids)
        totals = {
            sample_id: sum(
                float(dict(form.get("normalized_area_by_sample") or {}).get(sample_id) or 0.0)
                for form in denominator_forms
            )
            for sample_id in sample_ids
        }
        for form in forms:
            form["linked_feature_ids"] = sorted(form.get("linked_feature_ids") or [])
            areas = dict(form.get("normalized_area_by_sample") or {})
            included = form in denominator_forms
            form["included_in_denominator"] = included
            form["relative_level_by_sample"] = {
                sample_id: (
                    float(areas.get(sample_id) or 0.0) / totals[sample_id]
                    if included and quantifiable and totals[sample_id] > 0 else None
                )
                for sample_id in sample_ids
            }
        reference_sample = str(family.get("reference_sample") or (sample_ids[0] if sample_ids else ""))
        pairwise: list[dict[str, object]] = []
        for test_sample in sample_ids:
            if not reference_sample or test_sample == reference_sample:
                continue
            form_changes = []
            for form in denominator_forms:
                levels = dict(form.get("relative_level_by_sample") or {})
                reference_level = levels.get(reference_sample)
                test_level = levels.get(test_sample)
                if reference_level is None or test_level is None:
                    continue
                delta = float(test_level) - float(reference_level)
                form_changes.append({
                    "form_id": form["form_id"],
                    "label": form["label"],
                    "reference_level": reference_level,
                    "test_level": test_level,
                    "percentage_point_difference": delta,
                    "level_fold_test_over_reference": (
                        float(test_level) / float(reference_level)
                        if float(reference_level) > 0 else None
                    ),
                    "direction": "increased" if delta > 0.005 else ("decreased" if delta < -0.005 else "stable"),
                })
            strongest = max(form_changes, key=lambda item: abs(float(item["percentage_point_difference"])), default=None)
            pairwise.append({
                "reference_sample": reference_sample,
                "test_sample": test_sample,
                "form_changes": form_changes,
                "strongest_change": strongest,
            })
        formal_forms = [
            form for form in denominator_forms
            if str(form.get("ms2_confidence") or "").startswith("B_")
            and form.get("status") == "ms2_confirmed"
        ]
        if not quantifiable:
            quantitation_status = "insufficient_complementary_forms"
        elif len(formal_forms) == len(denominator_forms):
            quantitation_status = "formal_relative_quantitation"
        else:
            quantitation_status = "tentative_relative_quantitation"
        strongest_overall = max(
            (item["strongest_change"] for item in pairwise if item.get("strongest_change")),
            key=lambda item: abs(float(item["percentage_point_difference"])),
            default=None,
        )
        rows.append({
            **family,
            "family_id": family_key,
            "denominator_scope": denominator_scope,
            "quantifiable": quantifiable,
            "quantitation_status": quantitation_status,
            "relative_ms_response_only": True,
            "total_area_by_sample": totals,
            "forms": forms,
            "pairwise_comparisons": pairwise,
            "max_percentage_point_difference": (
                abs(float(strongest_overall["percentage_point_difference"]))
                if strongest_overall else None
            ),
            "high_value_candidate": bool(
                strongest_overall
                and quantitation_status == "formal_relative_quantitation"
                and abs(float(strongest_overall["percentage_point_difference"])) >= 0.05
            ),
            "linked_feature_ids": sorted(family["linked_feature_ids"]),
        })
    rows.sort(key=lambda row: (
        bool(row.get("high_value_candidate")),
        row.get("quantitation_status") == "formal_relative_quantitation",
        float(row.get("max_percentage_point_difference") or 0.0),
        float(row.get("feature_ranking") or 0.0),
    ), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
        row["quantitation_id"] = f"MODLEVEL_{rank:04d}"
    return rows


def build_proteolytic_modification_families(
    modification_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Summarize nested peptide backbones without adding unlike MS responses.

    Peptides sharing one terminus where one sequence contains the other are
    treated as alternative proteolytic forms.  Their raw areas remain separate;
    only within-backbone modification proportions and directions are combined
    into a site-level consensus.
    """
    eligible = [
        row for row in modification_rows
        if row.get("event_type") in {"site_modification_distribution", "glycoform_distribution"}
        and row.get("quantifiable")
        and row.get("chain")
        and row.get("start") is not None
        and row.get("end") is not None
        and row.get("sequence")
    ]
    parents = list(range(len(eligible)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(first: int, second: int) -> None:
        first_root, second_root = find(first), find(second)
        if first_root != second_root:
            parents[second_root] = first_root

    def nested(first: dict[str, object], second: dict[str, object]) -> bool:
        if str(first.get("chain")) != str(second.get("chain")):
            return False
        first_start, first_end = int(first["start"]), int(first["end"])
        second_start, second_end = int(second["start"]), int(second["end"])
        if first_start == second_start and first_end == second_end:
            return False
        contains = (
            first_start <= second_start <= second_end <= first_end
            or second_start <= first_start <= first_end <= second_end
        )
        if not contains or not (first_start == second_start or first_end == second_end):
            return False
        shorter, longer = sorted(
            (str(first.get("sequence") or ""), str(second.get("sequence") or "")),
            key=len,
        )
        return shorter in longer

    for first_index, first in enumerate(eligible):
        for second_index in range(first_index + 1, len(eligible)):
            if nested(first, eligible[second_index]):
                union(first_index, second_index)

    components: dict[int, list[dict[str, object]]] = {}
    for index, row in enumerate(eligible):
        components.setdefault(find(index), []).append(row)

    summaries: list[dict[str, object]] = []
    for members in components.values():
        if len(members) < 2:
            continue
        members.sort(key=lambda row: (int(row["start"]), int(row["end"])))
        reference_samples = {
            str(row.get("reference_sample") or "") for row in members
            if row.get("reference_sample")
        }
        if len(reference_samples) != 1:
            continue
        reference_sample = next(iter(reference_samples))
        sample_ids = sorted({
            str(sample_id)
            for row in members
            for sample_id in dict(row.get("total_area_by_sample") or {})
        })

        event_support: dict[str, list[dict[str, object]]] = {}
        member_summaries: list[dict[str, object]] = []
        for row in members:
            row_start = int(row["start"])
            sequence = str(row.get("sequence") or "")
            totals = dict(row.get("total_area_by_sample") or {})
            member_summaries.append({
                "quantitation_id": row.get("quantitation_id"),
                "rank": row.get("rank"),
                "start": row.get("start"),
                "end": row.get("end"),
                "sequence": sequence,
                "quantitation_status": row.get("quantitation_status"),
                "total_area_by_sample": totals,
            })
            pairwise_by_test = {
                str(pair.get("test_sample") or ""): pair
                for pair in row.get("pairwise_comparisons") or []
            }
            for form in row.get("forms") or []:
                if form.get("is_unmodified") or not form.get("included_in_denominator"):
                    continue
                modification = str(form.get("modification") or "").strip()
                # Ambiguous or compound forms remain in peptide details, but do
                # not create a cross-peptide localized site conclusion.
                if any(marker in modification for marker in (";", "；", "/", "位点未区分", "方案")):
                    continue
                match = MODIFICATION_SITE_PATTERN.match(modification)
                if not match:
                    continue
                local_position = int(match.group("position"))
                if not 1 <= local_position <= len(sequence):
                    continue
                absolute_position = row_start + local_position - 1
                residue = sequence[local_position - 1]
                modification_name = match.group("name")
                site_label = str(absolute_position) if modification_name.endswith(f"-{residue}") else f"{residue}{absolute_position}"
                event_label = f"{modification_name}@{site_label}"
                changes_by_test: dict[str, dict[str, object]] = {}
                for test_sample, pair in pairwise_by_test.items():
                    change = next((
                        item for item in pair.get("form_changes") or []
                        if item.get("form_id") == form.get("form_id")
                    ), None)
                    if change is not None:
                        changes_by_test[test_sample] = dict(change)
                event_support.setdefault(event_label, []).append({
                    "quantitation_id": row.get("quantitation_id"),
                    "sequence": sequence,
                    "form_label": form.get("label"),
                    "ms2_confidence": form.get("ms2_confidence"),
                    "status": form.get("status"),
                    "changes_by_test": changes_by_test,
                })

        site_conclusions: list[dict[str, object]] = []
        for event_label, support in event_support.items():
            unique_members = {str(item.get("quantitation_id") or "") for item in support}
            if len(unique_members) < 2:
                continue
            comparisons = []
            for test_sample in sample_ids:
                if test_sample == reference_sample:
                    continue
                evidence = []
                for item in support:
                    change = dict(dict(item.get("changes_by_test") or {}).get(test_sample) or {})
                    if not change:
                        continue
                    evidence.append({
                        "quantitation_id": item.get("quantitation_id"),
                        "sequence": item.get("sequence"),
                        "reference_level": change.get("reference_level"),
                        "test_level": change.get("test_level"),
                        "percentage_point_difference": change.get("percentage_point_difference"),
                        "direction": change.get("direction"),
                        "evidence_grade": (
                            "B" if str(item.get("ms2_confidence") or "").startswith("B_")
                            and item.get("status") == "ms2_confirmed" else "C"
                        ),
                    })
                informative = [item for item in evidence if item.get("direction") in {"increased", "decreased"}]
                directions = {str(item.get("direction")) for item in informative}
                consistent = len(evidence) >= 2 and len(informative) == len(evidence) and len(directions) == 1
                deltas = [float(item.get("percentage_point_difference") or 0.0) for item in evidence]
                comparisons.append({
                    "reference_sample": reference_sample,
                    "test_sample": test_sample,
                    "evidence": evidence,
                    "support_count": len(evidence),
                    "consistent": consistent,
                    "direction": next(iter(directions)) if consistent else "mixed_or_stable",
                    "min_percentage_point_difference": min(deltas) if deltas else None,
                    "max_percentage_point_difference": max(deltas) if deltas else None,
                    "median_percentage_point_difference": statistics.median(deltas) if deltas else None,
                    "all_b_grade": bool(evidence) and all(item["evidence_grade"] == "B" for item in evidence),
                })
            strongest = max(
                comparisons,
                key=lambda item: abs(float(item.get("median_percentage_point_difference") or 0.0)),
                default=None,
            )
            site_conclusions.append({
                "event_label": event_label,
                "support": support,
                "pairwise_comparisons": comparisons,
                "high_value_consensus": bool(
                    strongest
                    and strongest.get("consistent")
                    and strongest.get("all_b_grade")
                    and abs(float(strongest.get("median_percentage_point_difference") or 0.0)) >= 0.05
                ),
            })

        digestion_comparisons = []
        for test_sample in sample_ids:
            if test_sample == reference_sample:
                continue
            member_folds = []
            for member in member_summaries:
                totals = dict(member.get("total_area_by_sample") or {})
                reference_area = float(totals.get(reference_sample) or 0.0)
                test_area = float(totals.get(test_sample) or 0.0)
                fold = test_area / reference_area if reference_area > 0 and test_area > 0 else None
                member_folds.append({
                    "quantitation_id": member.get("quantitation_id"),
                    "sequence": member.get("sequence"),
                    "reference_area": reference_area,
                    "test_area": test_area,
                    "test_over_reference_fold": fold,
                })
            valid_folds = [float(item["test_over_reference_fold"]) for item in member_folds if item["test_over_reference_fold"]]
            opposite = any(fold >= 1.2 for fold in valid_folds) and any(fold <= 1 / 1.2 for fold in valid_folds)
            combined_reference = sum(float(item["reference_area"]) for item in member_folds)
            combined_test = sum(float(item["test_area"]) for item in member_folds)
            digestion_comparisons.append({
                "reference_sample": reference_sample,
                "test_sample": test_sample,
                "member_folds": member_folds,
                "opposite_member_directions": opposite,
                "suspected_digestion_redistribution": opposite,
                "uncorrected_combined_area_by_sample": {
                    reference_sample: combined_reference,
                    test_sample: combined_test,
                },
                "uncorrected_combined_test_over_reference_fold": (
                    combined_test / combined_reference if combined_reference > 0 and combined_test > 0 else None
                ),
            })

        if not site_conclusions and not any(
            item.get("suspected_digestion_redistribution") for item in digestion_comparisons
        ):
            continue
        chain = str(members[0].get("chain") or "")
        summaries.append({
            "family_id": f"proteolytic:{chain}:{min(int(row['start']) for row in members)}-{max(int(row['end']) for row in members)}",
            "chain": chain,
            "start": min(int(row["start"]) for row in members),
            "end": max(int(row["end"]) for row in members),
            "reference_sample": reference_sample,
            "members": member_summaries,
            "site_conclusions": site_conclusions,
            "digestion_comparisons": digestion_comparisons,
            "high_value_consensus": any(item.get("high_value_consensus") for item in site_conclusions),
            "relative_ms_response_only": True,
        })
    summaries.sort(key=lambda item: (
        bool(item.get("high_value_consensus")),
        len(item.get("members") or []),
    ), reverse=True)
    for rank, summary in enumerate(summaries, start=1):
        summary["rank"] = rank
    return summaries
