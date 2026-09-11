"""Explicit, offline cleavage rules used by both task creation and MS2 search.

Rules follow https://www.matrixscience.com/help/enzyme_help.html.
Glu-C variants are exposed separately because buffer-dependent specificity
must be selected from the experiment, not guessed from the enzyme name.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class EnzymeRule:
    label: str
    residues: str
    terminus: str
    blocked_by: str = ""


ENZYMES = {
    "trypsin": EnzymeRule("Trypsin（K/R 后切，后接 P 不切）", "KR", "C", "P"),
    "trypsin_p": EnzymeRule("Trypsin/P（K/R 后切，允许后接 P）", "KR", "C"),
    "lys_c": EnzymeRule("Lys-C（K 后切，后接 P 不切）", "K", "C", "P"),
    "lys_c_p": EnzymeRule("Lys-C/P（K 后切，允许后接 P）", "K", "C"),
    "arg_c": EnzymeRule("Arg-C（R 后切，后接 P 不切）", "R", "C", "P"),
    "glu_c": EnzymeRule("Glu-C / V8-E（E 后切，后接 P 不切）", "E", "C", "P"),
    "glu_c_de": EnzymeRule("Glu-C / V8-DE（D/E 后切，后接 P 不切）", "DE", "C", "P"),
    "asp_n": EnzymeRule("Asp-N（D 前切）", "D", "N"),
    "lys_n": EnzymeRule("Lys-N（K 前切）", "K", "N"),
    "chymotrypsin": EnzymeRule("Chymotrypsin（F/Y/W/L 后切，后接 P 不切）", "FYWL", "C", "P"),
}
DEFAULT_ENZYME = "trypsin"


def validate_enzyme(value: str | None) -> str:
    enzyme = str(value or DEFAULT_ENZYME).strip()
    if enzyme not in ENZYMES:
        raise ValueError(f"Unsupported enzyme / 不支持的酶类型: {enzyme}")
    return enzyme


def digest_enzyme(chains, max_missed_cleavages=2, min_length=6,
                  max_length=60, enzyme=DEFAULT_ENZYME):
    rule = ENZYMES[validate_enzyme(enzyme)]
    peptides = []
    for chain, sequence in chains.items():
        cuts = {0, len(sequence)}
        for index, residue in enumerate(sequence):
            if residue not in rule.residues:
                continue
            cut = index + 1 if rule.terminus == "C" else index
            neighbor = sequence[index + 1:index + 2] if rule.terminus == "C" else sequence[max(0, index - 1):index]
            if neighbor and neighbor in rule.blocked_by:
                continue
            cuts.add(cut)
        ordered = sorted(cuts)
        for left in range(len(ordered) - 1):
            for missed in range(max(0, max_missed_cleavages) + 1):
                right = left + missed + 1
                if right >= len(ordered):
                    break
                start, end = ordered[left], ordered[right]
                if min_length <= end - start <= max_length:
                    peptides.append((chain, start + 1, end, sequence[start:end]))
    return peptides
