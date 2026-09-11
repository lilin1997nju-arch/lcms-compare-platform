"""Conservative, versioned reagent exclusions for exploratory Unimod search.

IDs and specificity exceptions are reviewed against the bundled Unimod XML.
Do not infer reagent use from names, classify all chemical derivatives alike,
or turn these family-level switches into support for isobaric quantitation.
"""

PREP_FILTER_VERSION = 2
PREP_FIELDS = {
    "biotin": "生物素化／生物素探针",
    "tmt": "TMT 系列标记（含 iodoTMT、cysTMT）",
    "itraq": "iTRAQ 标记",
}
PREP_STATES = {"unknown": "不清楚／未填写", "no": "未使用", "yes": "使用过"}
REDUCTION_CHOICES = {
    "unknown": "不清楚（沿用线性还原肽模型）",
    "reduced": "已充分还原（DTT、TCEP 等）",
    "none": "未还原（仅搜索不含 Cys 的线性肽）",
}
ALKYLATION_CHOICES = {
    "unknown": "不清楚（沿用固定 CAM-Cys 假设）",
    "cam": "IAA/CAA 类（Cys +57.021464 Da）",
    "none": "未烷基化（Cys 不添加固定质量）",
}
ALL_PREP_FIELDS = {"reduction": "还原状态", "alkylation": "烷基化方式", **PREP_FIELDS}
PREP_CHOICES = {"reduction": REDUCTION_CHOICES, "alkylation": ALKYLATION_CHOICES,
                **{key: PREP_STATES for key in PREP_FIELDS}}

# Explicit IDs, not substring matching. Includes probes with non-obvious names
# (DBIA, DPIA, HNE-BAHAH, AMTzHexNAc2), but not unlabelled HNE itself.
REAGENT_IDS = {
    "biotin": frozenset({
        3, 20, 89, 92, 93, 112, 113, 114, 115, 116, 117, 118, 289, 290,
        294, 325, 332, 333, 343, 353, 357, 361, 522, 523, 538, 539, 774,
        800, 811, 884, 895, 912, 934, 993, 1012, 1031, 1037, 1039,
        1251, 1252, 1314, 1320, 1340, 1423, 1830, 1841, 2052, 2053,
        2062, 2067, 2068, 2069, 2106, 2126, 2127, 2128,
    }),
    "tmt": frozenset({737, 738, 739, 984, 985, 1341, 1342, 2015, 2016,
                      2017, 2050, 2122, 2123}),
    "itraq": frozenset({214, 532, 533, 730, 731}),
}


def normalize_sample_prep(value=None):
    if value is None:
        value = {}
    if not isinstance(value, dict) or set(value) - set(ALL_PREP_FIELDS):
        raise ValueError("制样条件格式无效")
    result = {}
    for key, label in ALL_PREP_FIELDS.items():
        state = value.get(key, "unknown")
        if not isinstance(state, str) or state not in PREP_CHOICES[key]:
            raise ValueError(f"{label}：请选择支持的制样选项")
        result[key] = state
    return result


def sample_prep_model(value=None):
    prep = normalize_sample_prep(value)
    cam = prep["alkylation"] != "none"
    warnings = []
    if prep["reduction"] == "none":
        warnings.append("未还原模式：仅搜索不含 Cys 的线性肽；不鉴定二硫键连接肽、含 Cys 肽或二硫键占有率。")
    elif prep["reduction"] == "unknown":
        warnings.append("还原条件未知：沿用线性还原肽假设，不支持二硫键连接肽鉴定。")
    if prep["alkylation"] == "unknown":
        warnings.append("烷基化条件未知：沿用旧版固定 CAM-Cys 假设（每个 Cys +57.021464 Da），请核对实验条件。")
    elif cam:
        warnings.append("IAA/CAA 按每个 Cys 完全烷基化建模；不同时枚举未反应 Cys 或其他烷基化试剂。")
    return {"carbamidomethyl_cys": cam,
            "fixed_modification": "Carbamidomethyl@C" if cam else "None",
            "fixed_cys_delta_da": 57.021463735 if cam else 0.0,
            "exclude_cys_peptides": prep["reduction"] == "none",
            "disulfide_search_supported": False, "warnings": warnings}


def filter_linear_candidates(candidates, value=None):
    if sample_prep_model(value)["exclude_cys_peptides"]:
        return [candidate for candidate in candidates if "C" not in candidate.sequence]
    return candidates


def sample_prep_cli_args(value=None):
    return [part for key, state in normalize_sample_prep(value).items()
            for part in (f"--prep-{key}", state)]


def filter_reagent_rules(record, sample_prep):
    """Preserve endogenous K biotinylation when artificial labelling is absent."""
    # Explicit absence of alkylation must not be undone by the rescue search.
    # UNIMOD:4 is the CAM adduct itself. Other derivative families are not
    # inferred absent simply from this one condition.
    if sample_prep.get("alkylation") == "none" and record.get("id") == 4:
        return []
    excluded = any(sample_prep[key] == "no" and record.get("id") in ids
                   for key, ids in REAGENT_IDS.items())
    if not excluded:
        return record["rules"]
    if record.get("id") == 3:
        return [rule for rule in record["rules"]
                if rule.get("site") == "K"
                and rule.get("position") == "Anywhere"
                and rule.get("classification") == "Post-translational"]
    return []
