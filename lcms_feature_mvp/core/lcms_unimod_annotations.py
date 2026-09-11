"""Selective display-only explanations; never used for candidate scoring.

Definitions paraphrase the bundled Unimod full_name/specificity/misc_notes.
Explicit IDs avoid interpreting an arbitrary candidate name as known chemistry.
"""
import re
from functools import lru_cache

from .lcms_unimod import load_catalog
from .lcms_sample_prep import REAGENT_IDS

# Routine modifications intentionally have no hover explanation.
COMMON_NAMES = {"Oxidation", "Dioxidation", "Deamidated", "Deamidation", "Methyl",
                "Dimethyl", "Trimethyl", "Acetyl", "Phospho", "Carbamidomethyl",
                "Dehydrated", "Glycation", "Gln->pyro-Glu", "Glu->pyro-Glu"}
EXPLANATIONS = {
    24: "丙烯酰胺加合物；属于试剂相关衍生物，需核对制样时是否接触丙烯酰胺。",
    43: "N-乙酰己糖胺（HexNAc）糖残基；该名称本身不能区分具体糖异构体。",
    53: "4-羟基壬烯醛（HNE）加合物，属于脂质相关修饰；不是普通的加氧修饰。",
    121: "二甘氨酸（GG）残基标签；数据库用于描述泛素化相关残留，命中不能单独证明发生了泛素化。",
    172: "NBS-12C 化学标记条目（Shimadzu）；不是肽段序列，需核对对应标记试剂。",
    178: "DAET 为 2-(二甲氨基)乙硫醇；条目描述磷酸基团经 β-消除后再进行硫醇衍生化，不是普通磷酸化。",
    187: "精氨酸与羟基苯基乙二醛相关的衍生物；需核对是否进行了对应化学处理。",
    200: "乙二硫醇（EDT）相关衍生化条目；需结合实验试剂判断。",
    272: "肽段 N 端磺化的化学衍生化条目（CAF），不是氨基酸序列。",
    314: "N-甲基马来酰亚胺相关加合物；需核对对应试剂，不能与 N-乙基马来酰亚胺混为一谈。",
    335: "还原型 HNE（4-羟基壬烯醛）相关加合物；名称中的 H(2) 表示相对质量组成变化。",
    378: "羧乙基化相关条目；描述所增加的化学基团，不是新增氨基酸。",
    428: "N-乙酰葡糖胺-1-磷酰基相关修饰，不等同于单独的磷酸化或 HexNAc 加成。",
    447: "去氧相关条目；不能把它等同于制样选项中的二硫键还原。",
    457: "萘-2,3-二甲醛（NDA）相关衍生化；需核对是否使用该化学试剂。",
    528: "数据库描述为脱酰胺后发生甲基化的组合条目；并非本程序分别确认了两个反应。",
    743: "脱水型 4-氧代壬烯醛（4-ONE）Michael 加合物；Delta:H(-2)O(-1) 表示扣除一分子水的组成。",
    830: "甲基乙二醛相关的二羟基加合物条目，不是普通的双氧化。",
    837: "精氨酸转变为硝基嘧啶基鸟氨酸的条目；不是普通氨基酸序列替换的确认。",
    876: "SUMO-1 修饰相关的 QEQTGG 残基标签；不能仅凭质量匹配确认 SUMO 化。",
    931: "数据库描述为脱酰胺后与乙醇发生酯化；需结合样品处理条件解释。",
    949: "3-脱氧葡糖酮相关缩合产物；属于特定反应产物解释，不能仅凭命中认定反应发生。",
    967: "膜蛋白提取相关的化学衍生化条目（thioacylPA）；数据库描述较简略，需要核对具体实验试剂。",
    1041: "脱氧 hypusine 相关修饰：数据库包含 K 位点的氨基丁基转移，以及 Q 位点腐胺相关衍生化规则；具体含义需结合位点。",
    1042: "乙酰化的脱氧 hypusine 相关条目，数据库描述与 eIF5A 的特殊赖氨酸修饰有关。",
    1290: "双 Carbamidomethyl 化相关条目；不同于每个 Cys 一次固定烷基化的默认模型。",
    1301: "转肽反应导致赖氨酸加成的修饰条目；这里的 Lys 不是简单标记序列中存在 K。",
    1355: "形成五元芳香杂环的修饰条目（azole）；该名称不代表已经确定具体环结构。",
    1913: "乙二醛来源的晚期糖基化终产物（AGE），不是酶促糖链修饰。",
    1922: "脯氨酸氧化为 5-羟基-2-氨基戊酸的条目；HAVA 是该产物的缩写。",
    1989: "O-甲基磷硫基相关加合物；数据库描述其来自 O-二甲基磷硫基加合物的脱烷基过程，需核对化学处理背景。",
    1992: "谷氨酰基与血清素相关的修饰（serotonylation），不是普通氨基酸替换。",
    2070: "谷氨酰胺加成的修饰条目；L-Gln 在这里不是简单标记序列中存在 Q。",
    2119: "UFMylation 相关的缬氨酸-甘氨酸（VG）残基标签；数据库描述胰酶消化后留在赖氨酸上的 VG，命中不等于确认 UFMylation。",
}


def modification_explanation(record):
    name, record_id = record["title"], record["id"]
    if name in COMMON_NAMES:
        return ""
    if record_id in EXPLANATIONS:
        return EXPLANATIONS[record_id]
    if record_id in REAGENT_IDS["biotin"] and record_id != 3:
        return "生物素／脱硫生物素探针相关衍生物；需核对具体试剂，不等同于天然生物素化。"
    if record_id in REAGENT_IDS["tmt"]:
        return "TMT 系列化学标记相关条目；需核对具体试剂和版本，不表示本程序已完成 TMT 报告离子定量。"
    if record_id in REAGENT_IDS["itraq"]:
        return "iTRAQ 同位素标记条目；需核对标记试剂，不表示本程序已完成 iTRAQ 专用分析。"
    if name.startswith("Cation:"):
        return "阳离子加合／置换质子的质量模型；不是普通共价修饰，[II]/[III] 等表示离子价态。"
    if re.fullmatch(r"(?:(?:dHex|HexNAc|HexA|Hex|NeuAc|NeuGc|Pent|Kdn|Sulf|Me)\(\d+\))+", name):
        return ("糖组成缩写，括号内为对应单元数。Hex=己糖；HexNAc=N-乙酰己糖胺；"
                "dHex=脱氧己糖；HexA=己糖醛酸；Pent=戊糖；NeuAc/NeuGc/Kdn 为不同唾液酸单元；"
                "Sulf/Me 表示硫酸化/甲基化。组成不等于糖链连接结构，也不能单凭质量确认糖型。")
    return ""


@lru_cache(maxsize=1)
def unimod_name_notes():
    notes = {}
    for record in load_catalog().records:
        explanation = modification_explanation(record)
        if explanation:
            notes[str(record["id"])] = (f"{explanation}\nUnimod 原文：{record['full_name']}"
                                       f"\nUNIMOD:{record['id']} · 数据库定义，仅供理解名称；不代表本样品已确认该修饰或位点。")
    return notes
