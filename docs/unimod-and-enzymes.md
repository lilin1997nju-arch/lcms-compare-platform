# 离线 Unimod 二阶段搜索与酶切类型

## 使用方式

新建任务时选择实验实际使用的酶切类型。当前支持 Trypsin、Trypsin/P、
Lys-C、Lys-C/P、Arg-C、Glu-C（仅 E 或 D/E 两种规则）、Asp-N、Lys-N 和
Chymotrypsin。界面同时显示切割残基、方向及 P 阻断规则。
未提供酶类型的旧任务仍按 Trypsin 处理；不是根据文件名猜测酶类型。
该选择会保存到任务参数，并传入常规候选和序列推断候选的生成流程。

提供 FASTA 的新任务默认执行离线 Unimod 第二阶段。选择一个 Feature 后，
在 MS/MS 图上方的“证据来源”切换原有证据或某个 Unimod 探索候选。
“显示样本”限制当前样本范围；不同样本的候选不互相冒用谱图。
谱图显示所选候选对应的真实扫描和定向 b/y 匹配。

旧任务结果不会因软件升级自动重新计算。没有第二阶段数据的报告不能仅靠
更新 HTML 获得新增证据，需要重新分析；现有结果应先保留或另建任务。

## 搜索范围与证据边界

- 完整 XML 在本地，包含单同位素质量差、残基与末端规则、分类、审核状态、
  隐藏标记、规则注释及中性丢失字段；运行期间不联网查询数据库。
- 仅处理尚无 B/C 级解释的显著差异 Feature，且每个样本复用第一阶段挑出的
  一个选中前体/同位素包络扫描。不搜索只有共隔离覆盖、无 MS2 或电荷未知的 Feature。
- 使用现有酶切/截短骨架，先匹配前体质量差和位点规则，再以快速 b/y 预筛保留
  每谱最多 12 个骨架，对保留骨架的匹配修饰及所有允许位点进行独立碎片验证。
  报告记录筛选量及被截断数量；这不是全库穷举。
- 一次只新增一个 Unimod 修饰。固定 Cys 烷基化仍保留，固定 CAM-Cys 不叠加
  相对于未修饰 Cys 定义的其他残基修饰；不支持两个新增修饰的任意组合。
- 蛋白 N/C 端规则使用 FASTA 中的真实起止位置，不能把任意酶切肽端当成蛋白端。
- AA 替换、同位素标记、交联/单端交联试剂、虚拟碎片离子、内部残基删除等
  需要不同模型的条目不参与本阶段。化学衍生化和制样伪影候选必须结合实验条件复核。
- 单独记录数据库中性丢失规则，但本阶段只评分经典 b/y，不用中性丢失或糖链诊断
  离子抬高评分。因此易丢失修饰可能漏检；糖基化命中不能代替原有专用糖肽鉴定。
- 默认质量差范围沿用原序列推断范围 -250 至 +2500 Da，前体容差沿用现有补充搜索
  的 20 ppm、碎片 20 ppm。至少 4 个 b/y 匹配且覆盖率至少 10% 才展示，
  每个样本最多展示 3 个候选。这些是探索展示门槛，不是 FDR 或定位置信度阈值。
- `q_value=null`，`exploratory_only=true`。候选以独立字段
  `feature_evidence[].unimod_rescue_candidates` 保存，不写入原有 accepted PSM、
  annotations 或修饰定量分母，也不把 Unknown 改名为正式肽段。
- 最高分位点仅用于绘制模型，所有其他允许位点及其分数保留，始终标注“位点未确认”。
  宽搜索中的偶然匹配仍可能很多，不应把新增候选数量当成新增可靠鉴定数量。

固定参数不变：`top_n_mz=200`、`candidate_max_per_tic=200`、`top_n_peaks=0`。

## 精简制样条件过滤

新建任务可选生物素化／生物素探针、TMT 系列（含 iodoTMT/cysTMT）和 iTRAQ
三项条件。每项为“使用过／未使用／不清楚”，默认不清楚，旧任务省略时亦如此。
这些条件适用于任务内全部样品；混合制样条件应选择不清楚，避免误排其他样品。

仅“未使用”会根据 `core/lcms_sample_prep.py` 的显式 Unimod ID 映射排除试剂规则，
在质量候选进入骨架预筛和 top-3 排序之前生效，不是事后隐藏结果。生物素条目
UNIMOD:3 只排除人工衍生化末端规则，保留天然 K 位点的 Post-translational 规则。
不会按名称子串或整个 Chemical derivative 分类批量删除，也不排除普通氧化、
脱酰胺、未标记 HNE 等无对应试剂依赖的条目。

“使用过”仅表示不应用该家庭的排除条件，不代表所有该类试剂都真实使用过，
不解除既有 Isotopic label 等模型限制，不新增 TMT/iTRAQ 专用鉴定或报告离子定量。
还原／烷基化设置独立于这三项试剂排除开关，见下节。

CLI 为 `--prep-biotin no --prep-tmt no --prep-itraq no`，各项接受
`unknown`/`no`/`yes`。任务参数、MS2 运行元数据和报告保存 `sample_prep`；
Unimod 汇总额外保存过滤版本、排除数据库条目 ID 和位点规则数（不是命中数）。
原有正式鉴定、评分和定量不变，历史结果不自动重算。

## 还原／烷基化模型（1.1.0）

- 还原状态：已充分还原（DTT、TCEP 等合并）、未还原、不清楚。
  充分还原后的线性肽使用相同质量模型，还原剂不作为固定质量附加到肽段。
- IAA（碘乙酰胺）和 CAA（氯乙酰胺）合并为一个选项，按所有 Cys 固定 Carbamidomethyl
  +57.021463735 Da 建模，常规、序列推断、诱饵、糖肽及 Unimod 的前体与碎片一致。
- 未烷基化：固定 Cys 增量为零，且排除 UNIMOD:4 CAM 探索规则，避免再次
  把已明确未使用的烷基化模型引入。不能据此排除全部化学衍生化。
- 未还原：保守地从常规及所有后续搜索候选池排除含 Cys 的序列（包括诱饵），
  仅提供不含 Cys 线性肽的分析，不宣称二硫键连接肽、游离 Cys 肽或二硫键占有率已被分析。
  MS1 Feature 检测不受此限制；未解释的 Feature 仍保留。
- 不清楚／旧任务缺省：保持原来的线性还原肽及固定 CAM-Cys 假设，界面和
  报告明确提示。不同样品的还原／烷基化方案应分别建任务，不用 unknown 代替混合方案。
- 尚不支持其他烷基化试剂、部分还原、部分烷基化或试剂副反应发生率推断。
  选用还原剂本身不会增加其副反应候选或改变评分。

CLI 增加 `--prep-reduction reduced|none|unknown` 和
`--prep-alkylation cam|none|unknown`。旧 `--no-carbamidomethyl-cys` 等价于
`--prep-alkylation none`，与显式 IAA/CAA 同时使用时报错，避免参数互相覆盖。
`parameters.sample_prep_model` 保存固定质量、非还原限制和警告。

依据：[IAA 的 Unimod 质量和位点](https://www.unimod.org/modifications_view.php?editid1=4)、
[还原与烷基化的系统评估](https://pmc.ncbi.nlm.nih.gov/articles/PMC5500753/)。

## 可重复的对照测试

`tools/benchmark_unimod.py --task <已有任务目录> --output <新的测试目录>`
按 baseline → unimod 顺序运行两个完整 MS2 流程。输入为同一批缓存 mzML，
各自使用任务 SQLite 的副本，不覆盖原始任务。基线通过 CLI
`--disable-unimod-rescue` 关闭新增阶段。

输出 `comparison.json` 包括总耗时、分阶段计时、候选计数，以及移除新增
Unimod 字段后整个报告（包含所有原有定量数据）的精确相等检查。
`tools/replay_unimod.py` 对同一基线结果重复三次新增阶段并验证候选完全一致。
原始日志、两份报告和 SQLite 一并留在测试目录。

这些计时不含 RAW 转换和 MS1 计算；顺序测试可能受缓存/系统负载影响。
总耗时差不能直接等同于新增阶段开销，也不能外推为其他数据集的固定百分比。

## 官方规则来源

- [Unimod 下载及隐藏位点标记](https://www.unimod.org/downloads.html)
- [Unimod 字段、质量差与末端定义](https://www.unimod.org/fields.html)
- [酶切规则表](https://www.matrixscience.com/help/enzyme_help.html)

数据库原始文件、来源校验和及许可证位于 `lcms_feature_mvp/core/data`，
随 Electron 安装包的 `core/**/*` 一起分发。仅维护者显式运行
`tools/update_unimod.py` 时才下载新数据库；桌面运行无此联网步骤。
