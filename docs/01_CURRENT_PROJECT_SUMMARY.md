# 当前 LC-MS 项目总结

## 1. 项目目标与当前状态

项目已把早期“选定 RT 区间后汇总 200-2000 m/z 总谱”的方式升级为 Peak-first RT-m/z Feature 比较。当前主工作流适用于 2 或 3 个同类型样品，例如 2 个 PTM 样品单独比较，或 3 个 SVA 样品单独比较。

当前系统已经能够：

- 读取 mzML 和 centroid CSV，保留 scan 级 RT、m/z、intensity。
- 计算 TIC、BPC、base peak m/z 和 scan 范围。
- 使用多个色谱特征峰估计全局 RT 平移，并支持用户选择 reference sample。
- 从 consensus TIC 和各样品 TIC 的并集中识别峰，避免只识别共同峰而漏掉差异 TIC 峰。
- 对每个 TIC 峰计算样品特异的局部 RT 漂移和积分范围。
- 对每个独立 TIC 峰比较峰形、面积、峰顶、峰宽、summed spectrum 和 apex spectrum。
- 从可疑 m/z 提取 XIC，形成每个样品的 LCMSFeature。
- 采用 MZmine Join aligner 思想，通过 m/z tolerance、RT tolerance 和权重进行跨样品 Feature 匹配。
- 对缺失样品进行 XIC gap filling，区分真实缺失与积分范围漂移。
- 对相近 RT-m/z 且来自相邻 TIC 峰的重复 Feature 进行全局去重。
- 用低相似性和丰度联合排序差异 Feature。
- 提供 TIC/BPC、summed spectrum、XIC、Feature 表和 Feature-level RT-m/z similarity heatmap。
- 保存确认后的 LCMSFeature，导出交互报告、静态离线报告包和 Feature 数据。
- 通过部门任务平台上传 RAW/mzML、排队转换、计算和查看报告。

## 2. 目录职责

### `lcms_feature_mvp/core`

- `lcms_models.py`：原始文件、scan、候选 m/z、XIC、Feature 和 FeatureGroup 数据类。
- `lcms_parser.py`：mzML 二进制数组解码、centroid CSV 读取、mock 数据生成、目录批量载入。
- `lcms_workbench.py`：TIC/BPC、全局峰 landmark 对齐、全局 RT-m/z 网格热图和旧版 Xcalibur workbench 数据准备。
- `lcms_peak_first.py`：当前核心算法，包括 TIC 峰识别、逐峰局部 RT 校正、谱图比较、XIC Feature、MZmine-style 对齐、gap filling、全局 Feature 去重和排序。
- `lcms_xic.py`：候选 m/z 筛选和基础 XIC 提取。
- `lcms_peak_detection.py`：轻量 XIC 峰识别和梯形积分。
- `lcms_feature_matching.py`：基础 RT+m/z Feature 匹配。
- `lcms_difference.py`：基础 FeatureGroup 差异分类。
- `lcms_export.py`：CSV/JSON/Feature Matrix 导出。

### `lcms_feature_mvp` 顶层

- `run_peak_first_compare.py`：当前主分析入口，同时包含 Peak-first HTML/Canvas/JavaScript 前端模板。
- `serve_peak_first_compare.py`：SQLite-backed 报告 API，支持多个 comparison 共用一个端口。
- `run_xcalibur_workbench.py`、`serve_xcalibur_workbench.py`：较早的全局 RT-m/z workbench，保留用于回溯和算法对照。
- `run_lcms_mvp.py`：最早的 scan/XIC/Feature Matrix MVP。
- `convert_raw_to_mzml.ps1`：Thermo RAW 转 mzML 包装脚本。
- `tests/test_lcms_mvp.py`：解析、XIC、Peak、Feature matching 和差异逻辑测试。

### `lcms_realdata_platform`

- `convert_raw.ps1`：真实数据转换。
- `build_peak_first_compare.ps1`：按 PTM、SVA 或显式样品 ID 建立 2/3 样品 comparison。
- `serve_peak_first_compare.ps1`：启动一个或多个 comparison 的报告服务。
- `build_workbench.ps1`、`serve.ps1`：旧版 workbench 包装脚本。
- `outputs*/lcms_peak_first_compare.html`：当前轻量前端示例，不含 SQLite 数据。

### `lcms_department_platform`

- `server.py`：上传表单、任务状态、单 worker 队列、RAW 转换、Peak-first 调用、报告反向路由及三种报告导出。
- `start_lcms_department_platform.ps1`：默认在 `127.0.0.1:8770` 启动。

## 3. 技术架构

当前分析与报告只依赖 Python 标准库：

- XML ElementTree 解析 mzML。
- `zlib/base64/struct` 解码 mzML binary array。
- `dataclasses/statistics/math` 实现算法。
- `sqlite3` 保存 bootstrap JSON、每个样品的 scan spectra 和人工确认 Feature。
- `http.server` 提供轻量本地 API。
- 原生 HTML/CSS/Canvas/JavaScript 绘图和交互，无 Node 构建步骤。

Thermo `.raw` 转换不是 Python 内置能力，依赖外部 `ThermoRawFileParser` 或等价转换器。

## 4. 当前输出

每个 comparison 的核心结果是：

```text
output/
  lcms_peak_first_compare.html
  lcms_peak_first_compare.sqlite
```

SQLite 的 `peak_first_artifacts` 表保存：

- `bootstrap`：样品、参数、色谱、TIC peaks、逐峰结果、Feature groups、总体 Feature 排名。
- `spectra:<sample_id>`：该样品用于交互谱图和 XIC 的 scan 级数据。

`saved_lcms_features` 表保存用户从页面确认的 Feature JSON。

## 5. 已知限制

- mzML 解析器针对当前 centroid MS1 工作流，尚未系统处理所有 vendor-specific mzML 变体。
- 当前主流程重点是 MS1 Feature，不包含 MS/MS 注释、同位素/电荷去卷积、adduct 注释和蛋白/肽段鉴定。
- 全局 RT 校正是 robust landmark shift，逐 TIC peak 再做局部平移；尚未做非线性 warping。
- Feature similarity 主要基于面积向量 CV，2 个样品时属于工程评分，不是统计显著性。
- 部门任务平台是单机单 worker MVP，未实现用户权限、断点续传、任务取消、分布式队列和审计。
- 前端模板目前嵌在 Python 字符串中。新项目若持续扩展，建议拆成独立静态资源和版本化 API。

## 6. 对下一项目最有价值的稳定资产

建议直接复用：

- `LCMSSpectrumScan` 和 mzML 解析。
- `PeakFirstParams` 集中参数对象。
- TIC peak union detection 与逐峰 local RT correction。
- summed spectrum、XIC、LCMSFeature 和 gap filling。
- FeatureGroup 的 `area_by_sample`、`features_by_sample`、`rt_correction_by_sample`。
- SQLite artifact + API 的轻量报告模式。

建议在新项目中重构后再扩展：

- 将内嵌 HTML 拆分为 `frontend/` 静态资源。
- 将大 JSON artifact 拆为明确的关系表或 Parquet/HDF5 对象存储。
- 将任务平台 worker 替换为可恢复的任务队列。
- 将差异 Feature 到二级色谱任务的交接定义为独立、版本化的数据对象。
