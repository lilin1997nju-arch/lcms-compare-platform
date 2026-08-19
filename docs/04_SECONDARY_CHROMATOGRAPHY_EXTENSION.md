# 二级色谱分析扩展指南

## 1. 建议的模块边界

下一项目不要把二级色谱逻辑直接塞入 `prepare_peak_first_payload`。推荐在 Peak-first comparison 完成并由分析人员确认差异 Feature 后创建独立任务：

```text
RAW/mzML
  -> Peak-first comparison
  -> reviewed FeatureHandoffRecord
  -> secondary chromatography job
  -> secondary chromatogram import/alignment/integration
  -> feature-to-fraction/component mapping
  -> review and report
```

这样可以独立重跑二级色谱，而不必重新执行 LC-MS 原始数据分析。

## 2. 新项目建议目录

```text
new_project/
  lcms_feature_mvp/              现有算法基线
  secondary_chromatography/
    models.py                    二级色谱 run、trace、peak、fraction、mapping
    parser.py                    二级色谱文件读取
    alignment.py                 RT/体积/馏分位置对齐
    peak_detection.py            二级色谱峰识别与积分
    feature_mapping.py           LCMSFeature 到二级峰/馏分映射
    service.py                   任务和查询服务
    export.py                    报告与表格导出
  frontend/
    lcms-review.html             可继续复用当前 Peak-first 页面
    secondary-analysis.html      二级色谱工作台
  api/
    server.py
  tests/
  data/                          不进入 Git
```

## 3. 建议新增数据对象

### `SecondaryChromatographyRun`

```text
run_id, project_id, sample_id, source_file,
detector_type, x_unit, y_unit,
gradient_method, column_info, created_at
```

### `SecondaryChromatogramTrace`

```text
run_id, channel_id, x_array, intensity_array,
raw_storage_uri, point_count
```

### `SecondaryChromatographyPeak`

```text
peak_id, run_id, rt_start, rt_apex, rt_end,
area, height, snr, fraction_start, fraction_end,
aligned_rt_apex, match_group_id
```

### `FeatureSecondaryMapping`

```text
mapping_id, source_feature_group_id,
secondary_peak_group_id, sample_id,
mapping_score, rt_consistency_score,
abundance_consistency_score,
mapping_status, reviewer, reviewed_at
```

数组继续使用 JSON/BLOB/Parquet/HDF5 或对象存储，不要将每个点拆成 SQL 行。

## 4. 第一阶段 MVP

1. 从现有 Peak-first SQLite 导入低相似性 Feature 候选。
2. 允许分析人员勾选并创建二级分析任务。
3. 上传二级色谱文件并展示多样品 overlay。
4. 识别二级色谱峰，支持框选放大、pan、手动边界修正和积分。
5. 按二级峰/馏分位置关联已确认 LCMSFeature。
6. 输出每个差异 Feature 在二级色谱中的峰、面积、样品和证据链接。
7. 保存人工复核状态并生成可分享报告。

第一版不建议立即做复杂结构鉴定、MS/MS 注释或机器学习分类。

## 5. 映射评分建议

Feature 到二级色谱峰的映射应保留可解释的分量，例如：

```text
mapping_score =
  0.45 * expected_fraction_or_rt_score
+ 0.35 * abundance_trend_score
+ 0.20 * presence_absence_agreement
```

如果二级色谱与 LC-MS 不是同一次运行或 RT 不可直接对应，应使用 fraction id、收集体积、上样顺序和样品标签，不要强行比较绝对 RT。

## 6. 复用当前前端的方式

当前前端的 Canvas 绘图代码已经支持：

- 框选阴影和缩放。
- 横向 pan。
- 多样品 offset。
- hover 与点击联动。
- TIC peak 的样品特异积分阴影。
- XIC 局部/完整范围切换。

可以抽取为通用模块：

```text
plot-core.js       坐标、缩放、pan、hover
trace-overlay.js   多曲线和 offset
peak-bands.js      样品特异积分区间
feature-link.js    表格、谱图、XIC 联动
```

当前实现嵌在 `run_peak_first_compare.py` 的 `PEAK_FIRST_TEMPLATE`。新项目第一次前端改造时应先做无行为变化的静态资源拆分，并用截图和交互回归测试保护现有功能。

## 7. 数据追溯

每个二级色谱结论至少应能回到：

- 原始 comparison/task id。
- 原始 mzML/RAW 文件哈希或受控数据 id。
- Peak-first 算法版本和参数。
- Feature group、parent TIC peak、m/z tolerance、global/local RT shift。
- XIC 原始点与积分边界。
- 二级色谱源文件、峰边界和人工修订记录。

正式系统不要仅保存截图或导出表格作为结论来源。

## 8. 测试优先级

建议优先补充：

- FeatureHandoffRecord schema round-trip。
- 相同 Feature 从相邻 TIC peaks 进入时的去重。
- 局部 RT shift 对 XIC 积分边界的影响。
- 二级峰连续/不连续边界和手动修订。
- 样品、Feature、二级峰关联的不可变追溯 id。
- 大文件上传失败、重复任务和任务恢复。
