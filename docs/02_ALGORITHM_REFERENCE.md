# Peak-first 算法说明

## 1. 数据读取

输入优先为 centroid mzML。每个 MS1 scan 解析为：

```text
scan_id, sample_id, rt, ms_level,
mz_array, intensity_array,
tic, base_peak_mz, base_peak_intensity
```

`mz_array` 与 `intensity_array` 始终成对保存，不拆成逐点 SQL 行。

## 2. 全局 RT 对齐

先在参考样品和其他样品中识别多个高丰度 chromatogram landmarks。候选峰按 RT 邻近进行一对一匹配，得到多组：

```text
rt_shift = reference_peak_rt - sample_peak_rt
```

使用中位数和 MAD 去除异常 shift，得到每个样品的 robust global shift：

```text
global_aligned_rt = raw_rt + global_rt_shift
```

如果 landmark 不足，则回退到主峰/最大面积峰对齐。reference sample 可以由用户指定。

## 3. TIC 峰识别

算法对 TIC 进行短窗口平滑和 rolling-min baseline 估计，再用低强度点的 median/MAD 估计噪声。候选峰需满足 S/N、prominence、面积、面积占比和峰宽阈值。

为避免峰尖噪声导致过窄切峰：

- 峰顶候选使用额外平滑后的 peak body 信号。
- 浅谷 shoulder 合并为一个峰，深谷保留为独立峰。
- 峰边界沿信号下降到噪声或峰高比例阈值处扩展。
- 连续峰在谷底共用边界，尽量不留空白。
- 孤立峰向空白区做小幅边缘延伸，避免起点过晚、终点过早。

候选源同时包括：

- 多样品对齐后的 consensus TIC peaks。
- 每个样品独立识别的 sample TIC peaks。

两者合并后可发现只出现在某个样品中的差异 TIC 峰。

## 4. 逐 TIC 峰局部 RT 校正

每个 TIC peak 独立进行局部校正。reference sample 固定不动，其他样品在受峰宽和绝对上限约束的搜索窗口内，使用基线校正、平滑后的完整 TIC 峰轮廓与 reference 做互相关匹配。评分同时考虑强度轮廓 cosine 和轮廓斜率一致性：

```text
local_shift(reference, peak) = 0
local_shift(sample, peak) = argmax(profile_similarity(reference, shifted_sample))
```

这避免复合峰中因主峰/肩峰相对高度改变而将不同子峰错误配对。若最佳轮廓分数过低，或位移带来的改善不足，自动拒绝该局部位移并保留 0；触及搜索边界的结果标记为 capped。最终用于该峰的 RT 为：

```text
peak_aligned_rt = raw_rt + global_shift(sample) + local_shift(sample, peak)
```

每个样品的实际切峰范围反向移动：

```text
sample_cut_start = aligned_peak_start - local_shift
sample_cut_end   = aligned_peak_end   - local_shift
```

这一步解决不同 TIC 峰漂移程度不同、复合峰最高点不对应，以及某样品 XIC 只落入积分框一小部分时产生的假差异。

## 5. TIC 峰一致性

每个 TIC peak 分别计算：

- shape score：归一化重采样峰形与 consensus shape 的 cosine。
- area score：跨样品面积一致性。
- apex score：局部校正前后峰顶位置一致性。
- width score：峰宽一致性。
- chromatogram score：上述分量的组合。
- spectrum score：summed spectrum 和 apex spectrum 的综合一致性。

总 peak consistency 当前偏重质谱：

```text
0.15 * shape
+ 0.10 * area
+ 0.05 * apex
+ 0.05 * width
+ 0.65 * spectrum
```

TIC peak 状态包括 `high_consistency`、`chromatogram_changed`、`chromatogram_consistent_spectrum_changed`、`possible_rt_shift`、`need_local_rt_mz_analysis` 和 `moderate_consistency`。

## 6. m/z 合并与可疑 m/z

在每个 TIC peak 的样品特异窗口内分别求 summed spectrum 和 apex spectrum。原始谱图中的质心点先按约 `10 ppm` 聚类和跨样品匹配；`mz_tolerance_da` 不再用于把相邻质心峰合成一个平均 m/z。

- 跨样品同峰匹配：当前默认 `10 ppm`。
- 单个质心峰的 XIC 搜索窗口：当前默认 `±0.16 Da`。
- Feature 的 `representative_mz`/`true_peak_mz` 取自实际观测质心峰，不再是整个 Da 窗口内多个峰的加权平均值。

changed m/z 候选仅按 `0.02 Da` 的最小间隔去重，使间距约 `0.5017 Da`（z=2）和 `0.3345 Da`（z=3）的同位素峰能够保留下来，供后续包络识别。

## 7. XIC 与 LCMSFeature

对每个候选 m/z，在 TIC peak 周围扩展上下文窗口提取 XIC。每个 scan 先在 `±0.16 Da` 内寻找最接近候选值的真实质心峰；只有该质心峰本身仍位于目标 m/z 的 `10 ppm` 内才接受，并只汇总该质心峰 `10 ppm` 范围内的信号。若窗口内只有邻近离子而目标离子缺失，该 scan 强度记为 0，不能让邻峰代替目标峰。因此，相邻同位素峰或其他离子既不会被平均，也不会制造错误的 XIC 峰型和 presence/absence。

XIC Feature 的峰顶为最高强度点，边界沿信号下降至噪声或峰高比例阈值处扩展，面积采用梯形积分。检测阈值不足时不直接丢弃，而是作为 `gap_filled` 保存，以支持跨样品追踪。

## 8. MZmine-style Feature alignment

每个 TIC peak 内，先由检测到的 Feature 得到 representative m/z 和 representative aligned RT。m/z 首先应用 `10 ppm` 硬门限：超出门限的候选不能作为同一真实峰匹配；门限内再计算匹配分数：

```text
mz_score = max(0, 1 - abs(mz - representative_mz) / mz_tolerance)
rt_score = max(0, 1 - abs(rt - representative_rt) / rt_tolerance)

match_score =
  (mz_weight * mz_score + rt_weight * rt_score)
  / (mz_weight + rt_weight)
```

当前默认：

- `cross_sample_mz_tolerance_ppm = 10`
- `feature_mz_weight = 0.7`
- `feature_rt_weight = 0.3`
- `feature_rt_tolerance_min = 0.12 min`
- 低于 0.30 且不能由 RT drift 解释的匹配标为 `low_score`

后续组分层会按同位素间距合并包络：

- z=2：理论间距约 `1.00335 / 2 = 0.5017 Da`
- z=3：理论间距约 `1.00335 / 3 = 0.3345 Da`
- 不同电荷态若可由同一中性质量解释，并且 RT、样品间变化趋势和包络证据一致，则合并为同一组分。

FeatureGroup 同时保存 `true_peak_mz`（实际质心峰）、`quantitation_mz`（XIC 定量锚点）和 `envelope_representative_mz`（同位素包络代表峰）。已合并组分还保存 `component_neutral_mass` 和代表电荷态。

## 9. Feature similarity、差异和丰度排序

每个 FeatureGroup 收集所有样品面积，先进行 local TIC normalization 和 `log1p`，再计算：

```text
CV = population_std(normalized_area) / mean(normalized_area)
similarity = 1 / (1 + CV)
difference = 1 - similarity
```

差异类型：

- 所有样品均缺失：`low_confidence`
- 部分样品检测到：`presence_absence`
- similarity < 0.60：`area_changed`
- similarity < 0.80：`moderate_difference`
- 其他：`common_feature`

有限的原始变化倍数保存在 `raw_max_fold_change`。用于报告、排序和可视化的 `max_fold_change` 最高记为 `1000`；原始面积不截断，因此仍可回溯实际定量值。

单 TIC peak 内排序使用：

```text
abundance_score = log1p(max_area)
ranking = difference * (1 + abundance_weight * abundance_score)
```

总体 Feature 表先将 abundance 归一化，再联合 difference 排序，使低丰度噪声不占据前列，同时保留真实低丰度差异的可见性。

## 10. 跨 TIC peak 去重

同一 RT-m/z Feature 可能被相邻 TIC peak 的上下文窗口重复识别。总体表按以下条件聚类：

```text
abs(mz1 - mz2) <= mz_tolerance * global_merge_factor
abs(rt1 - rt2) <= global_feature_rt_merge_tolerance
```

当前默认 RT 合并容差为 `0.8 min`，m/z 合并因子为 `2.0`。聚类保留 ranking 最高的代表 Feature，同时记录全部来源 Feature group 和 parent TIC peak IDs。

## 11. Feature heatmap

热图以去重后的 FeatureGroup/组分为绘制单元，不再绘制所有原始 RT-m/z 网格：

- 无 Feature 信号处保持空白。
- 两样品方向模式使用带符号的 `ln(test/reference)`：供试样品升高为红色、降低为蓝色。
- 倍数幅度模式使用 `ln(fold)` 连续着色。
- 两种模式都在 `1000×`（即 `|ln(fold)| = ln(1000)`）封顶，避免极端 presence/absence 比值拉爆比例尺。
- 点大小同时考虑 `ln(fold)`、abundance 和 ranking；原始面积不因封顶而改变。
- 重叠时差异优先级较高的 Feature 后绘制并位于顶层。
- hover 显示 m/z、RT、TIC peak、rank 等信息。
- 点击 Feature 可联动 TIC peak、summed spectrum、XIC 和表格。

## 12. 参数入口

集中参数定义在 `PeakFirstParams`。CLI 当前开放常用参数，其他参数可在新项目的配置层逐步开放。建议新项目不要在业务代码中散落常量，而是为每次 analysis run 保存完整参数快照和算法版本。
