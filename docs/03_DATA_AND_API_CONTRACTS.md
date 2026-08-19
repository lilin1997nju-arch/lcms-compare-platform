# 数据与 API 契约

## 1. 输入文件

### mzML

当前解析器读取 MS1 spectrum、scan start time、m/z array 和 intensity array，支持 32/64 bit float 和 zlib 压缩 binary array。RT 统一转换为分钟。

### centroid CSV

用于 mock 和轻量交换，基本字段为：

```text
sample_id, scan_id, rt, ms_level, mz, intensity
```

同一 `sample_id + scan_id + rt` 的多行聚合为一个 scan。

## 2. 核心 Python 对象

### `LCMSRawFile`

```text
raw_file_id, sample_id, project_id, file_name, file_path,
data_format, mz_min, mz_max, rt_min, rt_max, scan_count, created_at
```

### `LCMSSpectrumScan`

```text
scan_id, raw_file_id, sample_id, rt, ms_level,
mz_array, intensity_array, tic,
base_peak_mz, base_peak_intensity
```

### `LCMSFeature`

```text
feature_id, sample_id, raw_file_id, mz, mz_tolerance,
rt_start, rt_apex, rt_end, aligned_rt_apex,
area, height, signal_to_noise,
feature_group_id, match_status, annotation_status
```

### `LCMSFeatureGroup`

```text
feature_group_id, representative_mz, true_peak_mz, quantitation_mz,
envelope_representative_mz, representative_rt,
sample_count, missing_sample_count, cv, max_fold_change, raw_max_fold_change,
difference_type, confidence_score,
area_by_sample, normalized_area_by_sample, feature_ids
```

其中 `representative_mz` 保留为兼容字段，当前与 `true_peak_mz` 同义。`envelope_representative_mz` 表示同位素包络的代表 m/z，不应当替代真实质心峰用于 XIC 定量。

`max_fold_change` 是报告和热图使用的倍数，最大为 `1000`；`raw_max_fold_change` 保留未截断的有限倍数，原始面积字段始终不截断。

Peak-first payload 使用等价的 dict，并额外提供 `features_by_sample`、`match_score_by_sample`、`rt_correction_by_sample`，以及组分合并后的 `component_neutral_mass`、`component_representative_charge` 和子峰同位素序号。

## 3. SQLite

### `peak_first_artifacts`

```sql
artifact_key TEXT PRIMARY KEY,
payload_json TEXT NOT NULL,
updated_at TEXT NOT NULL
```

已使用的 key：

- `bootstrap`
- `spectra:<sample_id>`

### `saved_lcms_features`

```sql
feature_id TEXT PRIMARY KEY,
feature_json TEXT NOT NULL,
created_at TEXT NOT NULL,
updated_at TEXT NOT NULL
```

`bootstrap` 中最重要的路径：

```text
sample_ids
reference_sample
params
alignment.rt_shift_by_sample
chromatograms.<sample>.raw/aligned
tic_peaks[]
peak_results[]
peak_results[].sample_cut_bounds
peak_results[].local_alignment
peak_results[].feature_groups[]
global_feature_groups[]
```

## 4. 报告 API

### `GET /api/comparisons`

返回当前端口中的 comparison 列表和默认 comparison。

### `GET /api/bootstrap?comparison=<id>`

返回该 comparison 的完整 Peak-first bootstrap payload。

### `GET /api/xic`

参数：

```text
peak_id=<TICP_xxxx>
mz=<numeric>
comparison=<id>
full=0|1
```

返回：

```text
target_mz, mz_tolerance,
rt_start, rt_end,
integration_rt_start, integration_rt_end,
context_rt_start, context_rt_end,
sample_cut_bounds,
feature_bounds_by_sample,
xic_by_sample,
integration_by_sample
```

### `GET /api/features`

返回用户保存的 Feature。

### `POST /api/features`

请求：

```json
{"features": [{"feature_id": "..."}]}
```

当前实现为全量 replace/upsert 语义。

## 5. 建议新增的二级色谱交接对象

新项目建议定义 `FeatureHandoffRecord/v1`：

```json
{
  "schema_version": "FeatureHandoffRecord/v1",
  "source_result_id": "comparison-or-task-id",
  "source_sqlite": "relative/or/managed/object/id",
  "feature_group_id": "TICP_0074_FG_0001",
  "parent_tic_peak_ids": ["TICP_0073", "TICP_0074"],
  "representative_mz": 627.68738,
  "true_peak_mz": 627.68738,
  "quantitation_mz": 627.68738,
  "envelope_representative_mz": 627.68738,
  "component_neutral_mass": 1880.04023,
  "component_representative_charge": 3,
  "representative_aligned_rt": 40.010,
  "similarity_score": 0.685,
  "difference_score": 0.315,
  "ranking_score": 0.493,
  "difference_type": "moderate_difference",
  "max_fold_change": 2.70,
  "sample_measurements": {
    "sample-a": {
      "area": 1117333.2,
      "height": 0.0,
      "rt_start": 39.8,
      "rt_apex": 40.0,
      "rt_end": 40.2,
      "global_rt_shift": 0.0,
      "local_rt_shift": 0.01,
      "match_status": "matched",
      "gap_filled": false
    }
  },
  "review": {
    "status": "candidate",
    "reviewer": null,
    "reviewed_at": null,
    "notes": ""
  }
}
```

`source_sqlite` 在正式服务器上最好替换为受控的 result/object id，避免跨机器使用绝对路径。

## 6. API 版本化建议

新项目建议从 `/api/v1/` 开始，不直接修改旧报告 API：

```text
GET  /api/v1/comparisons/{id}/features
POST /api/v1/secondary-jobs
GET  /api/v1/secondary-jobs/{id}
GET  /api/v1/secondary-jobs/{id}/chromatograms
POST /api/v1/secondary-jobs/{id}/review
GET  /api/v1/secondary-jobs/{id}/report
```

Feature 交接请求应保存算法版本、参数快照和 source result id，以保证追溯。
