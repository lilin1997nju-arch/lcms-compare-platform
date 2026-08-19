# Feature 交接示例

从现有 Peak-first SQLite 导出低相似性 Feature：

```powershell
python .\export_differential_features.py `
  --sqlite D:\LCMS\result\lcms_peak_first_compare.sqlite `
  --output D:\LCMS\result\feature_handoff.json `
  --max-similarity 0.85 `
  --limit 100
```

输出是 `FeatureHandoffRecord/v1` 数组，可作为二级色谱任务的输入。默认读取已经跨 TIC peak 去重的 `global_feature_groups`。

这个示例只负责稳定交接，不修改原 SQLite，也不重新计算 Feature。
