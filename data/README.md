# 项目数据目录

这里用于存放可复用的小型参考资源和运行数据。

- `sequences/`：经过确认的 HC/LC FASTA 及来源说明，供 MS/MS 序列匹配使用。
- `structures/pdb_cache/`：已缓存的 PDB/mmCIF 结构，可离线用于结构展示。
- `mzML/`、`denosumab_mzML/` 等目录：本机分析数据，不应随代码迁移包复制。

部门电脑使用 `lcms_department_platform` 页面上传 RAW/mzML 后，任务数据会写入 `lcms_department_platform/jobs/`，不要求手动复制到本目录。
