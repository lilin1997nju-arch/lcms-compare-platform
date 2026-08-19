# 运行、部署与报告

## 1. 本地算法运行

将 mzML 放入独立输入目录后：

```powershell
python .\lcms_feature_mvp\run_peak_first_compare.py `
  --input-dir D:\LCMS\mzML `
  --output-dir D:\LCMS\result `
  --include-sample sample_a `
  --include-sample sample_b `
  --reference-sample sample_a
```

启动报告：

```powershell
python .\lcms_feature_mvp\serve_peak_first_compare.py `
  --output-dir D:\LCMS\result `
  --port 8768
```

## 2. 多 comparison 共用接口

```powershell
python .\lcms_feature_mvp\serve_peak_first_compare.py `
  --comparison "ptm|PTM comparison|D:\LCMS\ptm_result" `
  --comparison "sva|SVA comparison|D:\LCMS\sva_result" `
  --port 8768
```

同种比较方式可放在一个服务中通过 comparison 下拉列表切换，不需要为每组数据额外占用端口。

## 3. RAW 转换

迁移包已包含 ThermoRawFileParser 的完整自包含运行目录，上传 Thermo RAW 后由包装脚本自动生成 mzML；如需使用其他授权版本，可通过 `-ParserPath` 指定。转换器版本和命令行参数应写入任务日志，正式服务器还应保存 RAW 与 mzML 的 SHA-256。

## 4. 部门任务平台

三个源码目录必须同级，然后从项目根目录运行：

```powershell
.\lcms_department_platform\start_lcms_department_platform.ps1 -Port 8770
```

当前流程：

```text
上传 RAW/mzML
-> jobs/<task_id>/raw
-> RAW 转 mzML 或复制 mzML
-> run_peak_first_compare.py
-> jobs/<task_id>/output
-> 报告和导出
```

前端入口分成四块：

1. **数据上传**：任务名称、参考样品、多个 RAW/mzML 和核心比较参数；上传前显示文件名和大小，服务器按文件流写入 `jobs/<task_id>/raw`。
2. **运行序列排队**：按创建时间排队，显示排队位置、运行阶段、进度、耗时、上传体积和错误信息；等待中的任务可取消，失败任务可清理生成物后重试。
3. **结果查看**：已完成任务可直接打开在线报告，也可导出单文件 HTML、交互报告包或静态离线包。
4. **任务日志**：点击任务行查看转换器和分析命令的运行日志，便于判断是 RAW 转换失败、输入文件问题还是算法输出失败。

默认是单 worker 顺序执行。迁移到部门电脑时只需安装 Python 3.10+，迁移包已自带 RAW 转换器；如因授权或版本管理不随包部署，可在启动时通过 `-ParserPath` 指向其他转换器。只分析 `.mzML` 时不需要转换器。可双击根目录的 `一键启动部门平台.cmd`，停止时使用 `停止部门平台.cmd`。

当前版本不加入权限管理，适合单机或可信内网先行使用。若未来任务量明显增加，再考虑使用 Windows Service 或独立队列拆分 API、RAW 转换 worker 和分析 worker。

## 5. 报告类型

- 服务型报告：HTML + SQLite + 本地查看服务，保留任意 m/z XIC 互动。
- 单文件报告：压缩数据嵌入 HTML，方便携带但 HTML 较大。
- 静态离线包：轻量 `index.html` 调用同 zip 中的 `data/lcms_offline_data.js`，解压后双击查看，不需要端口。

静态离线包仍包含交互所需 scan spectra，因此体积主要取决于数据，而不是 HTML。若要进一步减小，可增加 lite 模式，只保留已计算 Feature 和预计算 XIC，但会失去任意 m/z 的即时 XIC 能力。

## 6. 服务器化建议

部门推广时建议拆为：

```text
Web/API service
Task queue
RAW conversion worker
LC-MS analysis worker
Result database
Object/file storage
Read-only report service
```

即使暂不加入权限管理，也建议先落实：文件类型/大小限制、任务日志、结果备份和数据保留策略。后续需要跨用户使用时再增加登录和项目权限。

三级结构模块默认把RCSB结构下载代理到后端，并缓存到 `data/structures/pdb_cache`。部门服务器需要允许到 `https://files.rcsb.org` 的出站HTTPS；若服务器不能联网，用户仍可在报告底部上传本地PDB/mmCIF文件完成结构映射。生产部署时建议：

- 对PDB ID接口保留严格格式校验，禁止任意URL代理。
- 对上传结构限制文件类型与大小，当前浏览器端上限为40 MB。
- 定期清理或归档PDB缓存，并记录结构ID、下载时间和文件哈希。
- MS2结果中应保存完整HC/LC FASTA序列，而不只保存链长度，确保报告迁移后仍能绘制全序列差异定位。
- 当前维得利珠单抗自动加载3V4P、地舒单抗自动加载5I1C；两者都是同源模板而非药物精确结构，页面必须保留模板提示和序列一致性说明。
- 抗体复合物结构先以完整HC/LC序列锁定抗体链，再在该链内进行短肽定位，避免把短肽误配到抗原或受体链。
- 把“结构定位”保持为证据展示层，不得自动提升MS2鉴定等级。

## 7. Git 与数据管理

Git 只保存：

- Python/PowerShell/HTML/JS/CSS 源码。
- 小型 mock 数据或 `.gitkeep`。
- 测试和文档。

不要提交：

- RAW、mzML、SQLite、CSV/JSON 大结果。
- `jobs/`、日志、缓存、临时导出包。
- 转换器二进制和授权相关文件。

新项目建议使用独立仓库，并在首个提交中记录本交接包版本和源 commit。
