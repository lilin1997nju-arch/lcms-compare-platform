# LC-MS 部门任务平台

这个目录提供一个轻量级部门内网入口，用于手动上传 RAW/mzML 文件、建立两样品或多样品比对任务、后台自动转换 mzML 并调用现有 Peak-first LC-MS 分析流程。

## 当前功能

- 新建 LC-MS 比对任务
- 上传至少 2 个 `.raw` 或 `.mzML` 文件
- `.raw` 自动调用迁移包内置的 `ThermoRawFileParser` 转换为 `.mzML`
- 后台单 worker 顺序执行任务队列
- 自动调用 `lcms_feature_mvp/run_peak_first_compare.py`
- 任务完成后提供报告入口
- 报告页面复用现有 Peak-first Compare 交互界面
- 导出可携带交互报告包，在其他电脑解压后用本地 Python 查看
- 上传页面显示文件清单、大小和格式校验
- 队列页面显示排队位置、运行阶段、进度、耗时和数据量
- 等待中的任务可以取消，失败任务可以清理结果后重试
- 结果支持在线查看、单文件 HTML、交互报告包和静态离线包

## 启动

从工作区根目录运行：

```powershell
.\lcms_department_platform\start_lcms_department_platform.ps1
```

也可以双击工作区根目录的 `一键启动部门平台.cmd`。它会优先寻找项目自带的虚拟环境或系统 Python，并自动打开浏览器。

默认访问：

```text
http://127.0.0.1:8770/
```

## 数据目录

每个任务会生成独立目录：

```text
lcms_department_platform/jobs/<task_id>/
  raw/      上传的原始文件
  mzML/     转换或复制后的 mzML
  output/   Peak-first HTML + SQLite 结果
  worker.log
```

## 注意

这是一个不包含权限管理的单机/部门内网版本，默认单 worker 顺序执行。当前上传上限为 8 GB，文件会流式写入任务目录，避免把整个 RAW 文件一次性读入内存。

如果后续任务量明显增加，再考虑把 worker 拆成 Windows 服务或队列服务；本次迁移包不引入权限管理，也不要求额外数据库。

建议服务器运维保留：

- `state/tasks.json`：任务状态，可随项目备份
- `jobs/<task_id>/raw`：原始上传数据
- `jobs/<task_id>/output`：HTML 和 SQLite 结果
- `logs/`：平台启动日志

不建议直接删除正在运行任务的目录；如需释放空间，优先从已导出的报告和已确认备份的任务开始归档。

## 导出报告

任务完成后，任务列表和报告页顶部都会出现导出入口。推荐优先使用：

```text
导出静态离线包
```

这个包解压后直接双击 `index.html` 即可查看，不需要启动端口服务。HTML 是轻量前端，交互数据保存在同包内的 `data/lcms_offline_data.js`。

也可以使用旧的服务型报告包：

```text
导出报告包 / 导出可交互报告包
```

服务型 zip 包不包含 RAW/mzML 原始数据，只包含：

- 报告 HTML
- Peak-first SQLite 结果库
- 轻量本地查看服务器
- Windows 启动脚本

在其他电脑解压后运行：

```powershell
powershell -ExecutionPolicy Bypass -File .\start_report_server.ps1
```

然后打开：

```text
http://127.0.0.1:8780/
```
