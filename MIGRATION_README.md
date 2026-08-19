# LC-MS 部门迁移包

这是可复制到部门其他 Windows 电脑的代码与运行包。

## 启动

1. 安装 Python 3.10 或更高版本。
2. 迁移包已包含 ThermoRawFileParser 的完整自包含运行目录，可直接处理 Thermo `.raw`；如需替换为其他授权版本，可在启动脚本中传入 `-ParserPath`。仅使用 `.mzML` 时不需要转换器。
3. 双击 `start_department_platform.cmd`（中文名称入口也一并提供）。
4. 浏览器打开 `http://127.0.0.1:8770/`。

## 验证

运行 `powershell -ExecutionPolicy Bypass -File .\tools\verify_migration.ps1`。

## 数据与结果

迁移包刻意不包含 RAW/mzML、任务 jobs、SQLite 结果和运行日志。到新电脑后，从页面上传数据；任务结果会写入 `lcms_department_platform\jobs\<task_id>`。
该包不包含权限管理，适合部门内网或单机使用。
