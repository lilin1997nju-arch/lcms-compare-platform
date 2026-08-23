# LC-MS 差异组分分析基础交接包

版本日期：2026-08-12  
算法基线：`lcms-compare-platform` main / commit `74ea4a9`  
用途：在新的工作目录中快速复用当前 LC-MS Peak-first 差异分析能力，并继续开发差异组分的二级色谱分析。

## 包含内容

- `lcms_feature_mvp/`：mzML/centroid CSV 解析、TIC/BPC、Peak-first、XIC、LCMSFeature、MZmine-style feature alignment、SQLite 报告接口与 Canvas 前端。
- `lcms_realdata_platform/`：RAW 转换、按 2 或 3 个样品构建比较、启动报告的 PowerShell 包装脚本，以及当前前端 HTML 示例。
- `lcms_department_platform/`：上传、排队、RAW 自动转换、分析任务、交互报告和离线报告包导出的单机部门平台 MVP。
- `examples/secondary_chromatography/`：从 Peak-first SQLite 提取低相似性 Feature 的交接示例。
- `docs/`：项目总结、算法说明、数据/API 契约、二级色谱扩展建议和部署说明。

包内不包含 RAW、mzML、SQLite 结果库、任务缓存和运行日志。原始数据与结果库体积大，且应由新项目自己的数据目录管理。

## 部门电脑快速启动

本项目根目录提供：

- `一键启动部门平台.cmd`：启动上传、排队和结果查看页面，默认端口 `8770`。
- `停止部门平台.cmd`：停止由启动脚本记录的部门平台进程。
- `打包部门迁移包.cmd`：生成不含原始数据和历史结果的代码迁移包。
- `tools/verify_migration.ps1`：检查 Python、源码和测试。

迁移包的完整说明见包内 `MIGRATION_README.md`；当前版本按要求不包含权限管理。

## Electron 单机离线版

Windows 桌面安装包位于：

```text
release/electron/LCMS-Desktop-1.0.0-Setup.exe
```

桌面版内置 Electron、便携 Python 和 ThermoRawFileParser，不要求同事单独安装运行环境。任务数据默认保存在安装目录下的 `data`，可在“设置与诊断”中更换到容量更充足的数据盘。升级或卸载会保留该数据目录。

安装目录中同时提供独立的 `LCMS-MCP.exe`，可直接配置为 Codex 的 STDIO MCP 服务。它与桌面版共用数据目录，不需要设置 TIC 或候选位数参数。

开发和重新构建说明见 `desktop/README.md`。TIC 峰与候选组分数量由质量阈值自动判定，桌面界面不提供固定数量参数。

## 在新工作目录中初始化

解压后，在本目录运行：

```powershell
.\initialize_new_project.ps1 -Destination "D:\NewLCMSProject"
```

或直接将以下三个目录保持为同级目录复制到新项目根目录：

```text
lcms_feature_mvp/
lcms_realdata_platform/
lcms_department_platform/
```

保持同级关系很重要，因为部门任务平台默认从项目根目录查找 `lcms_feature_mvp`。

## 快速验证

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\verify_migration.ps1
```

验证内容包括 Python 语法编译、核心单元测试和关键文件存在性检查。

## 推荐阅读顺序

1. `docs/01_CURRENT_PROJECT_SUMMARY.md`
2. `docs/02_ALGORITHM_REFERENCE.md`
3. `docs/03_DATA_AND_API_CONTRACTS.md`
4. `docs/04_SECONDARY_CHROMATOGRAPHY_EXTENSION.md`
5. `docs/05_OPERATIONS_AND_REPORTS.md`

## 最重要的扩展原则

不要让二级色谱模块重新解释网页状态。新模块应以 `FeatureHandoffRecord` 为输入，至少保存：

- comparison/project/sample 标识
- feature group、parent TIC peak、代表 m/z 和校正后 RT
- 每个样品的 raw area、height、feature RT 边界和 local RT shift
- similarity、difference、ranking、fold、match status
- 原始 Peak-first SQLite 路径或 immutable result id

这样新模块既能追溯原始 scan/XIC，又不会与当前页面实现绑定。
