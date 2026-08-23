# LC-MS/MS 分析平台 Electron 桌面壳

桌面壳负责启动本地 Python 服务、打开任务中心、管理数据目录，并在退出时回收 RAW 转换和分析进程。算法继续由现有 Python 代码执行。

## 开发运行

```powershell
cd .\desktop
npm install
npm start
```

开发模式优先使用 `release/LCMS_Department_Platform/runtime/python.exe`，并直接加载工作区中的 Python 源码。

## 数据目录

默认目录为：

```text
<安装目录>\data
```

可通过桌面版“设置与诊断”页面更改。配置保存在安装目录的 `lcms-desktop-config.json`，独立的 `LCMS-MCP.exe` 会读取同一配置。安装位置必须对当前用户可写。

升级或卸载时，安装器会保留安装目录下的 `data` 和数据目录配置，避免分析结果随程序文件一起删除。

## 独立 MCP

安装后，`LCMS-MCP.exe` 位于桌面主程序同级，可直接作为 Codex 的 STDIO MCP 命令，不需要附加 TIC 或候选位数参数。首次构建前安装构建依赖：

```powershell
..\release\LCMS_Department_Platform\runtime\python.exe -m pip install --target ..\release\build-tools -r requirements-build.txt
```

## 构建安装包

先确认 `release/LCMS_Department_Platform/runtime/python.exe` 存在，然后运行：

```powershell
npm run dist
```

安装包输出到 `release/electron/`。构建会先生成自包含的 `release/mcp/LCMS-MCP.exe`；Python、ThermoRawFileParser 和分析代码通过 `extraResources` 复制，不会放入 `app.asar`。
