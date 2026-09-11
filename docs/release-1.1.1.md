# LCMS Desktop 1.1.1

修复 1.1.0 定向 N-糖肽搜索的酶参数遗漏：任务所选酶现在从 CLI 传入
`search_feature_glycopeptide_scans`，再传入糖肽候选构建和 `_decoy`。
非 Trypsin 糖肽诱饵不再使用默认 Trypsin 末端保留规则。
默认未指定酶时仍为 Trypsin，不改变旧任务的默认行为。

## 验证

- 新增五项回归测试，覆盖全部十个酶选项、固定 CAM 开启/关闭，检查目标与诱饵
  的对应关系、酶切末端、修饰残基及质量一致性。
- Asp-N、Lys-N、Glu-C 合成糖肽谱通过实际定向搜索，正确命中目标序列和 G0F 糖型。
- 验证 CLI 调用、公开搜索函数和候选生成之间的参数传递，以及默认 Trypsin 的兼容性。
- 分析侧 85 项、平台/MCP 侧 7 项测试通过，共 92 项。

Electron 和独立 MCP.exe 版本同步为 1.1.1。
安装包：`release/electron/LCMS-Desktop-1.1.1-Setup.exe`（Windows x64，未数字签名）。
安装前关闭桌面程序并备份 `data` 和配置文件。没有自动替换正在运行的旧版，
也没有覆盖历史任务；非 Trypsin 的已有糖肽结果需重新分析才应用本次修复。

修复范围仅为酶参数传递，不意味着所有实验条件已经充分验证。
非还原/部分烷基化、探索候选及标记肽的既有限制仍见
[制样与搜索说明](unimod-and-enzymes.md)。
