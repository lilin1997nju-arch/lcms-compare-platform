将 ThermoRawFileParser.exe 放在此目录后，部门平台即可处理 Thermo .raw 文件。

也可以启动时通过：

  powershell -ExecutionPolicy Bypass -File ..\start_lcms_department_platform.ps1 -ParserPath C:\path\ThermoRawFileParser.exe

如果只上传 mzML，不需要提供转换器。
