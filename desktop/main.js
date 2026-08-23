"use strict";

const { app, BrowserWindow, dialog, ipcMain, Notification, powerSaveBlocker, session, shell } = require("electron");
const { spawn, spawnSync } = require("node:child_process");
const crypto = require("node:crypto");
const fs = require("node:fs");
const path = require("node:path");
const readline = require("node:readline");

const APP_TITLE = "LC-MS/MS 分析平台";
const BACKEND_START_TIMEOUT_MS = 45_000;
const BACKEND_STOP_TIMEOUT_MS = 12_000;

let mainWindow = null;
let backendProcess = null;
let backendOrigin = "";
let backendToken = "";
let backendLogPath = "";
let dataRoot = "";
let isQuitting = false;
let closePromptOpen = false;
let notificationTimer = null;
let knownTaskStates = new Map();
let sleepBlockerId = null;

function existingPath(candidates) {
  return candidates.find(candidate => candidate && fs.existsSync(candidate)) || "";
}

function getRuntimeLayout() {
  if (app.isPackaged) {
    const platformRoot = path.join(process.resourcesPath, "platform");
    return {
      backendRoot: path.join(platformRoot, "app"),
      pythonExe: path.join(platformRoot, "runtime", "python.exe")
    };
  }
  const workspaceRoot = path.resolve(__dirname, "..");
  return {
    backendRoot: workspaceRoot,
    pythonExe: existingPath([
      path.join(workspaceRoot, "release", "LCMS_Department_Platform", "runtime", "python.exe"),
      path.join(workspaceRoot, ".venv", "Scripts", "python.exe"),
      "python.exe"
    ])
  };
}

function getDataConfigPath() {
  return path.join(getInstallRoot(), "lcms-desktop-config.json");
}

function getInstallRoot() {
  return app.isPackaged ? path.dirname(app.getPath("exe")) : path.resolve(__dirname, "..");
}

function getLegacyDataConfigPath() {
  return path.join(app.getPath("userData"), "desktop-config.json");
}

function readDataConfig(configPath) {
  try {
    const saved = JSON.parse(fs.readFileSync(configPath, "utf8"));
    if (saved && typeof saved.dataRoot === "string" && saved.dataRoot.trim()) return path.resolve(saved.dataRoot);
  } catch (_error) {
    // Missing or malformed configuration falls through to the next source.
  }
  return "";
}

function readDataRoot() {
  const configuredByEnvironment = String(process.env.LCMS_DATA_ROOT || "").trim();
  if (configuredByEnvironment) return path.resolve(configuredByEnvironment);
  const configured = readDataConfig(getDataConfigPath());
  if (configured) return configured;
  const legacyConfigured = readDataConfig(getLegacyDataConfigPath());
  if (legacyConfigured) {
    saveDataRoot(legacyConfigured);
    return legacyConfigured;
  }
  return path.join(getInstallRoot(), "data");
}

function saveDataRoot(nextRoot) {
  fs.mkdirSync(getInstallRoot(), { recursive: true });
  fs.writeFileSync(getDataConfigPath(), JSON.stringify({ dataRoot: nextRoot }, null, 2), "utf8");
}

function ensureWritableDirectory(target) {
  fs.mkdirSync(target, { recursive: true });
  const probe = path.join(target, `.lcms-write-test-${process.pid}-${crypto.randomBytes(6).toString("hex")}`);
  try {
    fs.writeFileSync(probe, "ok", { encoding: "utf8", flag: "wx" });
  } finally {
    if (fs.existsSync(probe)) fs.unlinkSync(probe);
  }
}

function sendStartupStatus(message, kind = "info") {
  if (mainWindow && !mainWindow.isDestroyed()) {
    mainWindow.webContents.send("desktop:backend-status", { message, kind });
  }
}

function appendBackendLog(prefix, line) {
  if (!backendLogPath) return;
  fs.appendFileSync(backendLogPath, `[${new Date().toISOString()}] ${prefix}${line}\n`, "utf8");
}

async function waitForBackend(origin) {
  const deadline = Date.now() + BACKEND_START_TIMEOUT_MS;
  while (Date.now() < deadline) {
    try {
      const response = await fetch(`${origin}/api/status`, { cache: "no-store" });
      if (response.ok) return await response.json();
    } catch (_error) {
      // The process may have printed its port before the HTTP server is ready.
    }
    await new Promise(resolve => setTimeout(resolve, 250));
  }
  throw new Error("本地分析服务启动超时，请查看桌面日志。 ");
}

async function startBackend() {
  const runtime = getRuntimeLayout();
  const serverPath = path.join(runtime.backendRoot, "lcms_department_platform", "server.py");
  const parserPath = path.join(runtime.backendRoot, "lcms_department_platform", "tools", "ThermoRawFileParser", "ThermoRawFileParser.exe");
  for (const requiredPath of [runtime.pythonExe, serverPath, parserPath]) {
    if (!requiredPath || !fs.existsSync(requiredPath)) {
      throw new Error(`缺少运行文件：${requiredPath || "Python 运行时"}`);
    }
  }

  fs.mkdirSync(dataRoot, { recursive: true });
  const desktopLogDir = path.join(dataRoot, "logs", "desktop");
  fs.mkdirSync(desktopLogDir, { recursive: true });
  backendLogPath = path.join(desktopLogDir, "backend.log");
  backendToken = crypto.randomBytes(32).toString("hex");

  const args = [
    serverPath,
    "--host", "127.0.0.1",
    "--port", "0",
    "--root", dataRoot,
    "--parser-path", parserPath,
    "--python", runtime.pythonExe,
    "--control-token", backendToken
  ];

  sendStartupStatus("正在启动 Python 分析服务…");
  backendProcess = spawn(runtime.pythonExe, args, {
    cwd: runtime.backendRoot,
    windowsHide: true,
    shell: false,
    env: { ...process.env, PYTHONUTF8: "1", PYTHONUNBUFFERED: "1" },
    stdio: ["ignore", "pipe", "pipe"]
  });

  const port = await new Promise((resolve, reject) => {
    let settled = false;
    const timeout = setTimeout(() => {
      if (!settled) {
        settled = true;
        reject(new Error("未收到本地服务端口，启动失败。"));
      }
    }, BACKEND_START_TIMEOUT_MS);
    const output = readline.createInterface({ input: backendProcess.stdout });
    output.on("line", line => {
      appendBackendLog("OUT ", line);
      const match = /^LCMS_PORT=(\d+)$/.exec(line.trim());
      if (match && !settled) {
        settled = true;
        clearTimeout(timeout);
        resolve(Number(match[1]));
      }
    });
    const errors = readline.createInterface({ input: backendProcess.stderr });
    errors.on("line", line => appendBackendLog("ERR ", line));
    backendProcess.once("error", error => {
      appendBackendLog("SPAWN ", error.stack || String(error));
      if (!settled) {
        settled = true;
        clearTimeout(timeout);
        reject(error);
      }
    });
    backendProcess.once("exit", code => {
      appendBackendLog("EXIT ", `code=${code}`);
      if (!settled) {
        settled = true;
        clearTimeout(timeout);
        reject(new Error(`本地分析服务提前退出，代码 ${code}。`));
      }
      backendProcess = null;
      backendOrigin = "";
    });
  });

  backendOrigin = `http://127.0.0.1:${port}`;
  await waitForBackend(backendOrigin);
  sendStartupStatus("分析服务已就绪，正在打开任务中心…");
  return backendOrigin;
}

async function backendJson(route, options = {}) {
  if (!backendOrigin) throw new Error("本地分析服务尚未启动。");
  const response = await fetch(`${backendOrigin}${route}`, { cache: "no-store", ...options });
  if (!response.ok) throw new Error(await response.text());
  return await response.json();
}

async function activeTasks() {
  try {
    const payload = await backendJson("/api/tasks");
    return (payload.tasks || []).filter(task => ["waiting", "running", "canceling"].includes(task.status));
  } catch (_error) {
    return [];
  }
}

async function waitForProcessExit(child, timeoutMs) {
  if (!child || child.exitCode !== null) return true;
  return await new Promise(resolve => {
    const timeout = setTimeout(() => resolve(false), timeoutMs);
    child.once("exit", () => {
      clearTimeout(timeout);
      resolve(true);
    });
  });
}

async function stopBackend() {
  const child = backendProcess;
  if (!child || child.exitCode !== null) return;
  try {
    await fetch(`${backendOrigin}/api/control/shutdown`, {
      method: "POST",
      headers: { "X-LCMS-Control-Token": backendToken }
    });
  } catch (error) {
    appendBackendLog("SHUTDOWN ", String(error));
  }
  if (await waitForProcessExit(child, BACKEND_STOP_TIMEOUT_MS)) return;
  if (process.platform === "win32") {
    spawnSync("taskkill.exe", ["/PID", String(child.pid), "/T", "/F"], { windowsHide: true });
  } else {
    child.kill("SIGKILL");
  }
}

async function quitApplication() {
  if (isQuitting) return;
  isQuitting = true;
  if (notificationTimer) clearInterval(notificationTimer);
  if (sleepBlockerId !== null && powerSaveBlocker.isStarted(sleepBlockerId)) powerSaveBlocker.stop(sleepBlockerId);
  await stopBackend();
  if (mainWindow && !mainWindow.isDestroyed()) mainWindow.destroy();
  app.quit();
}

async function handleWindowClose(event) {
  if (isQuitting || closePromptOpen) {
    if (!isQuitting) event.preventDefault();
    return;
  }
  event.preventDefault();
  closePromptOpen = true;
  try {
    const tasks = await activeTasks();
    if (tasks.length) {
      const result = await dialog.showMessageBox(mainWindow, {
        type: "warning",
        title: "分析任务仍在运行",
        message: `当前有 ${tasks.length} 个任务正在运行或等待。`,
        detail: "退出会取消正在运行的转换或分析进程。",
        buttons: ["最小化并继续运行", "取消任务并退出", "返回"],
        defaultId: 0,
        cancelId: 2,
        noLink: true
      });
      if (result.response === 0) mainWindow.minimize();
      if (result.response === 1) await quitApplication();
      return;
    }
    await quitApplication();
  } finally {
    closePromptOpen = false;
  }
}

function safeTaskDirectory(taskId) {
  if (!/^[\p{L}\p{N}_-]+$/u.test(taskId)) return "";
  const jobsRoot = path.resolve(dataRoot, "jobs");
  const target = path.resolve(jobsRoot, taskId);
  return target.startsWith(`${jobsRoot}${path.sep}`) ? target : "";
}

async function pollTaskNotifications() {
  try {
    const payload = await backendJson("/api/tasks");
    const tasks = payload.tasks || [];
    const hasActive = tasks.some(task => ["waiting", "running", "canceling"].includes(task.status));
    if (hasActive && sleepBlockerId === null) sleepBlockerId = powerSaveBlocker.start("prevent-app-suspension");
    if (!hasActive && sleepBlockerId !== null) {
      if (powerSaveBlocker.isStarted(sleepBlockerId)) powerSaveBlocker.stop(sleepBlockerId);
      sleepBlockerId = null;
    }
    for (const task of tasks) {
      const previous = knownTaskStates.get(task.task_id);
      if (previous && previous !== task.status && ["finished", "failed", "canceled"].includes(task.status)) {
        const body = task.status === "finished" ? "分析已完成，可以查看报告。" : task.status === "failed" ? "分析失败，请查看任务日志。" : "任务已取消。";
        new Notification({ title: task.task_name || APP_TITLE, body }).show();
      }
      knownTaskStates.set(task.task_id, task.status);
    }
  } catch (_error) {
    // Temporary polling failures are surfaced by the web UI status panel.
  }
}

function registerIpcHandlers() {
  ipcMain.handle("desktop:get-app-info", () => ({
    version: app.getVersion(),
    dataRoot,
    backendOrigin,
    backendLogPath
  }));
  ipcMain.handle("desktop:open-data-directory", async () => {
    fs.mkdirSync(dataRoot, { recursive: true });
    const error = await shell.openPath(dataRoot);
    return { ok: !error, error };
  });
  ipcMain.handle("desktop:open-task-directory", async (_event, taskId) => {
    const target = safeTaskDirectory(taskId);
    if (!target || !fs.existsSync(target)) return { ok: false, error: "任务目录不存在。" };
    const error = await shell.openPath(target);
    return { ok: !error, error };
  });
  ipcMain.handle("desktop:choose-data-directory", async () => {
    if ((await activeTasks()).length) return { ok: false, error: "有任务正在运行或排队，暂时不能更改数据目录。" };
    const selection = await dialog.showOpenDialog(mainWindow, {
      title: "选择 LC-MS 数据目录",
      defaultPath: dataRoot,
      properties: ["openDirectory", "createDirectory"]
    });
    if (selection.canceled || !selection.filePaths[0]) return { ok: false, canceled: true };
    const nextRoot = path.resolve(selection.filePaths[0]);
    try {
      ensureWritableDirectory(nextRoot);
    } catch (error) {
      return { ok: false, error: `所选目录不可写：${error.message || error}` };
    }
    saveDataRoot(nextRoot);
    await dialog.showMessageBox(mainWindow, {
      type: "info",
      title: "数据目录已更新",
      message: "应用将重新启动并使用新的数据目录。",
      detail: nextRoot,
      buttons: ["重新启动"],
      noLink: true
    });
    app.relaunch();
    await quitApplication();
    return { ok: true, dataRoot: nextRoot };
  });
}

async function createMainWindow() {
  mainWindow = new BrowserWindow({
    title: APP_TITLE,
    icon: path.join(__dirname, "assets", "icon.png"),
    width: 1400,
    height: 900,
    minWidth: 1180,
    minHeight: 720,
    show: false,
    backgroundColor: "#edf4fb",
    autoHideMenuBar: true,
    webPreferences: {
      preload: path.join(__dirname, "preload.js"),
      contextIsolation: true,
      nodeIntegration: false,
      sandbox: true,
      webSecurity: true
    }
  });
  mainWindow.webContents.setWindowOpenHandler(() => ({ action: "deny" }));
  mainWindow.webContents.on("will-navigate", (event, url) => {
    if (url.startsWith("file:") || (backendOrigin && url.startsWith(`${backendOrigin}/`))) return;
    event.preventDefault();
  });
  mainWindow.on("close", handleWindowClose);
  mainWindow.once("ready-to-show", () => mainWindow.show());
  await mainWindow.loadFile(path.join(__dirname, "startup.html"));
  mainWindow.show();
}

async function bootstrap() {
  try {
    dataRoot = readDataRoot();
    ensureWritableDirectory(dataRoot);
    registerIpcHandlers();
    session.defaultSession.setPermissionRequestHandler((_webContents, _permission, callback) => callback(false));
    await createMainWindow();
    const origin = await startBackend();
    await mainWindow.loadURL(`${origin}/`);
    await pollTaskNotifications();
    notificationTimer = setInterval(pollTaskNotifications, 5000);
  } catch (error) {
    sendStartupStatus(`${error.message || error}\n日志：${backendLogPath || "尚未创建"}`, "error");
    dialog.showErrorBox(
      "LC-MS/MS 分析平台启动失败",
      `${String(error.stack || error)}\n\n默认数据目录位于安装目录下。若安装目录不可写，请将程序安装到当前用户可写目录。`
    );
  }
}

const hasSingleInstanceLock = app.requestSingleInstanceLock();
if (!hasSingleInstanceLock) {
  app.quit();
} else {
  app.on("second-instance", () => {
    if (!mainWindow) return;
    if (mainWindow.isMinimized()) mainWindow.restore();
    mainWindow.show();
    mainWindow.focus();
  });
  app.whenReady().then(bootstrap);
  app.on("window-all-closed", () => {
    if (process.platform !== "darwin" && !isQuitting) void quitApplication();
  });
  app.on("before-quit", event => {
    if (isQuitting) return;
    event.preventDefault();
    void quitApplication();
  });
}
