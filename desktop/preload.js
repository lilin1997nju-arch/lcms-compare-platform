"use strict";

const { contextBridge, ipcRenderer } = require("electron");

contextBridge.exposeInMainWorld("lcmsDesktop", Object.freeze({
  isDesktop: true,
  getAppInfo: () => ipcRenderer.invoke("desktop:get-app-info"),
  openDataDirectory: () => ipcRenderer.invoke("desktop:open-data-directory"),
  openTaskDirectory: taskId => ipcRenderer.invoke("desktop:open-task-directory", String(taskId || "")),
  chooseDataDirectory: () => ipcRenderer.invoke("desktop:choose-data-directory"),
  onBackendStatus: callback => {
    if (typeof callback !== "function") return () => {};
    const listener = (_event, payload) => callback(payload);
    ipcRenderer.on("desktop:backend-status", listener);
    return () => ipcRenderer.removeListener("desktop:backend-status", listener);
  }
}));
