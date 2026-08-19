param(
    [string]$InstallDir = "",
    [switch]$NoLaunch
)

$ErrorActionPreference = "Stop"
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing

$SourceRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$DefaultInstallDir = Join-Path $env:LOCALAPPDATA "LCMS Department Platform"

function Get-PayloadRoot {
    $payloadArchive = Join-Path $SourceRoot "LCMS-Payload.zip"
    if (-not (Test-Path -LiteralPath $payloadArchive -PathType Leaf)) {
        return $SourceRoot
    }

    $payloadRoot = Join-Path $SourceRoot "payload_extract"
    if (Test-Path -LiteralPath $payloadRoot) {
        Remove-Item -LiteralPath $payloadRoot -Recurse -Force
    }
    Expand-Archive -LiteralPath $payloadArchive -DestinationPath $payloadRoot -Force
    return $payloadRoot
}

function New-DesktopShortcut([string]$TargetInstallDir) {
    $desktop = [Environment]::GetFolderPath("Desktop")
    if (-not $desktop) { return }

    try {
        $shortcutPath = Join-Path $desktop "LC-MS Department Platform.lnk"
        $shell = New-Object -ComObject WScript.Shell
        $shortcut = $shell.CreateShortcut($shortcutPath)
        $shortcut.TargetPath = Join-Path $TargetInstallDir "Start-LCMS-Platform.cmd"
        $shortcut.WorkingDirectory = $TargetInstallDir
        $shortcut.Description = "LC-MS Department Platform"
        $shortcut.Save()
    } catch {
        Write-Warning ("Could not create a desktop shortcut: {0}" -f $_.Exception.Message)
    }
}

function Invoke-PlatformInstall(
    [string]$TargetInstallDir,
    [scriptblock]$ProgressCallback
) {
    $TargetInstallDir = [System.IO.Path]::GetFullPath($TargetInstallDir)
    $payloadRoot = Get-PayloadRoot
    try {
        if ($ProgressCallback) { & $ProgressCallback 0 "Preparing installation..." $true }

        New-Item -ItemType Directory -Force -Path $TargetInstallDir | Out-Null
        $sourceFiles = @(Get-ChildItem -LiteralPath $payloadRoot -Recurse -File)
        if (-not $sourceFiles.Count) {
            throw "The installation payload is empty."
        }

        $total = $sourceFiles.Count
        $index = 0
        if ($ProgressCallback) { & $ProgressCallback 0 "Copying platform files..." $false }

        foreach ($sourceFile in $sourceFiles) {
            $index++
            $relative = $sourceFile.FullName.Substring($payloadRoot.Length).TrimStart("\", "/")
            $destination = Join-Path $TargetInstallDir $relative
            $destinationParent = Split-Path -Parent $destination
            New-Item -ItemType Directory -Force -Path $destinationParent | Out-Null
            Copy-Item -LiteralPath $sourceFile.FullName -Destination $destination -Force

            if ($ProgressCallback) {
                $percent = [int][Math]::Floor(($index * 100.0) / $total)
                & $ProgressCallback $percent ("Installing {0} of {1}: {2}" -f $index, $total, $relative) $false
            }
        }

        New-DesktopShortcut $TargetInstallDir
        if ($ProgressCallback) { & $ProgressCallback 100 "Installation completed." $false }

        if (-not $NoLaunch) {
            $startScript = Join-Path $TargetInstallDir "Start-LCMS-Platform.cmd"
            Start-Process -FilePath $startScript -WorkingDirectory $TargetInstallDir
        }

        return $TargetInstallDir
    }
    finally {
        if ($payloadRoot -ne $SourceRoot -and (Test-Path -LiteralPath $payloadRoot)) {
            Remove-Item -LiteralPath $payloadRoot -Recurse -Force -ErrorAction SilentlyContinue
        }
    }
}

# A supplied path is used for unattended installation and test automation.
if ($InstallDir) {
    Invoke-PlatformInstall -TargetInstallDir $InstallDir
    Write-Host "LC-MS Department Platform installed to: $InstallDir"
    exit 0
}

$form = New-Object System.Windows.Forms.Form
$form.Text = "LC-MS Department Platform Installer"
$form.StartPosition = "CenterScreen"
$form.Size = New-Object System.Drawing.Size(720, 430)
$form.MinimumSize = New-Object System.Drawing.Size(720, 430)
$form.FormBorderStyle = "FixedDialog"
$form.MaximizeBox = $false
$form.MinimizeBox = $false

$title = New-Object System.Windows.Forms.Label
$title.Text = "LC-MS Department Platform"
$title.Font = New-Object System.Drawing.Font("Segoe UI", 16, [System.Drawing.FontStyle]::Bold)
$title.Location = New-Object System.Drawing.Point(28, 22)
$title.AutoSize = $true
$form.Controls.Add($title)

$subtitle = New-Object System.Windows.Forms.Label
$subtitle.Text = "Choose an installation folder, then click Install."
$subtitle.Location = New-Object System.Drawing.Point(30, 58)
$subtitle.AutoSize = $true
$form.Controls.Add($subtitle)

$pathLabel = New-Object System.Windows.Forms.Label
$pathLabel.Text = "Installation folder:"
$pathLabel.Location = New-Object System.Drawing.Point(30, 105)
$pathLabel.AutoSize = $true
$form.Controls.Add($pathLabel)

$pathBox = New-Object System.Windows.Forms.TextBox
$pathBox.Text = $DefaultInstallDir
$pathBox.Location = New-Object System.Drawing.Point(30, 130)
$pathBox.Size = New-Object System.Drawing.Size(545, 28)
$form.Controls.Add($pathBox)

$browseButton = New-Object System.Windows.Forms.Button
$browseButton.Text = "Browse..."
$browseButton.Location = New-Object System.Drawing.Point(585, 129)
$browseButton.Size = New-Object System.Drawing.Size(100, 30)
$form.Controls.Add($browseButton)

$launchCheck = New-Object System.Windows.Forms.CheckBox
$launchCheck.Text = "Start the platform after installation"
$launchCheck.Checked = $true
$launchCheck.Location = New-Object System.Drawing.Point(30, 175)
$launchCheck.AutoSize = $true
$form.Controls.Add($launchCheck)

$statusLabel = New-Object System.Windows.Forms.Label
$statusLabel.Text = "Ready to install."
$statusLabel.Location = New-Object System.Drawing.Point(30, 225)
$statusLabel.Size = New-Object System.Drawing.Size(655, 42)
$statusLabel.AutoEllipsis = $true
$form.Controls.Add($statusLabel)

$progressBar = New-Object System.Windows.Forms.ProgressBar
$progressBar.Minimum = 0
$progressBar.Maximum = 100
$progressBar.Value = 0
$progressBar.Location = New-Object System.Drawing.Point(30, 275)
$progressBar.Size = New-Object System.Drawing.Size(655, 24)
$form.Controls.Add($progressBar)

$installButton = New-Object System.Windows.Forms.Button
$installButton.Text = "Install"
$installButton.DialogResult = [System.Windows.Forms.DialogResult]::None
$installButton.Location = New-Object System.Drawing.Point(485, 335)
$installButton.Size = New-Object System.Drawing.Size(95, 32)
$form.Controls.Add($installButton)
$form.AcceptButton = $installButton

$cancelButton = New-Object System.Windows.Forms.Button
$cancelButton.Text = "Cancel"
$cancelButton.Location = New-Object System.Drawing.Point(590, 335)
$cancelButton.Size = New-Object System.Drawing.Size(95, 32)
$form.Controls.Add($cancelButton)
$form.CancelButton = $cancelButton

$browseButton.Add_Click({
    $dialog = New-Object System.Windows.Forms.FolderBrowserDialog
    $dialog.Description = "Select the LC-MS Department Platform installation folder"
    $dialog.SelectedPath = $pathBox.Text
    if ($dialog.ShowDialog() -eq [System.Windows.Forms.DialogResult]::OK) {
        $pathBox.Text = $dialog.SelectedPath
    }
    $dialog.Dispose()
})

$cancelButton.Add_Click({ $form.Close() })

$installButton.Add_Click({
    $selectedPath = $pathBox.Text.Trim()
    if (-not $selectedPath) {
        [System.Windows.Forms.MessageBox]::Show("Please choose an installation folder.", "Installation", "OK", "Warning") | Out-Null
        return
    }

    $installButton.Enabled = $false
    $browseButton.Enabled = $false
    $cancelButton.Enabled = $false
    $pathBox.Enabled = $false
    $launchCheck.Enabled = $false
    $form.UseWaitCursor = $true

    $callback = {
        param($percent, $message, $marquee)
        $statusLabel.Text = $message
        if ($marquee) {
            $progressBar.Style = "Marquee"
        } else {
            if ($progressBar.Style -ne "Blocks") { $progressBar.Style = "Blocks" }
            $progressBar.Value = [Math]::Max(0, [Math]::Min(100, [int]$percent))
        }
        [System.Windows.Forms.Application]::DoEvents()
    }

    try {
        # The selected launch setting is applied only after the copy succeeds.
        $NoLaunch = -not $launchCheck.Checked
        Invoke-PlatformInstall -TargetInstallDir $selectedPath -ProgressCallback $callback | Out-Null
        $form.UseWaitCursor = $false
        $progressBar.Style = "Blocks"
        $progressBar.Value = 100
        $statusLabel.Text = "Installation completed."
        [System.Windows.Forms.MessageBox]::Show(
            "LC-MS Department Platform was installed successfully.`n`nLocation: $selectedPath`n`nWeb interface: http://127.0.0.1:8770/",
            "Installation complete",
            "OK",
            "Information"
        ) | Out-Null
        $form.Close()
    } catch {
        $form.UseWaitCursor = $false
        $installButton.Enabled = $true
        $browseButton.Enabled = $true
        $cancelButton.Enabled = $true
        $pathBox.Enabled = $true
        $launchCheck.Enabled = $true
        $progressBar.Style = "Blocks"
        [System.Windows.Forms.MessageBox]::Show($_.Exception.Message, "Installation failed", "OK", "Error") | Out-Null
    }
})

[void]$form.ShowDialog()
$form.Dispose()
