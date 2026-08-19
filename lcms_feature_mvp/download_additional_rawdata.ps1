param(
    [string]$OutputDir = "data/raw/zenodo_5005513",
    [int]$Retry = 8,
    [int]$RetryDelaySeconds = 10,
    [switch]$Force
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$ResolvedOutput = Join-Path $Root $OutputDir
New-Item -ItemType Directory -Force -Path $ResolvedOutput | Out-Null
Add-Type -AssemblyName System.Net.Http

$RecordUrl = "https://zenodo.org/api/records/5005513"
$Targets = @(
    "SCX-HPLC-MS_Intact_MabThera_Untreated_3.raw",
    "SCX-HPLC-MS_Intact_Reditux_Untreated_3.raw",
    "SCX-HPLC-MS_Intact_MabThera_Glycated.raw",
    "SCX-HPLC-MS_Intact_Reditux_Glycated.raw",
    "SCX-HPLC-MS_Intact_MabThera_Deamidated.raw",
    "SCX-HPLC-MS_Intact_Reditux_Deamidated.raw"
)

function Get-FileMd5 {
    param([string]$Path)
    return (Get-FileHash -LiteralPath $Path -Algorithm MD5).Hash.ToLowerInvariant()
}

function Test-DownloadedFile {
    param(
        [string]$Path,
        [int64]$ExpectedSize,
        [string]$ExpectedMd5
    )
    if (-not (Test-Path -LiteralPath $Path)) {
        return $false
    }
    $item = Get-Item -LiteralPath $Path
    if ($item.Length -ne $ExpectedSize) {
        return $false
    }
    if ($ExpectedMd5 -and (Get-FileMd5 -Path $Path) -ne $ExpectedMd5) {
        return $false
    }
    return $true
}

function Invoke-ResumableDownload {
    param(
        [string]$Uri,
        [string]$Destination,
        [int64]$ExistingBytes
    )
    $handler = [System.Net.Http.HttpClientHandler]::new()
    $client = [System.Net.Http.HttpClient]::new($handler)
    try {
        $request = [System.Net.Http.HttpRequestMessage]::new([System.Net.Http.HttpMethod]::Get, $Uri)
        if ($ExistingBytes -gt 0) {
            $request.Headers.Range = [System.Net.Http.Headers.RangeHeaderValue]::new($ExistingBytes, $null)
        }
        $response = $client.SendAsync($request, [System.Net.Http.HttpCompletionOption]::ResponseHeadersRead).GetAwaiter().GetResult()
        if (-not $response.IsSuccessStatusCode) {
            throw "HTTP $([int]$response.StatusCode) $($response.ReasonPhrase)"
        }
        $append = $ExistingBytes -gt 0 -and $response.StatusCode -eq [System.Net.HttpStatusCode]::PartialContent
        if ($ExistingBytes -gt 0 -and -not $append) {
            Write-Host "Server did not honor Range; restarting $Destination"
        }
        $mode = if ($append) { [System.IO.FileMode]::Append } else { [System.IO.FileMode]::Create }
        $inputStream = $response.Content.ReadAsStreamAsync().GetAwaiter().GetResult()
        $outputStream = [System.IO.FileStream]::new($Destination, $mode, [System.IO.FileAccess]::Write, [System.IO.FileShare]::None)
        try {
            $buffer = New-Object byte[] (1024 * 1024)
            while ($true) {
                $read = $inputStream.Read($buffer, 0, $buffer.Length)
                if ($read -le 0) {
                    break
                }
                $outputStream.Write($buffer, 0, $read)
            }
        } finally {
            $outputStream.Dispose()
            $inputStream.Dispose()
        }
    } finally {
        $client.Dispose()
    }
}

Write-Host "Reading Zenodo record $RecordUrl"
$record = Invoke-RestMethod -Uri $RecordUrl

foreach ($name in $Targets) {
    $file = $record.files | Where-Object { $_.key -eq $name } | Select-Object -First 1
    if (-not $file) {
        throw "Missing in Zenodo record: $name"
    }
    $expectedMd5 = ""
    if ($file.checksum -like "md5:*") {
        $expectedMd5 = $file.checksum.Substring(4).ToLowerInvariant()
    }
    $dest = Join-Path $ResolvedOutput $name
    if (-not $Force -and (Test-DownloadedFile -Path $dest -ExpectedSize ([int64]$file.size) -ExpectedMd5 $expectedMd5)) {
        Write-Host "Complete: $name"
        continue
    }
    if ($Force -and (Test-Path -LiteralPath $dest)) {
        Remove-Item -LiteralPath $dest
    }

    for ($attempt = 1; $attempt -le $Retry; $attempt++) {
        $currentSize = if (Test-Path -LiteralPath $dest) { (Get-Item -LiteralPath $dest).Length } else { 0 }
        Write-Host "Downloading/resuming $name attempt $attempt/$Retry current=$currentSize expected=$($file.size)"
        try {
            Invoke-ResumableDownload -Uri $file.links.self -Destination $dest -ExistingBytes $currentSize
        } catch {
            Write-Host "Attempt failed: $($_.Exception.Message)"
        }
        if (Test-DownloadedFile -Path $dest -ExpectedSize ([int64]$file.size) -ExpectedMd5 $expectedMd5) {
            Write-Host "Verified: $name"
            break
        }
        if ($attempt -eq $Retry) {
            $actualSize = if (Test-Path -LiteralPath $dest) { (Get-Item -LiteralPath $dest).Length } else { 0 }
            throw "Download failed verification for $name current=$actualSize expected=$($file.size)"
        }
        Start-Sleep -Seconds $RetryDelaySeconds
    }
}

Write-Host "Additional RAW download complete: $ResolvedOutput"
Get-ChildItem -LiteralPath $ResolvedOutput -Filter "*.raw" | Sort-Object Name | Select-Object Name,Length
