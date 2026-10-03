param(
    [ValidateSet('build','run','bench','check')][string]$Action = 'run',
    [ValidateRange(1,1000)][int]$Count = 50,
    [string]$Compiler = '',
    [string]$InputPath = ''
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
$buildPath = Join-Path $projectRoot 'build'
$resultDir = Join-Path $projectRoot 'results'
$exePath = Join-Path $buildPath 'portfolio.exe'
if ($Action -eq 'build') {
    if (-not $Compiler) {
        $command = Get-Command zig -ErrorAction SilentlyContinue
        if ($command) { $Compiler = $command.Source }
    }
    if (-not $Compiler) { throw 'Specify -Compiler with the path to zig.exe, or add Zig to PATH.' }
    [void][IO.Directory]::CreateDirectory($buildPath)
    $env:ZIG_LOCAL_CACHE_DIR = Join-Path $buildPath 'zig-cache'
    $env:ZIG_GLOBAL_CACHE_DIR = Join-Path $buildPath 'zig-global-cache'
    & $Compiler c++ -O3 -std=c++17 -mcpu=baseline -ffp-contract=off -flto -s '-Wl,--subsystem,windows' (Join-Path $projectRoot 'src/portfolio.cpp') -o $exePath
    if ($LASTEXITCODE -ne 0) { throw 'Build failed.' }
    Write-Host "Built: $exePath"
    return
}
if (-not (Test-Path -LiteralPath $exePath)) { throw 'Run -Action build first.' }
if (-not $InputPath) { $InputPath = Join-Path $projectRoot 'examples/synthetic.xls' }
$InputPath = (Resolve-Path -LiteralPath $InputPath).Path
if ($Action -eq 'check') {
    & $exePath --self-test $InputPath | Out-Host
    if ($LASTEXITCODE -ne 0) { throw 'Self-test failed.' }
    return
}
[void][IO.Directory]::CreateDirectory($resultDir)
$resultPath = Join-Path $resultDir 'results.json'
$iterations = if ($Action -eq 'bench') { $Count } else { 1 }
$samples = @()
for ($i = 0; $i -lt $iterations; $i++) {
    $startInfo = New-Object Diagnostics.ProcessStartInfo
    $startInfo.FileName = $exePath
    $startInfo.Arguments = '"' + $InputPath + '" "' + $resultPath + '"'
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $process = New-Object Diagnostics.Process
    $process.StartInfo = $startInfo
    $watch = [Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $process.WaitForExit()
    $output = $stdoutTask.GetAwaiter().GetResult()
    $errorOutput = $stderrTask.GetAwaiter().GetResult()
    $exitCode = $process.ExitCode
    $watch.Stop()
    $process.Dispose()
    if ($exitCode -ne 0) { throw "Solver failed: $errorOutput" }
    $samples += $watch.Elapsed.TotalMilliseconds
    if ($Action -eq 'run') { $output | Out-Host }
}
if ($Action -eq 'bench') {
    $sorted = @($samples | Sort-Object)
    $median = if ($iterations % 2) { $sorted[[int][Math]::Floor($iterations/2)] } else { ($sorted[$iterations/2-1] + $sorted[$iterations/2])/2 }
    [ordered]@{
        count = $iterations
        median_ms = $median
        min_ms = $sorted[0]
        max_ms = $sorted[-1]
        samples_ms = $samples
        boundary = 'PowerShell invocation through captured output and process exit; includes input, solves and export; no warmup discarded.'
    } | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $resultDir 'benchmark.json') -Encoding UTF8
    Write-Host ('Process median: {0:F3} ms ({1} samples)' -f $median,$iterations)
}
Write-Host "Results: $resultPath"
