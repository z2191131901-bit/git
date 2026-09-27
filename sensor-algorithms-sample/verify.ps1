# 编译和运行本目录全部独立 C 示例；失败即停止。
$ErrorActionPreference = 'Stop'
try {
    $compiler = Get-Command gcc -ErrorAction Stop
    $outputDir = Join-Path ([System.IO.Path]::GetTempPath()) ('sensor-algorithms-' + [Guid]::NewGuid().ToString('N'))
    [System.IO.Directory]::CreateDirectory($outputDir) | Out-Null
    $sources = @(Get-ChildItem -LiteralPath $PSScriptRoot -Recurse -Filter '*.c' -File | Sort-Object FullName)
    if ($sources.Count -eq 0) { throw 'No C examples found.' }
    foreach ($source in $sources) {
        $output = Join-Path $outputDir ($source.BaseName + '.exe')
        Write-Output ('Checking: ' + $source.Name)
        & $compiler.Source -std=c99 -Wall -Wextra -Werror -pedantic $source.FullName -lm -o $output
        if ($LASTEXITCODE -ne 0) { throw ('Compilation failed: ' + $source.Name) }
        & $output
        if ($LASTEXITCODE -ne 0) { throw ('Example check failed: ' + $source.Name) }
    }
    Write-Output ('PASS: ' + $sources.Count + ' examples. Build output: ' + $outputDir)
} catch {
    Write-Error $_
    exit 1
}
