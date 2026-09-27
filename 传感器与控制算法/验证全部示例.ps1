# 中文文件保留原名，为旧版 MinGW 建立临时 ASCII 编译副本。
param([string]$示例='')
$ErrorActionPreference='Stop'
try {
 $compiler=Get-Command gcc -ErrorAction Stop
 $outputDir=Join-Path ([IO.Path]::GetTempPath()) ('sensor-algorithms-'+[Guid]::NewGuid().ToString('N'))
 [IO.Directory]::CreateDirectory($outputDir) | Out-Null
 $files=@(Get-ChildItem -LiteralPath $PSScriptRoot -Recurse -File | Where-Object {$_.Extension -in '.c','.h'} | Sort-Object FullName)
 $mapping=@{}
 $number=0
 foreach($file in $files) {$mapping[$file.FullName]=('unit_'+$number+$file.Extension); $number++}
 foreach($file in $files) {
  $content=[IO.File]::ReadAllText($file.FullName)
  # 只改临时副本的本地 include，不改系统头文件和算法逻辑。
  $content=[regex]::Replace($content,'(?m)^(\s*#\s*include\s*)"([^"]+)"',{
   param($match)
   $header=[IO.Path]::GetFullPath((Join-Path $file.DirectoryName $match.Groups[2].Value))
   if(-not $mapping.ContainsKey($header)) {throw ('Unknown header: '+$header)}
   return $match.Groups[1].Value+'"'+$mapping[$header]+'"'
  })
  [IO.File]::WriteAllText((Join-Path $outputDir $mapping[$file.FullName]),$content,[Text.UTF8Encoding]::new($false))
 }
 $sources=@($files | Where-Object {$_.Extension -eq '.c' -and (!$示例 -or $_.Name -eq $示例)})
 if($sources.Count -eq 0) {throw 'No matching C examples.'}
 foreach($source in $sources) {
  $inputPath=Join-Path $outputDir $mapping[$source.FullName]
  $output=Join-Path $outputDir ([IO.Path]::GetFileNameWithoutExtension($inputPath)+'.exe')
  Write-Output ('Checking: '+$source.Name)
  & $compiler.Source -std=c99 -Wall -Wextra -Werror -pedantic $inputPath -lm -o $output
  if($LASTEXITCODE -ne 0) {throw ('Compilation failed: '+$source.Name)}
  & $output
  if($LASTEXITCODE -ne 0) {throw ('Example check failed: '+$source.Name)}
 }
 Write-Output ('PASS: '+$sources.Count+' examples. Build output: '+$outputDir)
} catch {Write-Error $_; exit 1}
