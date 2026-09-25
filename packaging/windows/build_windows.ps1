param([switch]$SkipGates)
$ErrorActionPreference = 'Stop'
$root = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..\..'))
$python = Join-Path $root '.venv\Scripts\python.exe'
$iscc = Join-Path $env:LOCALAPPDATA 'Programs\Inno Setup 6\ISCC.exe'
$innoVersion = (Get-ItemProperty 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\*' -ErrorAction SilentlyContinue | Where-Object { $_.DisplayName -like 'Inno Setup*' -and $_.InstallLocation -and $iscc.StartsWith($_.InstallLocation.TrimEnd('\'), [StringComparison]::OrdinalIgnoreCase) } | Select-Object -First 1 -ExpandProperty DisplayVersion)
$build = Join-Path $root 'build\release'
$dist = Join-Path $root 'dist\release'
$artifacts = Join-Path $root 'release_artifacts'
if (-not (Test-Path -LiteralPath $python)) { throw 'Validated .venv Python is missing.' }
if (-not (Test-Path -LiteralPath $iscc)) { throw 'Inno Setup ISCC.exe is missing.' }
if ([Environment]::Is64BitProcess -ne $true) { throw 'An x86_64 build host is required.' }
$version = (& $python -c 'from pichanalysis.version import __version__; print(__version__)').Trim()
if ($LASTEXITCODE -ne 0 -or $version -notmatch '^\d+\.\d+\.\d+$') { throw 'Invalid canonical release version.' }
foreach ($path in @($build,$dist)) {
    $resolved = [IO.Path]::GetFullPath($path)
    if (-not $resolved.StartsWith($root + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) { throw "Unsafe build path: $resolved" }
    if (Test-Path -LiteralPath $resolved) { Remove-Item -LiteralPath $resolved -Recurse -Force }
    New-Item -ItemType Directory -Path $resolved | Out-Null
}
New-Item -ItemType Directory -Path $artifacts -Force | Out-Null
if (-not $SkipGates) {
    & $python -m pytest (Join-Path $root 'tests')
    if ($LASTEXITCODE) { throw 'Python regressions failed.' }
    & $python -m compileall -q (Join-Path $root 'src')
    if ($LASTEXITCODE) { throw 'Compile check failed.' }
}
$spec = Join-Path $PSScriptRoot 'PichAnalysis.spec'
# Windows PowerShell 5.1 treats redirected native stderr as a PowerShell error.
# Exit codes below remain the authoritative gates.
$ErrorActionPreference = 'Continue'
$originalPath = $env:PATH
$env:PATH = "$(Join-Path $root '.venv\Scripts');$($env:SystemRoot)\System32;$($env:SystemRoot)"
try {
    & $python -m PyInstaller --noconfirm --distpath $dist --workpath $build $spec *> (Join-Path $artifacts 'pyinstaller-build.log')
    $pyinstallerExit = $LASTEXITCODE
} finally {
    $env:PATH = $originalPath
}
if ($pyinstallerExit) { throw 'PyInstaller build failed; inspect pyinstaller-build.log.' }
$bundle = Join-Path $dist 'PichAnalysis'
& $python (Join-Path $PSScriptRoot 'verify_release.py') $bundle (Join-Path $build 'bundle-smoke') *> (Join-Path $artifacts 'bundle-verification.log')
if ($LASTEXITCODE) { throw 'Bundle smoke failed; inspect bundle-verification.log.' }
$zip = Join-Path $artifacts "PichAnalysis-$version-Windows-x86_64-portable.zip"
if (Test-Path -LiteralPath $zip) { Remove-Item -LiteralPath $zip -Force }
Compress-Archive -LiteralPath $bundle -DestinationPath $zip -CompressionLevel Optimal
$definitionVersion = "/DAppVersion=$version"
$definitionBundle = "/DBundleDir=$bundle"
$definitionOutput = "/DArtifactDir=$artifacts"
& $iscc $definitionVersion $definitionBundle $definitionOutput (Join-Path $PSScriptRoot 'PichAnalysis.iss') *> (Join-Path $artifacts 'inno-build.log')
if ($LASTEXITCODE) { throw 'Inno Setup compilation failed; inspect inno-build.log.' }
$setup = Join-Path $artifacts "PichAnalysis-$version-Windows-x86_64-Setup.exe"
if (-not (Test-Path -LiteralPath $setup)) { throw "Installer was not produced: $setup" }
& $python (Join-Path $PSScriptRoot 'finalize_release.py') --root $root --installer $setup --portable $zip --output $artifacts --inno-version $innoVersion
if ($LASTEXITCODE) { throw 'Release manifest generation failed.' }
Write-Output "Windows release artifacts built in $artifacts"
