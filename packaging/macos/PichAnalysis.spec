"""macOS windowed app specification; execute only on macOS."""

from pathlib import Path

ROOT = Path(SPECPATH).resolve().parents[1]
datas = [(str(path), "r_scripts") for path in sorted((ROOT / "r_scripts").glob("*.R"))]
datas += [(str(path), "r_scripts/lib") for path in sorted((ROOT / "r_scripts" / "lib").glob("*.R"))]
datas += [(str(ROOT / "docs" / name), "docs") for name in ("installation-windows.md", "runtime-dependencies.md")]
datas += [(str(ROOT / "THIRD_PARTY_NOTICES.txt"), ".")]

a = Analysis([str(ROOT / "packaging" / "windows" / "launcher.py")],
    pathex=[str(ROOT / "src")], binaries=[], datas=datas, hiddenimports=[],
    hookspath=[], hooksconfig={}, runtime_hooks=[], excludes=["PyQt5", "PyQt6", "pytest"], noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="PichAnalysis",
    debug=False, bootloader_ignore_signals=False, strip=False, upx=False,
    console=False, target_arch=None)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="PichAnalysis")
app = BUNDLE(coll, name="PichAnalysis.app", icon=None, bundle_identifier="org.pichanalysis.app")
