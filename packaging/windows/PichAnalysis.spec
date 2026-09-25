"""Windows x86_64 onedir bundle from the repository's validated environment."""

from pathlib import Path

from PyInstaller.utils.win32.versioninfo import FixedFileInfo, StringFileInfo, StringStruct, StringTable, VarFileInfo, VarStruct, VSVersionInfo


ROOT = Path(SPECPATH).resolve().parents[1]
version_namespace = {}
exec((ROOT / "src" / "pichanalysis" / "version.py").read_text(encoding="utf-8"), version_namespace)
VERSION = version_namespace["__version__"]
parts = tuple(int(part) for part in VERSION.split(".")) + (0,)

# Only application-owned runtime R scripts are distributed; never tests or snapshots.
datas = [(str(path), "r_scripts") for path in sorted((ROOT / "r_scripts").glob("*.R"))]
datas += [(str(path), "r_scripts/lib") for path in sorted((ROOT / "r_scripts" / "lib").glob("*.R"))]
datas += [(str(ROOT / "docs" / name), "docs") for name in ("installation-windows.md", "runtime-dependencies.md")]
datas += [(str(ROOT / "docs" / "releases" / "0.1.0.md"), "docs/releases")]
datas += [(str(ROOT / "THIRD_PARTY_NOTICES.txt"), ".")]

version_info = VSVersionInfo(
    ffi=FixedFileInfo(filevers=parts, prodvers=parts, mask=0x3F, flags=0, OS=0x40004, fileType=1, subtype=0),
    kids=[
        StringFileInfo([StringTable("040904B0", [
            StringStruct("FileDescription", "PichAnalysis"),
            StringStruct("FileVersion", ".".join(map(str, parts))),
            StringStruct("ProductName", "PichAnalysis"),
            StringStruct("ProductVersion", VERSION),
        ])]),
        VarFileInfo([VarStruct("Translation", [1033, 1200])]),
    ],
)

a = Analysis(
    [str(ROOT / "packaging" / "windows" / "launcher.py")],
    pathex=[str(ROOT / "src")],
    binaries=[],
    datas=datas,
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["PyQt5", "PyQt6", "pytest"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True, name="PichAnalysis",
    debug=False, bootloader_ignore_signals=False, strip=False, upx=False,
    console=False, version=version_info,
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="PichAnalysis")
