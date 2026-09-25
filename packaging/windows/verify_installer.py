"""Controlled per-user install, launch, uninstall and reinstall test."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
import winreg
from pathlib import Path

from verify_release import verify


def registered_installations() -> list[dict[str, str]]:
    found = []
    key = winreg.OpenKey(winreg.HKEY_CURRENT_USER,
        r"Software\Microsoft\Windows\CurrentVersion\Uninstall")
    with key:
        for index in range(winreg.QueryInfoKey(key)[0]):
            name = winreg.EnumKey(key, index)
            with winreg.OpenKey(key, name) as entry:
                try:
                    display = winreg.QueryValueEx(entry, "DisplayName")[0]
                except FileNotFoundError:
                    continue
                if display == "PichAnalysis":
                    try:
                        location = winreg.QueryValueEx(entry, "InstallLocation")[0]
                    except FileNotFoundError:
                        location = ""
                    found.append({"key": name, "location": location})
    return found


def run_checked(command: list[str], log: Path) -> None:
    result = subprocess.run(command + [f"/LOG={log}"], timeout=240, check=False)
    if result.returncode:
        raise RuntimeError(f"Installer lifecycle command exited {result.returncode}; see {log}")


def verify_installer(installer: Path, output: Path) -> dict:
    installer, output = installer.resolve(), output.resolve()
    if not installer.is_file():
        raise FileNotFoundError(installer)
    existing = registered_installations()
    if existing:
        raise RuntimeError(f"Refusing to overwrite an existing PichAnalysis installation: {existing}")
    output.mkdir(parents=True, exist_ok=True)
    install_dir = output / "Installed App With Spaces"
    project_dir = output / "project-check"
    shortcut = Path(os.environ["APPDATA"]) / "Microsoft/Windows/Start Menu/Programs/PichAnalysis.lnk"
    base = [str(installer), "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART", f"/DIR={install_dir}"]
    history = []
    try:
        run_checked(base, output / "install-first.log")
        if not (install_dir / "PichAnalysis.exe").is_file() or not shortcut.is_file():
            raise AssertionError("Installed executable or Start Menu shortcut is missing")
        history.append({"phase": "fresh_install", "smoke": verify(install_dir, project_dir)})
        project = project_dir / "Release Smoke Project" / "project.json"
        if not project.is_file():
            raise AssertionError("External test project was not created")
        uninstaller = install_dir / "unins000.exe"
        if not uninstaller.is_file():
            raise AssertionError("Inno uninstaller is missing")
        run_checked([str(uninstaller), "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART"], output / "uninstall-first.log")
        if (install_dir / "PichAnalysis.exe").exists() or not project.is_file():
            raise AssertionError("Uninstall failed or removed the external user project")
        history.append({"phase": "uninstall", "project_preserved": True})
        run_checked(base, output / "install-second.log")
        history.append({"phase": "reinstall", "smoke": verify(install_dir, project_dir)})
        run_checked(base, output / "install-same-version.log")
        if not shortcut.is_file():
            raise AssertionError("Start Menu shortcut missing after same-version reinstall")
        history.append({"phase": "same_version_reinstall", "smoke": verify(install_dir, project_dir)})
        uninstallers = sorted(install_dir.glob("unins*.exe"))
        if not uninstallers:
            raise AssertionError("No Inno uninstaller exists after same-version reinstall")
        run_checked([str(uninstallers[-1]), "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART"],
            output / "uninstall-final.log")
        if (install_dir / "PichAnalysis.exe").exists() or not project.is_file():
            raise AssertionError("Final uninstall failed or removed the external user project")
        history.append({"phase": "final_cleanup", "project_preserved": True})
        return {"status": "passed", "install_directory": str(install_dir),
            "external_project": str(project.parent), "history": history}
    finally:
        (output / "installer-verification.json").write_text(json.dumps({
            "completed_phases": [item["phase"] for item in history],
            "external_project_exists": (project_dir / "Release Smoke Project" / "project.json").is_file(),
            "installed_executable_exists": (install_dir / "PichAnalysis.exe").is_file(),
        }, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("installer", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify_installer(args.installer, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
