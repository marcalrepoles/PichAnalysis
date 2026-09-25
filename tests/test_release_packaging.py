"""Static release invariants plus manifest generation over temporary artifacts."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from pichanalysis.version import __version__, windows_file_version
from pichanalysis.core.release_preflight import CRITICAL_SCRIPTS


ROOT = Path(__file__).parents[1]


def test_single_version_source_and_metadata():
    assert __version__ == "0.1.0"
    assert windows_file_version() == "0.1.0.0"
    assert "dynamic = [\"version\"]" in (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "from .version import __version__" in (ROOT / "src/pichanalysis/__init__.py").read_text(encoding="utf-8")


def test_runtime_resource_and_bundle_rules():
    spec = (ROOT / "packaging/windows/PichAnalysis.spec").read_text(encoding="utf-8")
    resolver = (ROOT / "src/pichanalysis/core/resources.py").read_text(encoding="utf-8")
    for name in CRITICAL_SCRIPTS:
        assert (ROOT / "r_scripts" / name).is_file()
    assert 'glob("*.R")' in spec and 'r_scripts/lib' in spec
    assert 'glob("*.R")' in spec and 'r_scripts/tests' not in spec
    assert "console=False" in spec and "COLLECT(" in spec
    assert "if not getattr(sys, \"frozen\", False):" in resolver
    assert "release_artifacts/" in (ROOT / ".gitignore").read_text(encoding="utf-8")


def test_installer_is_per_user_and_preserves_external_projects():
    installer = (ROOT / "packaging/windows/PichAnalysis.iss").read_text(encoding="utf-8")
    assert "PrivilegesRequired=lowest" in installer
    assert "{localappdata}\\Programs\\PichAnalysis" in installer
    assert "{autoprograms}\\PichAnalysis" in installer
    assert "[UninstallDelete]" not in installer
    assert "[Run]" in installer


def test_release_manifest_schema_and_checksums(tmp_path):
    script = ROOT / "packaging/windows/finalize_release.py"
    specification = importlib.util.spec_from_file_location("release_finalizer", script)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    installer = tmp_path / f"PichAnalysis-{__version__}-Windows-x86_64-Setup.exe"
    portable = tmp_path / f"PichAnalysis-{__version__}-Windows-x86_64-portable.zip"
    installer.write_bytes(b"installer test fixture")
    portable.write_bytes(b"portable test fixture")
    manifest = module.finalize(ROOT, installer, portable, tmp_path, "test compiler")
    assert manifest["version"] == __version__
    assert manifest["architecture"] == "x86_64"
    assert manifest["working_tree_state"] == "modified"
    assert manifest["r_runtime_bundled"] is False
    assert manifest["scientific_databases_bundled"] is False
    assert (tmp_path / "SHA256SUMS.txt").is_file()
    assert json.loads((tmp_path / "release_manifest.json").read_text())["artifacts"] == manifest["artifacts"]

def test_branding_assets_and_minimum_splash(monkeypatch):
    from PySide6.QtCore import QEventLoop, QTimer
    from PySide6.QtWidgets import QApplication, QMainWindow
    from pichanalysis import app as app_module
    from pichanalysis.core.resources import branding_asset

    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    assert branding_asset("icon.ico").is_file()
    assert branding_asset("splash.png").is_file()
    application = QApplication.instance() or QApplication([])
    monkeypatch.setattr(app_module, "MainWindow", QMainWindow)
    window, splash = app_module.show_branded_window(application, minimum_splash_ms=3000)
    assert splash.isVisible() and not window.isVisible()
    assert not application.windowIcon().isNull()
    observations = []
    loop = QEventLoop()
    QTimer.singleShot(2500, lambda: observations.append((splash.isVisible(), window.isVisible())))
    QTimer.singleShot(3300, loop.quit)
    loop.exec()
    assert observations == [(True, False)]
    assert window.isVisible() and not splash.isVisible()
    window.close()
