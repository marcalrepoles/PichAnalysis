"""Controlled installed-bundle self-check used by release verification only."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication

from .core.consolidated_report import ReportConfig, ReportSection, SourceAdapter, generate_report
from .core.importer import import_into_project
from .core.project import create_project, open_project
from .core.r_runtime import RRuntime
from .core.release_preflight import run_preflight
from .core.resources import r_script, r_scripts_dir
from .ui.main_window import MainWindow


def run_packaged_smoke(output: Path, *, simulate_missing_r: bool = False) -> int:
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    result_path = output / ("smoke_missing_r.json" if simulate_missing_r else "smoke.json")
    result: dict[str, object] = {"status": "failed", "frozen": bool(getattr(sys, "frozen", False))}
    window = None
    try:
        app = QApplication.instance() or QApplication([])
        window = MainWindow()
        app.processEvents()
        result["window_title"] = window.windowTitle()
        result["cross_module_loaded"] = window.analyses_page.cross_module_explorer is not None
        root = output / "Release Smoke Project"
        if (root / "project.json").is_file():
            project = open_project(root)
        else:
            project = create_project(output, root.name)
            source = output / "tiny_input.csv"
            source.write_text("accession,LFQ Control 1,LFQ Treatment 1\nP04637,12,13\n", encoding="utf-8")
            import_into_project(project, source)
        project = open_project(root)
        result["project_reopened"] = project.root == root
        result["imported_rows"] = project.config["input"]["rows"]
        scripts = r_scripts_dir()
        result["r_scripts_dir"] = str(scripts)
        result["r_script_count"] = len(tuple(scripts.glob("[0-9][0-9]_*.R")))
        runtime = RRuntime(configured_executable=str(output / "missing" / "Rscript.exe")) if simulate_missing_r else RRuntime()
        if runtime.available:
            from .core.external_process import run_external
            probe = run_external((runtime.executable, "--vanilla", "-e", "cat('probe\\tTRUE\\t1.0\\n')"),
                capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=15, shell=False)
            result["r_probe"] = {"returncode": probe.returncode, "stdout": probe.stdout, "stderr": probe.stderr}
        preflight = run_preflight(project.root, runtime)
        result["preflight"] = {
            "rscript": preflight.rscript.installed,
            "resources": all(item.installed for item in preflight.resources),
            "python": all(item.installed for item in preflight.python),
            "writable": all(item.installed for item in preflight.writable),
            "required_r_packages": all(item.installed for item in preflight.r_packages
                                       if item.name not in {"org.Hs.eg.db", "org.Mm.eg.db", "MsCoreUtils", "imputeLCMD", "impute", "statmod"}),
        }
        result["missing_required_r_packages"] = [item.name for item in preflight.r_packages
            if not item.installed and item.name not in {"org.Hs.eg.db", "org.Mm.eg.db", "MsCoreUtils", "imputeLCMD", "impute", "statmod"}]
        result["r_package_failure_details"] = [item.detail for item in preflight.r_packages if item.detail][:2]
        if not simulate_missing_r:
            if not runtime.available:
                raise RuntimeError("Rscript is not available in the normal release-test environment.")
            r_output = output / "r_runtime_result.txt"
            execution = runtime.run(r_script("00_runtime_test.R"), str(r_output))
            if execution.returncode or "status=ok" not in r_output.read_text(encoding="utf-8"):
                raise RuntimeError(f"Bundled R resource execution failed: {execution.stderr}")
            result["r_resource_execution"] = True
        else:
            if preflight.rscript.installed:
                raise RuntimeError("Missing-R simulation unexpectedly detected Rscript.")

        folder = root / "analyses" / "release_fixture" / "runs" / "one"
        tables = folder / "tables"
        tables.mkdir(parents=True, exist_ok=True)
        (tables / "summary.csv").write_text("protein,count\nP04637,1\n", encoding="utf-8")
        adapter = SourceAdapter("release_fixture", "Release fixture", lambda _project: ["one"],
            lambda _project, _run_id: {"run_root": folder, "metadata": {"run_id": "one", "input_lineage": "release-smoke"}},
            lambda _project, _run_id, _outputs: folder, ("tables",))
        report = generate_report(project, ReportConfig(sections=(ReportSection(
            "release_fixture", "one", ("tables/summary.csv",)),)), (adapter,))
        result["report_html"] = (report / "report.html").is_file()
        result["report_pdf"] = (report / "report.pdf").read_bytes().startswith(b"%PDF-")
        result["status"] = "passed"
        return 0
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
        return 1
    finally:
        if window is not None:
            window.close()
        result_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
