"""Regressions for reproduced application and release-readiness failures."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

from pichanalysis.core.importer import ImportError as TableImportError, read_table, import_into_project
from pichanalysis.core.project import ProjectError, create_project, open_project
from pichanalysis.core.r_runtime import discover_rscript, RRuntime
from pichanalysis.core.resources import r_script, r_scripts_dir
from pichanalysis.core.release_preflight import run_preflight, CRITICAL_SCRIPTS, Check, Preflight, R_REQUIRED_BY, R_OPTIONAL_BY


def test_resource_override_and_different_cwd(tmp_path, monkeypatch):
    root = tmp_path / "Packaged App"
    scripts = root / "r_scripts"
    scripts.mkdir(parents=True)
    (scripts / "00_runtime_test.R").write_text("cat('ok')\n")
    (scripts / "01_mapping_annotation.R").write_text("cat('ok')\n")
    monkeypatch.setenv("PICHANALYSIS_RESOURCE_ROOT", str(root))
    monkeypatch.chdir(tmp_path)
    assert r_scripts_dir() == scripts.resolve()
    assert r_script("01_mapping_annotation.R") == (scripts / "01_mapping_annotation.R").resolve()
    with pytest.raises(ValueError):
        r_script("../wrong.R")


def test_explicit_rscript_configuration_wins_over_path(tmp_path):
    configured = tmp_path / "R With Spaces" / "Rscript.exe"
    on_path = tmp_path / "other" / "Rscript.exe"
    configured.parent.mkdir(); on_path.parent.mkdir()
    configured.touch(); on_path.touch()
    assert discover_rscript(str(configured), platform="win32", which=lambda _: str(on_path)) == str(configured)


@pytest.mark.parametrize("suffix,contents", [
    ("csv", "protein,protein\nP1,P2\n"),
    ("tsv", "protein\tprotein\nP1\tP2\n"),
])
def test_duplicate_input_columns_are_rejected(tmp_path, suffix, contents):
    path = tmp_path / f"duplicate.{suffix}"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(TableImportError, match="Duplicate column"):
        read_table(path)


def test_duplicate_xlsx_columns_are_rejected(tmp_path):
    path = tmp_path / "duplicate.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame([["protein", "protein"], ["P1", "P2"]]).to_excel(
            writer, index=False, header=False)
    with pytest.raises(TableImportError, match="Duplicate column"):
        read_table(path)


@pytest.mark.parametrize("contents", ["protein,value\n", "protein,value\n,\n"])
def test_header_only_or_blank_rows_are_rejected(tmp_path, contents):
    path = tmp_path / "empty.csv"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(TableImportError, match="no usable data rows"):
        read_table(path)


def test_unicode_project_relocation_and_invalid_metadata(tmp_path):
    project = create_project(tmp_path, "Análise Proteômica With Spaces")
    assert open_project(project.root).name == project.name
    moved = tmp_path / "Moved Project"
    shutil.copytree(project.root, moved)
    assert open_project(moved).name == project.name
    (moved / "project.json").write_text("{broken", encoding="utf-8")
    with pytest.raises(ProjectError, match="project.json"):
        open_project(moved)
    (moved / "project.json").write_text(json.dumps({"schema_version": 99,
        "application": "PichAnalysis", "project_name": "x", "input": {}}), encoding="utf-8")
    with pytest.raises(ProjectError, match="compat"):
        open_project(moved)


def test_release_preflight_reports_r_packages_and_missing_resources(tmp_path, monkeypatch):
    class MissingR:
        available = False
        executable = None

    monkeypatch.setattr("pichanalysis.core.release_preflight.r_scripts_dir",
                        lambda: tmp_path)
    status = run_preflight(project_directory=tmp_path, runtime=MissingR())
    assert not status.ready
    assert len(status.resources) == len(CRITICAL_SCRIPTS)
    assert all(not item.installed for item in status.resources)
    assert all(not item.installed for item in status.r_packages)
    assert all(item.installed for item in status.writable)

def test_bom_semicolon_xlsx_sheet_and_moderate_input(tmp_path):
    project = create_project(tmp_path, "Import With Spaces")
    csv_path = tmp_path / "input semicolon.csv"
    csv_path.write_text("\ufeffprotein;value\nP1;1\n", encoding="utf-8")
    frame, _ = read_table(csv_path)
    assert frame.loc[0, "protein"] == "P1"
    workbook = tmp_path / "two sheets.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        pd.DataFrame({"protein": ["WRONG"]}).to_excel(writer, sheet_name="First", index=False)
        pd.DataFrame({"protein": ["RIGHT"]}).to_excel(writer, sheet_name="Second", index=False)
    with pytest.raises(TableImportError, match="worksheet"):
        read_table(workbook)
    selected, sheet = read_table(workbook, "Second")
    assert sheet == "Second" and selected.loc[0, "protein"] == "RIGHT"
    large = tmp_path / "moderate.csv"
    large.write_text("protein,value\n" + "".join(f"P{i},{i}\n" for i in range(2000)),
                     encoding="utf-8")
    imported = import_into_project(project, large)
    assert len(imported.dataframe) == 2000
    assert imported.processed_path.is_file()
    assert imported.original_path.read_bytes() == large.read_bytes()

def test_r_process_failure_keeps_context_and_output(tmp_path):
    script = tmp_path / "script with spaces.py"
    script.write_text("import sys\nprint('context from stdout')\n"
                      "print('failure from stderr', file=sys.stderr)\nsys.exit(7)\n",
                      encoding="utf-8")
    runtime = RRuntime(sys.executable)
    result = runtime.run(script)
    assert result.returncode == 7
    assert "context from stdout" in result.stdout
    assert "failure from stderr" in result.stderr
    assert str(script) in result.command

def test_preflight_does_not_require_both_species_orgdb_packages():
    status = Preflight(
        python=(Check("PySide6", "Python", True),),
        rscript=Check("Rscript", "R", True),
        r_packages=tuple(Check(name, "R", True) for name in R_REQUIRED_BY)
            + tuple(Check(name, "conditional", False) for name in R_OPTIONAL_BY),
        resources=(Check("00_runtime_test.R", "resources", True),),
        writable=(Check("Temporary directory", "output", True),),
    )
    assert status.ready