"""Fresh-process offline application lifecycle smoke; no scientific run is started."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

from PySide6.QtWidgets import QApplication

from pichanalysis.core.importer import import_into_project
from pichanalysis.core.project import ProjectError, create_project, open_project
from pichanalysis.core.resources import r_script
from pichanalysis.ui.main_window import MainWindow


def main() -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    application = QApplication.instance() or QApplication([])
    with tempfile.TemporaryDirectory(prefix="pich_hardening_") as temporary:
        parent = Path(temporary)
        project = create_project(parent, "Análise Proteômica With Spaces")
        window = MainWindow()
        window._activate_project(project, "Smoke project created")
        for index in range(window.navigation.count()):
            window.navigation.setCurrentRow(index)
            assert window.pages.currentIndex() == index
        assert window.report_page.project.root == project.root
        input_file = parent / "input with spaces.csv"
        input_file.write_text("protein,LFQ A,LFQ B\nP12345,4,5\nQ12345,6,7\n", encoding="utf-8")
        imported = import_into_project(project, input_file)
        assert imported.processed_path.is_file()
        window._activate_project(open_project(project.root), "Smoke project reopened")
        assert window.data_page is not None
        assert window.report_page.history.rowCount() == 0
        assert r_script("00_runtime_test.R").is_file()
        window.close()
        assert not window.isVisible()
        imported.processed_path.unlink()
        incomplete_window = MainWindow()
        incomplete_window._activate_project(open_project(project.root), "Missing processed input")
        assert "Imported project artifact is missing" in incomplete_window.data_page.summary.text()
        incomplete_window.close()
        from pichanalysis.core import r_runtime
        original_discovery = r_runtime.discover_rscript
        try:
            r_runtime.discover_rscript = lambda *args, **kwargs: None
            missing_r_window = MainWindow()
            assert not missing_r_window.runtime.available
            assert "Not found" in missing_r_window.scripts_page.runtime_status.text()
            missing_r_window.close()
        finally:
            r_runtime.discover_rscript = original_discovery

        invalid = parent / "invalid"
        invalid.mkdir()
        try:
            open_project(invalid)
        except ProjectError:
            pass
        else:
            raise AssertionError("Missing project.json was accepted")
        (invalid / "project.json").write_text("{invalid", encoding="utf-8")
        try:
            open_project(invalid)
        except ProjectError:
            pass
        else:
            raise AssertionError("Corrupted project.json was accepted")
    print("Application hardening offline smoke: PASS")


if __name__ == "__main__":
    main()
