from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox, QLabel, QPlainTextEdit, QPushButton, QVBoxLayout, QWidget,
)

from ..core.r_runtime import RRuntime


class ScriptsPage(QWidget):
    test_requested = Signal()

    def __init__(self, scripts_dir: Path, runtime: RRuntime) -> None:
        super().__init__()
        self.scripts_dir = scripts_dir
        self.project_scripts_dir: Path | None = None
        self.runtime_status = QLabel()
        self.test_button = QPushButton("Test R")
        self.script_picker = QComboBox()
        self.script_view = QPlainTextEdit()
        self.script_view.setReadOnly(True)
        self.result_view = QPlainTextEdit()
        self.result_view.setReadOnly(True)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("R Runtime"))
        layout.addWidget(self.runtime_status)
        layout.addWidget(self.test_button)
        layout.addWidget(QLabel("Application R scripts"))
        layout.addWidget(self.script_picker)
        layout.addWidget(self.script_view, 2)
        layout.addWidget(QLabel("Result / Session logs"))
        layout.addWidget(self.result_view, 1)
        self.test_button.clicked.connect(self.test_requested)
        self.script_picker.currentTextChanged.connect(self._show_script)
        self.refresh_runtime(runtime)
        self.refresh_scripts()

    def refresh_runtime(self, runtime: RRuntime) -> None:
        if runtime.available:
            self.runtime_status.setText(
                f"Found\nExecutable: {runtime.executable}\nVersion: {runtime.version() or 'unknown'}"
            )
        else:
            self.runtime_status.setText(
                "Not found. Scientific analyses require R and Rscript on PATH."
            )
        self.test_button.setEnabled(runtime.available)

    def refresh_scripts(self) -> None:
        self.script_picker.clear()
        paths = list(self.scripts_dir.rglob("*.R")) if self.scripts_dir.is_dir() else []
        if self.project_scripts_dir and self.project_scripts_dir.is_dir():
            paths.extend(self.project_scripts_dir.rglob("*.R"))
        for path in sorted(paths):
            try:
                label = str(path.relative_to(self.scripts_dir))
            except ValueError:
                label = "snapshot/" + str(path.relative_to(self.project_scripts_dir))
            self.script_picker.addItem(label, str(path))

    def set_project_scripts(self, directory: Path | None) -> None:
        self.project_scripts_dir = directory
        self.refresh_scripts()

    def _show_script(self, name: str) -> None:
        if not name:
            self.script_view.clear()
            return
        try:
            path = Path(self.script_picker.currentData())
            self.script_view.setPlainText(path.read_text(encoding="utf-8"))
        except OSError as error:
            self.script_view.setPlainText(f"Could not open the script: {error}")
