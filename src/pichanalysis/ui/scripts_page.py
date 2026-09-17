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
        self.runtime_status = QLabel()
        self.test_button = QPushButton("Testar R")
        self.script_picker = QComboBox()
        self.script_view = QPlainTextEdit()
        self.script_view.setReadOnly(True)
        self.result_view = QPlainTextEdit()
        self.result_view.setReadOnly(True)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("R Runtime"))
        layout.addWidget(self.runtime_status)
        layout.addWidget(self.test_button)
        layout.addWidget(QLabel("Scripts R do aplicativo"))
        layout.addWidget(self.script_picker)
        layout.addWidget(self.script_view, 2)
        layout.addWidget(QLabel("Resultado / Logs da sessão"))
        layout.addWidget(self.result_view, 1)
        self.test_button.clicked.connect(self.test_requested)
        self.script_picker.currentTextChanged.connect(self._show_script)
        self.refresh_runtime(runtime)
        self.refresh_scripts()

    def refresh_runtime(self, runtime: RRuntime) -> None:
        if runtime.available:
            self.runtime_status.setText(
                f"Encontrado\nExecutável: {runtime.executable}\nVersão: {runtime.version() or 'não identificada'}"
            )
        else:
            self.runtime_status.setText(
                "Não encontrado. As análises científicas exigirão R e Rscript no PATH."
            )
        self.test_button.setEnabled(runtime.available)

    def refresh_scripts(self) -> None:
        self.script_picker.clear()
        if self.scripts_dir.is_dir():
            self.script_picker.addItems([path.name for path in sorted(self.scripts_dir.glob("*.R"))])

    def _show_script(self, name: str) -> None:
        if not name:
            self.script_view.clear()
            return
        try:
            self.script_view.setPlainText((self.scripts_dir / name).read_text(encoding="utf-8"))
        except OSError as error:
            self.script_view.setPlainText(f"Não foi possível abrir o script: {error}")

