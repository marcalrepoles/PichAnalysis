from __future__ import annotations

import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QFileDialog, QInputDialog, QListWidget, QMainWindow, QMessageBox,
    QStackedWidget, QHBoxLayout, QWidget,
)

from ..core.importer import ImportError, import_into_project, read_table, xlsx_sheets
from ..core.project import Project, ProjectError, create_project, open_project
from ..core.r_runtime import RRuntime
from .analyses_page import AnalysesPage
from .data_page import DataPage
from .project_page import ProjectPage
from .scripts_page import ScriptsPage


APPLICATION_ROOT = Path(__file__).resolve().parents[3]


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("PichAnalysis")
        self.resize(1100, 720)
        self.project: Project | None = None
        self.runtime = RRuntime()
        self.logger = logging.getLogger("pichanalysis")
        self.logger.setLevel(logging.INFO)
        self.project_page = ProjectPage()
        self.data_page = DataPage()
        self.analyses_page = AnalysesPage()
        self.scripts_page = ScriptsPage(APPLICATION_ROOT / "r_scripts", self.runtime)
        self.navigation = QListWidget()
        self.navigation.addItems(["Projeto", "Dados", "Análises", "Scripts / Logs"])
        self.navigation.setFixedWidth(170)
        self.pages = QStackedWidget()
        for page in (self.project_page, self.data_page, self.analyses_page, self.scripts_page):
            self.pages.addWidget(page)
        container = QWidget()
        layout = QHBoxLayout(container)
        layout.addWidget(self.navigation)
        layout.addWidget(self.pages, 1)
        self.setCentralWidget(container)
        self.navigation.currentRowChanged.connect(self.pages.setCurrentIndex)
        self.navigation.setCurrentRow(0)
        self.project_page.new_button.clicked.connect(self._new_project)
        self.project_page.open_button.clicked.connect(self._open_project)
        self.project_page.folder_button.clicked.connect(self._open_folder)
        self.data_page.import_button.clicked.connect(self._import_data)
        self.scripts_page.test_requested.connect(self._test_r)
        self._set_project_enabled(False)

    def _set_project_enabled(self, enabled: bool) -> None:
        self.data_page.import_button.setEnabled(enabled)
        self.navigation.item(2).setFlags(
            self.navigation.item(2).flags() | self.navigation.item(2).flags().ItemIsEnabled
            if enabled else self.navigation.item(2).flags() & ~self.navigation.item(2).flags().ItemIsEnabled
        )

    def _configure_log(self) -> None:
        for handler in tuple(self.logger.handlers):
            handler.close()
            self.logger.removeHandler(handler)
        if self.project:
            handler = logging.FileHandler(self.project.root / "logs" / "pichanalysis.log", encoding="utf-8")
            handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
            self.logger.addHandler(handler)

    def _activate_project(self, project: Project, event: str) -> None:
        self.project = project
        self._configure_log()
        self.logger.info(event)
        self.project_page.show_project(project)
        self._set_project_enabled(True)

    def _new_project(self) -> None:
        name, accepted = QInputDialog.getText(self, "Novo projeto", "Nome do projeto:")
        if not accepted or not name.strip():
            return
        parent = QFileDialog.getExistingDirectory(self, "Escolha a pasta onde criar o projeto")
        if not parent:
            return
        try:
            self._activate_project(create_project(Path(parent), name), "Projeto criado")
        except ProjectError as error:
            QMessageBox.warning(self, "Não foi possível criar", str(error))

    def _open_project(self) -> None:
        directory = QFileDialog.getExistingDirectory(self, "Abrir projeto PichAnalysis")
        if not directory:
            return
        try:
            self._activate_project(open_project(Path(directory)), "Projeto aberto")
        except ProjectError as error:
            QMessageBox.warning(self, "Projeto inválido", str(error))

    def _open_folder(self) -> None:
        if self.project and not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project.root))):
            QMessageBox.warning(self, "Pasta", "Não foi possível abrir a pasta do projeto.")

    def _import_data(self) -> None:
        if not self.project:
            return
        filename, _ = QFileDialog.getOpenFileName(
            self, "Importar tabela", "", "Tabelas (*.xlsx *.csv *.tsv)"
        )
        if not filename:
            return
        source = Path(filename)
        sheet = None
        try:
            if source.suffix.lower() == ".xlsx":
                sheets = xlsx_sheets(source)
                if len(sheets) > 1:
                    sheet, accepted = QInputDialog.getItem(
                        self, "Escolher planilha", "Planilha:", sheets, 0, False
                    )
                    if not accepted:
                        return
            result = import_into_project(self.project, source, sheet)
            self.data_page.show_result(result)
            self.project_page.show_project(self.project)
            self.logger.info("Importação concluída: %s", result.original_path.name)
            self.navigation.setCurrentRow(1)
        except (ImportError, ProjectError) as error:
            self.logger.exception("Erro de importação")
            QMessageBox.warning(self, "Falha na importação", str(error))

    def _test_r(self) -> None:
        script = APPLICATION_ROOT / "r_scripts" / "00_runtime_test.R"
        output_dir = self.project.root / "logs" if self.project else Path(tempfile.gettempdir())
        output = output_dir / "r_runtime_test.txt"
        try:
            result = self.runtime.run(script, str(output))
            detail = output.read_text(encoding="utf-8") if output.exists() else result.stdout
            if result.returncode != 0:
                raise RuntimeError(result.stderr or "Rscript terminou com erro.")
            self.scripts_page.result_view.setPlainText(detail)
            self.logger.info("Teste do R concluído")
        except (RuntimeError, OSError) as error:
            self.logger.exception("Teste do R falhou")
            QMessageBox.warning(self, "Teste do R", str(error))

