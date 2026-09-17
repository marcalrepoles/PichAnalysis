from __future__ import annotations

import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QUrl
from PySide6.QtGui import QCloseEvent, QDesktopServices
from PySide6.QtWidgets import (
    QFileDialog, QInputDialog, QListWidget, QMainWindow, QMessageBox,
    QStackedWidget, QHBoxLayout, QWidget,
)

from ..core.column_mapping import load_mapping, save_mapping, validate_mapping
from ..core.importer import ImportError, ImportResult, import_into_project, xlsx_sheets
from ..core.mapping_analysis import (
    build_mapping_arguments, cache_path, export_result, new_run_id, read_mapping_outputs,
)
from ..core.organism import get_organism, has_biological_results, set_organism
from ..core.go_analysis import export_go, prepare_go_arguments, read_go_outputs
from ..core.presence_analysis import (
    build_presence_arguments, export_presence, read_presence_outputs,
)
from ..core.project import Project, ProjectError, create_project, open_project
from ..core.r_runtime import RRuntime
from .analyses_page import AnalysesPage
from .data_page import DataPage
from .database_manager_page import DatabaseManagerPage
from .mapping_worker import MappingWorker
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
        self.mapping_worker: MappingWorker | None = None
        self.presence_worker: MappingWorker | None = None
        self.go_worker: MappingWorker | None = None
        self.project_page = ProjectPage()
        self.data_page = DataPage()
        self.database_manager_page = DatabaseManagerPage()
        self.analyses_page = AnalysesPage()
        self.scripts_page = ScriptsPage(APPLICATION_ROOT / "r_scripts", self.runtime)
        self.navigation = QListWidget()
        self.navigation.addItems(["Project", "Data", "Analyses", "Database Manager", "Scripts / Logs"])
        self.navigation.setFixedWidth(170)
        self.pages = QStackedWidget()
        for page in (self.project_page, self.data_page, self.analyses_page, self.database_manager_page, self.scripts_page):
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
        self.data_page.save_mapping_requested.connect(self._save_mapping)
        self.scripts_page.test_requested.connect(self._test_r)
        self.analyses_page.organism_requested.connect(self._save_organism)
        self.analyses_page.run_requested.connect(self._run_mapping)
        self.analyses_page.export_table_requested.connect(lambda: self._export_mapping("protein_catalog.csv"))
        self.analyses_page.export_workbook_requested.connect(lambda: self._export_mapping("protein_mapping.xlsx"))
        self.analyses_page.open_results_requested.connect(self._open_results)
        presence = self.analyses_page.presence_page
        presence.run_requested.connect(self._run_presence)
        presence.export_table_requested.connect(lambda: self._export_presence("tables/classification.csv"))
        presence.export_workbook_requested.connect(lambda: self._export_presence("tables/presence_absence.xlsx"))
        presence.export_graph_requested.connect(self._export_presence_graph)
        presence.open_graphs_requested.connect(self._open_presence_graphs)
        go = self.analyses_page.go_page
        go.run_requested.connect(self._run_go)
        go.export_table_requested.connect(self._export_go_path)
        go.export_workbook_requested.connect(lambda:self._export_go_path(str(self.project.root/"analyses"/"GO"/"GO_analysis.xlsx") if self.project else ""))
        go.export_graph_requested.connect(self._export_go_path)
        go.open_folder_requested.connect(self._open_go)
        self._set_project_enabled(False)

    def _set_project_enabled(self, enabled: bool) -> None:
        self.data_page.import_button.setEnabled(enabled)
        self.navigation.item(2).setFlags(
            self.navigation.item(2).flags() | self.navigation.item(2).flags().ItemIsEnabled
            if enabled else self.navigation.item(2).flags() & ~self.navigation.item(2).flags().ItemIsEnabled
        )

    def closeEvent(self, event: QCloseEvent) -> None:
        if ((self.mapping_worker and self.mapping_worker.isRunning()) or
                (self.presence_worker and self.presence_worker.isRunning()) or
                (self.go_worker and self.go_worker.isRunning()) or self.database_manager_page.is_running()):
            QMessageBox.information(
                self, "Operation in progress",
                "Wait for the current analysis or database download before closing PichAnalysis.",
            )
            event.ignore()
            return
        event.accept()

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
        self._restore_data()
        self.analyses_page.set_project(project)
        self.scripts_page.set_project_scripts(project.root / "scripts" / "runs")
        try:
            if (project.root / "mapping" / "latest_metadata.json").is_file():
                self.analyses_page.show_outputs(read_mapping_outputs(project))
        except RuntimeError:
            self.logger.exception("Não foi possível restaurar resultados de mapeamento")
        try:
            if (project.root / "analyses" / "presence_absence" / "latest_metadata.json").is_file():
                self.analyses_page.presence_page.show_outputs(read_presence_outputs(project))
        except RuntimeError:
            self.logger.exception("Não foi possível restaurar presença/ausência")
        try:
            if (project.root/"analyses"/"GO"/"latest_metadata.json").is_file(): self.analyses_page.go_page.show_outputs(read_go_outputs(project))
        except RuntimeError: self.logger.exception("Não foi possível restaurar GO")

    def _restore_data(self) -> None:
        if not self.project:
            return
        input_config = self.project.config.get("input", {})
        processed_name = input_config.get("processed_file")
        original_name = input_config.get("original_file")
        if not processed_name or not original_name:
            return
        try:
            processed = self.project.root / processed_name
            original = self.project.root / original_name
            frame = pd.read_csv(processed, encoding="utf-8")
            result = ImportResult(
                frame,
                original,
                processed,
                str(input_config.get("format", "csv")),
                input_config.get("sheet"),
            )
            mapping = load_mapping(self.project)
            self.data_page.show_result(result, mapping)
            if mapping:
                self.data_page.show_validation(validate_mapping(mapping))
        except Exception:
            self.logger.exception("Não foi possível restaurar os dados processados")

    def _new_project(self) -> None:
        name, accepted = QInputDialog.getText(self, "Novo projeto", "Nome do projeto:")
        if not accepted or not name.strip():
            return
        parent = QFileDialog.getExistingDirectory(self, "Choose the folder where the project will be created")
        if not parent:
            return
        try:
            self._activate_project(create_project(Path(parent), name), "Project created")
        except ProjectError as error:
            QMessageBox.warning(self, "Could not create project", str(error))

    def _open_project(self) -> None:
        directory = QFileDialog.getExistingDirectory(self, "Open PichAnalysis project")
        if not directory:
            return
        try:
            self._activate_project(open_project(Path(directory)), "Project opened")
        except ProjectError as error:
            QMessageBox.warning(self, "Invalid project", str(error))

    def _open_folder(self) -> None:
        if self.project and not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project.root))):
            QMessageBox.warning(self, "Folder", "Could not open the project folder.")

    def _import_data(self) -> None:
        if not self.project:
            return
        filename, _ = QFileDialog.getOpenFileName(
            self, "Import table", "", "Tables (*.xlsx *.csv *.tsv)"
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
                        self, "Choose worksheet", "Worksheet:", sheets, 0, False
                    )
                    if not accepted:
                        return
            result = import_into_project(self.project, source, sheet)
            self.data_page.show_result(result, load_mapping(self.project))
            self.project_page.show_project(self.project)
            self.logger.info("Importação concluída: %s", result.original_path.name)
            self.navigation.setCurrentRow(1)
        except (ImportError, ProjectError) as error:
            self.logger.exception("Erro de importação")
            QMessageBox.warning(self, "Import failed", str(error))

    def _save_mapping(self, columns: dict) -> None:
        if not self.project:
            return
        try:
            validation = save_mapping(self.project, columns)
            self.data_page.show_validation(validation)
            self.logger.info(
                "Configuração de colunas salva: %s",
                "válida" if validation.valid else "inválida",
            )
            self.analyses_page.set_project(self.project)
        except ProjectError as error:
            self.logger.exception("Erro ao salvar configuração de colunas")
            QMessageBox.warning(self, "Configuration", str(error))

    def _save_organism(self, name: str, tax_id: str) -> None:
        if not self.project:
            return
        current = get_organism(self.project)
        changing = current and (current.name != name.strip() or current.tax_id != tax_id.strip())
        if changing and has_biological_results(self.project):
            answer = QMessageBox.question(self, "Change organism",
                "Biological results already exist. Changing the organism invalidates them. Continue?")
            if answer != QMessageBox.StandardButton.Yes:
                return
        try:
            set_organism(self.project, name, tax_id)
            self.logger.info("Organismo configurado: %s (%s)", name.strip(), tax_id.strip())
            self.analyses_page.set_project(self.project)
        except (ValueError, ProjectError) as error:
            QMessageBox.warning(self, "Organism", str(error))

    def _run_mapping(self, refresh: bool) -> None:
        if not self.project or self.mapping_worker:
            return
        run_id = new_run_id()
        try:
            arguments = build_mapping_arguments(self.project, cache_path(), refresh, run_id)
        except ValueError as error:
            QMessageBox.warning(self, "Mapping", str(error))
            return
        self.logger.info("Mapeamento iniciado run_id=%s organismo=%s tipo=%s cache=%s",
            run_id, self.project.config.get("organism_tax_id"),
            arguments[arguments.index("--id-type") + 1], not refresh)
        worker = MappingWorker(self.runtime, APPLICATION_ROOT / "r_scripts" / "01_mapping_annotation.R", arguments)
        self.mapping_worker = worker
        self.analyses_page.set_running(True)
        worker.succeeded.connect(lambda stdout, stderr: self._mapping_finished(run_id, stdout, stderr))
        worker.failed.connect(lambda message, stderr: self._mapping_failed(run_id, message, stderr))
        worker.finished.connect(worker.deleteLater)
        worker.start()

    def _mapping_finished(self, run_id: str, stdout: str, stderr: str) -> None:
        self.mapping_worker = None
        self.analyses_page.set_running(False)
        try:
            outputs = read_mapping_outputs(self.project)
            self.analyses_page.show_outputs(outputs)
            self.scripts_page.set_project_scripts(self.project.root / "scripts" / "runs")
            meta = outputs.metadata
            self.logger.info("Mapeamento finalizado run_id=%s mapped=%s ambiguous=%s unmapped=%s exit_code=0",
                run_id, meta.get("mapped_unique_count"), meta.get("ambiguous_count"), meta.get("unmapped_count"))
        except RuntimeError as error:
            QMessageBox.warning(self, "Resultados", str(error))

    def _mapping_failed(self, run_id: str, message: str, stderr: str) -> None:
        self.mapping_worker = None
        self.analyses_page.set_running(False)
        self.logger.error("Mapeamento falhou run_id=%s: %s", run_id, message)
        self.scripts_page.result_view.setPlainText(stderr)
        QMessageBox.warning(self, "Mapping failed", message)

    def _export_mapping(self, filename: str) -> None:
        if not self.project:
            return
        source = self.project.root / "mapping" / "tables" / filename
        destination, _ = QFileDialog.getSaveFileName(self, "Export result", filename)
        if destination:
            try:
                export_result(source, Path(destination))
            except OSError as error:
                QMessageBox.warning(self, "Export", str(error))

    def _open_results(self) -> None:
        if self.project:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project.root / "mapping" / "tables")))

    def _run_presence(self, parameters: dict) -> None:
        if not self.project or self.presence_worker:
            return
        run_id = new_run_id()
        try:
            arguments = build_presence_arguments(self.project, run_id=run_id, **parameters)
        except ValueError as error:
            QMessageBox.warning(self, "Presence / absence", str(error)); return
        self.logger.info("Presença/ausência iniciada run_id=%s tipo=%s condições=%s threshold=%s zero_ausência=%s",
            run_id, parameters["quantification_type"], len(parameters["selected_conditions"]),
            parameters["threshold"], parameters["zero_is_missing"])
        worker = MappingWorker(self.runtime, APPLICATION_ROOT / "r_scripts" / "02_presence_absence.R", arguments)
        self.presence_worker = worker; self.analyses_page.presence_page.set_running(True)
        worker.succeeded.connect(lambda stdout,stderr:self._presence_finished(run_id,stdout,stderr))
        worker.failed.connect(lambda message,stderr:self._presence_failed(run_id,message,stderr))
        worker.finished.connect(worker.deleteLater); worker.start()

    def _presence_finished(self, run_id: str, stdout: str, stderr: str) -> None:
        self.presence_worker=None; self.analyses_page.presence_page.set_running(False)
        try:
            outputs=read_presence_outputs(self.project); self.analyses_page.presence_page.show_outputs(outputs)
            self.scripts_page.set_project_scripts(self.project.root / "scripts" / "runs")
            self.logger.info("Presença/ausência finalizada run_id=%s entidades=%s exit_code=0",run_id,outputs.metadata.get("total_entities"))
        except RuntimeError as error: QMessageBox.warning(self,"Resultados",str(error))

    def _presence_failed(self, run_id: str, message: str, stderr: str) -> None:
        self.presence_worker=None; self.analyses_page.presence_page.set_running(False)
        self.logger.error("Presença/ausência falhou run_id=%s: %s",run_id,message)
        self.scripts_page.result_view.setPlainText(stderr); QMessageBox.warning(self,"Analysis failed",message)

    def _export_presence(self, relative: str) -> None:
        if not self.project:return
        source=self.project.root/"analyses"/"presence_absence"/relative
        destination,_=QFileDialog.getSaveFileName(self,"Export result",source.name)
        if destination:
            try: export_presence(source,Path(destination))
            except OSError as error: QMessageBox.warning(self,"Export",str(error))

    def _export_presence_graph(self, source: str) -> None:
        if not source:return
        destination,_=QFileDialog.getSaveFileName(self,"Export graph",Path(source).name)
        if destination:
            try: export_presence(Path(source),Path(destination))
            except OSError as error: QMessageBox.warning(self,"Export",str(error))

    def _open_presence_graphs(self) -> None:
        if self.project: QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project.root/"analyses"/"presence_absence"/"graphs")))

    def _run_go(self, parameters: dict) -> None:
        if not self.project or self.go_worker:return
        run_id=new_run_id()
        try: arguments=prepare_go_arguments(self.project,run_id=run_id,**parameters)
        except ValueError as error:
            if "do not belong to the selected background" in str(error):
                answer=QMessageBox.question(self,"Target outside background",str(error)+"\nRestrict the target to the background and continue?")
                if answer!=QMessageBox.StandardButton.Yes:return
                try:arguments=prepare_go_arguments(self.project,run_id=run_id,allow_target_outside_background=True,**parameters)
                except ValueError as second:QMessageBox.warning(self,"Gene Ontology",str(second));return
            else:QMessageBox.warning(self,"Gene Ontology",str(error));return
        self.logger.info("GO iniciado run_id=%s target=%s background=%s ontologias=%s",run_id,parameters["target_selection"],parameters["background_selection"],",".join(parameters["ontologies"]))
        worker=MappingWorker(self.runtime,APPLICATION_ROOT/"r_scripts"/"03_go_analysis.R",arguments);self.go_worker=worker;self.analyses_page.go_page.set_running(True)
        worker.succeeded.connect(lambda stdout,stderr:self._go_finished(run_id,stdout,stderr));worker.failed.connect(lambda message,stderr:self._go_failed(run_id,message,stderr));worker.finished.connect(worker.deleteLater);worker.start()

    def _go_finished(self,run_id:str,stdout:str,stderr:str)->None:
        self.go_worker=None;self.analyses_page.go_page.set_running(False)
        try:
            outputs=read_go_outputs(self.project);self.analyses_page.go_page.show_outputs(outputs);self.scripts_page.set_project_scripts(self.project.root/"scripts"/"runs")
            self.logger.info("GO finalizado run_id=%s anotadas=%s não_anotadas=%s exit_code=0",run_id,outputs.metadata.get("annotated_count"),outputs.metadata.get("unannotated_count"))
        except RuntimeError as error:QMessageBox.warning(self,"Resultados GO",str(error))

    def _go_failed(self,run_id:str,message:str,stderr:str)->None:
        self.go_worker=None;self.analyses_page.go_page.set_running(False);self.logger.error("GO failed run_id=%s: %s",run_id,message);self.scripts_page.result_view.setPlainText(stderr);QMessageBox.warning(self,"GO failed",message)

    def _export_go_path(self,source:str)->None:
        if not source:return
        destination,_=QFileDialog.getSaveFileName(self,"Export GO result",Path(source).name)
        if destination:
            try:export_go(Path(source),Path(destination))
            except OSError as error:QMessageBox.warning(self,"Export",str(error))

    def _open_go(self)->None:
        if self.project:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project.root/"analyses"/"GO")))

    def _test_r(self) -> None:
        script = APPLICATION_ROOT / "r_scripts" / "00_runtime_test.R"
        output_dir = self.project.root / "logs" if self.project else Path(tempfile.gettempdir())
        output = output_dir / "r_runtime_test.txt"
        try:
            result = self.runtime.run(script, str(output))
            detail = output.read_text(encoding="utf-8") if output.exists() else result.stdout
            if result.returncode != 0:
                raise RuntimeError(result.stderr or "Rscript exited with an error.")
            self.scripts_page.result_view.setPlainText(detail)
            self.logger.info("Teste do R concluído")
        except (RuntimeError, OSError) as error:
            self.logger.exception("Teste do R falhou")
            QMessageBox.warning(self, "Teste do R", str(error))
