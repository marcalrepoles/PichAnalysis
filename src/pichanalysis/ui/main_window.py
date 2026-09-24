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
from ..core.kegg_analysis import prepare_kegg_arguments,read_kegg_outputs
from ..core.reactome_analysis import ReactomeParameters
from ..core.mitocarta_analysis import MitoCartaParameters
from ..core.interpro_pfam_analysis import InterProPfamParameters
from ..core.string_analysis import StringParameters
from ..core.complex_analysis import ComplexParameters
from ..core.mtdna_analysis import MtdnaParameters
from ..core.presence_analysis import (
    build_presence_arguments, export_presence, read_presence_outputs,
)
from ..core.project import Project, ProjectError, create_project, open_project
from ..core.r_runtime import RRuntime
from .analyses_page import AnalysesPage
from .data_page import DataPage
from .database_manager_page import DatabaseManagerPage
from .mapping_worker import MappingWorker
from .reactome_page import ReactomeAnalysisWorker
from .mitocarta_page import MitoCartaAnalysisWorker
from .interpro_pfam_page import InterProPfamAnalysisWorker
from .string_page import StringAnalysisWorker
from .complexes_page import ComplexAnalysisWorker
from .mtdna_page import MtdnaAnalysisWorker
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
        self.kegg_worker: MappingWorker | None = None
        self.reactome_worker: ReactomeAnalysisWorker | None = None
        self.mitocarta_worker: MitoCartaAnalysisWorker | None = None
        self.interpro_pfam_worker: InterProPfamAnalysisWorker | None = None
        self.string_worker: StringAnalysisWorker | None = None
        self._string_threads: list[StringAnalysisWorker] = []
        self.complex_worker: ComplexAnalysisWorker | None = None
        self._complex_threads: list[ComplexAnalysisWorker] = []
        self.mtdna_worker: MtdnaAnalysisWorker | None = None
        self._mtdna_threads: list[MtdnaAnalysisWorker] = []
        self.project_page = ProjectPage()
        self.data_page = DataPage()
        self.database_manager_page = DatabaseManagerPage()
        self.analyses_page = AnalysesPage(self.database_manager_page.manager)
        self.analyses_page.proteomics_qc_page.runtime = self.runtime
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
        kegg=self.analyses_page.kegg_page
        kegg.run_requested.connect(self._run_kegg)
        kegg.open_database_requested.connect(lambda:self.navigation.setCurrentRow(3))
        reactome=self.analyses_page.reactome_page
        reactome.run_requested.connect(self._run_reactome)
        reactome.open_database_requested.connect(lambda:self.navigation.setCurrentRow(3))
        mitocarta=self.analyses_page.mitocarta_page
        mitocarta.run_requested.connect(self._run_mitocarta)
        mitocarta.open_database_requested.connect(lambda:self.navigation.setCurrentRow(3))
        interpro=self.analyses_page.interpro_pfam_page
        interpro.run_requested.connect(self._run_interpro_pfam)
        interpro.open_database_requested.connect(lambda:self.navigation.setCurrentRow(3))
        string_page=self.analyses_page.string_page
        string_page.run_requested.connect(self._run_string)
        string_page.open_database_requested.connect(lambda:self.navigation.setCurrentRow(3))
        complex_page=self.analyses_page.complex_page
        complex_page.run_requested.connect(self._run_complex)
        complex_page.open_database_requested.connect(lambda:self.navigation.setCurrentRow(3))
        mtdna_page=self.analyses_page.mtdna_page
        mtdna_page.run_requested.connect(self._run_mtdna)
        mtdna_page.open_database_requested.connect(lambda:self.navigation.setCurrentRow(3))
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
                (self.go_worker and self.go_worker.isRunning()) or (self.kegg_worker and self.kegg_worker.isRunning()) or
                (self.reactome_worker and self.reactome_worker.isRunning()) or
                (self.mitocarta_worker and self.mitocarta_worker.isRunning()) or
                (self.interpro_pfam_worker and self.interpro_pfam_worker.isRunning()) or
                (self.string_worker and self.string_worker.isRunning()) or
                bool(self._string_threads) or
                (self.complex_worker and self.complex_worker.isRunning()) or
                bool(self._complex_threads) or
                (self.mtdna_worker and self.mtdna_worker.isRunning()) or
                bool(self._mtdna_threads) or
                self.analyses_page.interpro_pfam_page.is_running() or
                self.analyses_page.proteomics_qc_page.is_running() or self.database_manager_page.is_running()):
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

    def _run_kegg(self,parameters:dict,allow_target_outside_background:bool=False)->None:
        if not self.project or self.kegg_worker:return
        run_id=new_run_id()
        try:arguments=prepare_kegg_arguments(self.project,self.database_manager_page.manager,run_id=run_id,allow_target_outside_background=allow_target_outside_background,**parameters)
        except ValueError as error:QMessageBox.warning(self,"KEGG Pathways",str(error));return
        worker=MappingWorker(self.runtime,APPLICATION_ROOT/"r_scripts"/"04_kegg_analysis.R",arguments);self.kegg_worker=worker;self.analyses_page.kegg_page.set_running(True)
        worker.succeeded.connect(lambda stdout,stderr:self._kegg_finished(run_id));worker.failed.connect(lambda message,stderr:self._kegg_failed(message,stderr,parameters));worker.finished.connect(worker.deleteLater);worker.start()
    def _kegg_finished(self,run_id:str)->None:
        self.kegg_worker=None;self.analyses_page.kegg_page.set_running(False)
        try:self.analyses_page.kegg_page.show_outputs(read_kegg_outputs(self.project));self.scripts_page.set_project_scripts(self.project.root/"scripts"/"runs")
        except RuntimeError as error:QMessageBox.warning(self,"KEGG results",str(error))
    def _kegg_failed(self,message:str,stderr:str,parameters:dict)->None:
        self.kegg_worker=None;self.analyses_page.kegg_page.set_running(False);self.scripts_page.result_view.setPlainText(stderr)
        import re
        match=re.search(r"(\d+) target KEGG gene\(s\) are outside",message+stderr)
        if match:
            if self._confirm_target_adjustment(int(match.group(1))):self._run_kegg(parameters,True)
        else:QMessageBox.warning(self,"KEGG analysis failed",message)
    def _confirm_target_adjustment(self,count:int)->bool:
        box=QMessageBox(self);box.setWindowTitle("Target genes outside background");box.setText(f"{count} target genes are not present in the selected background.\n\nPathway enrichment requires the target set to be contained within the background.\n\nContinue using only target genes present in the background?")
        proceed=box.addButton("Continue",QMessageBox.ButtonRole.AcceptRole);box.addButton("Cancel",QMessageBox.ButtonRole.RejectRole);box.exec();return box.clickedButton()==proceed

    def _run_reactome(self,parameters:dict,allow_target_outside_background:bool=False)->None:
        if not self.project or self.reactome_worker:return
        run_id=new_run_id();options=ReactomeParameters(**parameters,allow_target_outside_background=allow_target_outside_background)
        worker=ReactomeAnalysisWorker(self.project,self.database_manager_page.manager,self.runtime,run_id,options)
        self.reactome_worker=worker;self.analyses_page.reactome_page.set_running(True)
        worker.succeeded.connect(self._reactome_finished)
        worker.target_outside_background.connect(lambda details,p=parameters:self._reactome_target_outside(details,p))
        worker.failed.connect(self._reactome_failed);worker.finished.connect(worker.deleteLater);worker.start()

    def _reactome_finished(self,outputs)->None:
        self.reactome_worker=None;self.analyses_page.reactome_page.set_running(False)
        self.analyses_page.reactome_page.show_outputs(outputs)
        if self.project:self.scripts_page.set_project_scripts(self.project.root/"scripts"/"runs")

    def _reactome_failed(self,message:str)->None:
        self.reactome_worker=None;self.analyses_page.reactome_page.set_running(False)
        QMessageBox.warning(self,"Reactome analysis failed",message)

    def _reactome_target_outside(self,details:dict,parameters:dict)->None:
        self.reactome_worker=None;self.analyses_page.reactome_page.set_running(False)
        if self._confirm_reactome_adjustment(details):self._run_reactome(parameters,True)

    def _confirm_reactome_adjustment(self,details:dict)->bool:
        count=len(details.get("entities_outside_background",[]));box=QMessageBox(self)
        box.setWindowTitle("Target outside background")
        box.setText(f"{count} canonical Reactome entity/entities are outside the selected background.\n\nThe analysis cannot continue without a decision. Continue will restrict the target to entities present in the selected background and record that adjustment in the run.")
        proceed=box.addButton("Continue",QMessageBox.ButtonRole.AcceptRole);box.addButton("Cancel",QMessageBox.ButtonRole.RejectRole);box.exec()
        return box.clickedButton()==proceed

    def _run_mitocarta(self,parameters:dict,allow_target_outside_background:bool=False)->None:
        if not self.project or self.mitocarta_worker:return
        run_id=new_run_id();options=MitoCartaParameters(**parameters,allow_target_outside_background=allow_target_outside_background)
        worker=MitoCartaAnalysisWorker(self.project,self.database_manager_page.manager,self.runtime,run_id,options)
        self.mitocarta_worker=worker;self.analyses_page.mitocarta_page.set_running(True)
        worker.succeeded.connect(self._mitocarta_finished)
        worker.target_outside_background.connect(lambda details,p=parameters:self._mitocarta_target_outside(details,p))
        worker.failed.connect(self._mitocarta_failed);worker.finished.connect(worker.deleteLater);worker.start()

    def _mitocarta_finished(self,outputs)->None:
        self.mitocarta_worker=None;self.analyses_page.mitocarta_page.set_running(False)
        self.analyses_page.mitocarta_page.show_outputs(outputs)
        if self.project:self.scripts_page.set_project_scripts(self.project.root/"scripts"/"runs")

    def _mitocarta_failed(self,message:str)->None:
        self.mitocarta_worker=None;self.analyses_page.mitocarta_page.set_running(False)
        QMessageBox.warning(self,"MitoCarta analysis failed",message)

    def _mitocarta_target_outside(self,details:dict,parameters:dict)->None:
        self.mitocarta_worker=None;self.analyses_page.mitocarta_page.set_running(False)
        if self._confirm_mitocarta_adjustment(details):self._run_mitocarta(parameters,True)

    def _confirm_mitocarta_adjustment(self,details:dict)->bool:
        count=len(details.get("entities_outside_background",[]));target=details.get("initial_target_size",details.get("target_size","unknown"));background=details.get("initial_background_size",details.get("background_size","unknown"));box=QMessageBox(self)
        box.setWindowTitle("Target outside background")
        box.setText(f"{count} gene-resolved target entity/entities are outside the selected background (target size: {target}; background size: {background}).\n\nMitoCarta enrichment requires the target to be contained within the background. Continue will restrict the target to entities present in the background and record that adjustment in the run.")
        proceed=box.addButton("Continue",QMessageBox.ButtonRole.AcceptRole);box.addButton("Cancel",QMessageBox.ButtonRole.RejectRole);box.exec()
        return box.clickedButton()==proceed

    def _run_interpro_pfam(self,parameters:dict,allow_target_outside_background:bool=False)->None:
        if not self.project or self.interpro_pfam_worker:return
        options=InterProPfamParameters(**parameters,allow_target_outside_background=allow_target_outside_background);worker=InterProPfamAnalysisWorker(self.project,self.database_manager_page.manager,self.runtime,new_run_id(),options);self.interpro_pfam_worker=worker;self.analyses_page.interpro_pfam_page.set_running(True)
        worker.succeeded.connect(self._interpro_pfam_finished);worker.target_outside_background.connect(lambda details,p=parameters:self._interpro_pfam_target_outside(details,p));worker.failed.connect(self._interpro_pfam_failed);worker.finished.connect(worker.deleteLater);worker.start()
    def _interpro_pfam_finished(self,outputs):
        self.interpro_pfam_worker=None;self.analyses_page.interpro_pfam_page.set_running(False);self.analyses_page.interpro_pfam_page.show_outputs(outputs)
        if self.project:self.scripts_page.set_project_scripts(self.project.root/"scripts/runs")
    def _interpro_pfam_failed(self,message):self.interpro_pfam_worker=None;self.analyses_page.interpro_pfam_page.set_running(False);QMessageBox.warning(self,"InterPro/Pfam analysis failed",message)
    def _interpro_pfam_target_outside(self,details,parameters):
        self.interpro_pfam_worker=None;self.analyses_page.interpro_pfam_page.set_running(False)
        if self._confirm_interpro_pfam_adjustment(details):self._run_interpro_pfam(parameters,True)
    def _confirm_interpro_pfam_adjustment(self,details):
        count=len(details.get("entities_outside_background",[]));target=details.get("initial_target_size","unknown");background=details.get("initial_background_size","unknown");box=QMessageBox(self);box.setWindowTitle("Target outside background");box.setText(f"{count} target protein(s) are outside the selected background (target size: {target}; background size: {background}).\n\nEnrichment requires the target to be contained within the background. Continue will authorize the backend to restrict the target and record the decision.");proceed=box.addButton("Continue",QMessageBox.ButtonRole.AcceptRole);box.addButton("Cancel",QMessageBox.ButtonRole.RejectRole);box.exec();return box.clickedButton()==proceed

    def _run_string(self, parameters: dict, allow_target_outside_background: bool = False) -> None:
        if not self.project or self.string_worker:
            return
        try:
            options = StringParameters(**parameters, allow_target_outside_background=allow_target_outside_background)
        except (TypeError, ValueError) as error:
            QMessageBox.warning(self, "STRING analysis", str(error))
            return
        worker = StringAnalysisWorker(self.project, self.database_manager_page.manager, self.runtime, new_run_id(), options)
        self.string_worker = worker
        self._string_threads.append(worker)
        self.analyses_page.string_page.set_running(True)
        worker.succeeded.connect(self._string_finished)
        worker.failed.connect(self._string_failed)
        worker.expansion_limit.connect(self._string_expansion_limit)
        worker.target_outside.connect(lambda details, p=parameters: self._string_target_outside(details, p))
        worker.finished.connect(lambda w=worker: self._string_threads.remove(w))
        worker.finished.connect(worker.deleteLater)
        worker.start()

    def _string_finished(self, outputs) -> None:
        self.string_worker = None
        self.analyses_page.string_page.set_running(False)
        self.analyses_page.string_page.show_outputs(outputs)
        if self.project:
            self.scripts_page.set_project_scripts(self.project.root / "scripts" / "runs")

    def _string_failed(self, message: str) -> None:
        self.string_worker = None
        self.analyses_page.string_page.set_running(False)
        QMessageBox.warning(self, "STRING analysis failed", message)

    def _string_expansion_limit(self, details: dict) -> None:
        self.string_worker = None
        self.analyses_page.string_page.set_running(False)
        QMessageBox.warning(self, "STRING expansion safety limit",
            f"The network reaches {details.get('observed')} external proteins at hop {details.get('hop')}, "
            f"above the configured limit of {details.get('limit')} (score ≥ {details.get('threshold')}, "
            f"{details.get('network_type')} network). Raise the confidence threshold, reduce the maximum hop "
            "or narrow the seed set, then run again. No network was silently truncated.")

    def _string_target_outside(self, details: dict, parameters: dict) -> None:
        self.string_worker = None
        self.analyses_page.string_page.set_running(False)
        box = QMessageBox(self)
        box.setWindowTitle("STRING target outside background")
        box.setText(f"{len(details.get('entities_outside_background', []))} target proteins are outside the selected "
            f"background (target {details.get('target_size')}; background {details.get('background_size')}). "
            "Continue will authorize the backend to restrict the target and record the adjustment.")
        proceed = box.addButton("Continue", QMessageBox.ButtonRole.AcceptRole)
        box.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
        box.exec()
        if box.clickedButton() == proceed:
            self._run_string(parameters, True)
    def _run_complex(self, parameters: dict, allow_target_outside_background: bool = False) -> None:
        if not self.project or self.complex_worker:
            return
        try:
            options = ComplexParameters(**parameters, allow_target_outside_background=allow_target_outside_background)
        except (TypeError, ValueError) as error:
            QMessageBox.warning(self, "Complex analysis", str(error))
            return
        worker = ComplexAnalysisWorker(self.project, self.database_manager_page.manager, self.runtime, new_run_id(), options)
        self.complex_worker = worker
        self._complex_threads.append(worker)
        self.analyses_page.complex_page.set_running(True)
        worker.succeeded.connect(self._complex_finished)
        worker.failed.connect(self._complex_failed)
        worker.target_outside.connect(lambda details, p=parameters: self._complex_target_outside(details, p))
        worker.finished.connect(lambda w=worker: self._complex_threads.remove(w))
        worker.finished.connect(worker.deleteLater)
        worker.start()

    def _complex_finished(self, outputs) -> None:
        self.complex_worker = None
        self.analyses_page.complex_page.set_running(False)
        self.analyses_page.complex_page.show_outputs(outputs)
        if self.project:
            self.scripts_page.set_project_scripts(self.project.root / "scripts" / "runs")

    def _complex_failed(self, message: str) -> None:
        self.complex_worker = None
        self.analyses_page.complex_page.set_running(False)
        QMessageBox.warning(self, "Complex analysis failed", message)

    def _complex_target_outside(self, details: dict, parameters: dict) -> None:
        self.complex_worker = None
        self.analyses_page.complex_page.set_running(False)
        box = QMessageBox(self)
        box.setWindowTitle("Complex target outside background")
        box.setText(f"{len(details.get('entities_outside_background', []))} target proteins are outside the selected "
            f"background (target size: {details.get('initial_target_size')}; background size: "
            f"{details.get('initial_background_size')}). Continue authorizes the backend to restrict the target.")
        proceed = box.addButton("Continue", QMessageBox.ButtonRole.AcceptRole)
        box.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
        box.exec()
        if box.clickedButton() == proceed:
            self._run_complex(parameters, True)
    def _run_mtdna(self, parameters: dict, allow_target_outside_background: bool = False) -> None:
        if not self.project or self.mtdna_worker:
            return
        try:
            options = MtdnaParameters(**parameters,
                allow_target_outside_background=allow_target_outside_background)
        except (TypeError, ValueError) as error:
            QMessageBox.warning(self, "mtDNA Evidence analysis", str(error))
            return
        worker = MtdnaAnalysisWorker(self.project, self.database_manager_page.manager,
            self.runtime, new_run_id(), options)
        self.mtdna_worker = worker
        self._mtdna_threads.append(worker)
        self.analyses_page.mtdna_page.set_running(True)
        worker.succeeded.connect(self._mtdna_finished)
        worker.failed.connect(self._mtdna_failed)
        worker.target_outside.connect(lambda details, p=parameters:self._mtdna_target_outside(details,p))
        worker.finished.connect(lambda w=worker:self._mtdna_threads.remove(w))
        worker.finished.connect(worker.deleteLater)
        worker.start()

    def _mtdna_finished(self, outputs) -> None:
        self.mtdna_worker = None
        self.analyses_page.mtdna_page.set_running(False)
        self.analyses_page.mtdna_page.show_outputs(outputs)
        if self.project:self.scripts_page.set_project_scripts(self.project.root / "scripts" / "runs")

    def _mtdna_failed(self, message: str) -> None:
        self.mtdna_worker = None
        self.analyses_page.mtdna_page.set_running(False)
        QMessageBox.warning(self, "mtDNA Evidence analysis failed", message)

    def _mtdna_target_outside(self, details: dict, parameters: dict) -> None:
        self.mtdna_worker = None
        self.analyses_page.mtdna_page.set_running(False)
        box = QMessageBox(self)
        box.setWindowTitle("mtDNA target outside background")
        box.setText(f"{len(details.get('entities_outside_background', []))} target entities are outside the selected "
            f"background (target: {details.get('initial_target_size')}; background: "
            f"{details.get('initial_background_size')}). Continue will restrict the target and record the adjustment.")
        proceed = box.addButton("Continue", QMessageBox.ButtonRole.AcceptRole)
        box.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
        box.exec()
        if box.clickedButton() == proceed:self._run_mtdna(parameters,True)
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
