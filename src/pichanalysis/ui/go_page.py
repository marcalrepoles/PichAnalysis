from __future__ import annotations

from pathlib import Path

import pandas as pd
from PySide6.QtCore import Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox, QHBoxLayout,
    QLabel, QLineEdit, QPlainTextEdit, QProgressBar, QPushButton, QSpinBox,
    QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from ..core.go_analysis import GOOutputs, available_sets, go_readiness, list_go_runs
from ..core.project import Project


class GOPage(QWidget):
    run_requested=Signal(dict); export_table_requested=Signal(str); export_workbook_requested=Signal(); export_graph_requested=Signal(str); open_folder_requested=Signal()

    def __init__(self) -> None:
        super().__init__(); self.project:Project|None=None; self.outputs:GOOutputs|None=None
        self.target=QComboBox(); self.background=QComboBox(); self.manual_rows=QLineEdit(); self.manual_rows.setPlaceholderText("Linhas separadas por vírgula")
        self.ontology=QComboBox(); self.ontology.addItem("Todos — BP, MF e CC",["BP","MF","CC"])
        for value,label in (("BP","BP — Biological Process"),("MF","MF — Molecular Function"),("CC","CC — Cellular Component")): self.ontology.addItem(label,[value])
        self.evidence=QComboBox(); self.evidence.addItem("Todas as evidências","all"); self.evidence.addItem("Excluir IEA","exclude_iea"); self.evidence.addItem("Somente evidências experimentais","experimental")
        self.fdr=QDoubleSpinBox(); self.fdr.setRange(.000001,1); self.fdr.setValue(.05); self.fdr.setDecimals(4)
        self.pvalue=QDoubleSpinBox(); self.pvalue.setRange(.000001,1); self.pvalue.setValue(1); self.pvalue.setDecimals(4)
        self.min_count=QSpinBox(); self.min_count.setRange(1,100000); self.min_count.setValue(3)
        self.top_n=QSpinBox(); self.top_n.setRange(1,200); self.top_n.setValue(20)
        self.simplify=QCheckBox("Simplificar termos redundantes (preserva originais)"); self.simplify_cutoff=QDoubleSpinBox(); self.simplify_cutoff.setRange(.1,1); self.simplify_cutoff.setValue(.7)
        form=QFormLayout(); form.addRow("1. Conjunto analisado",self.target); form.addRow("Seleção manual",self.manual_rows); form.addRow("2. Background",self.background); form.addRow("3. Ontologia",self.ontology); form.addRow("4. Evidências",self.evidence)
        advanced=QGroupBox("5. Opções avançadas de enrichment"); advanced.setCheckable(True); advanced.setChecked(False); advanced_layout=QFormLayout(advanced)
        advanced_layout.addRow("FDR cutoff",self.fdr); advanced_layout.addRow("p-value cutoff",self.pvalue); advanced_layout.addRow("Minimum gene count",self.min_count); advanced_layout.addRow("Top N",self.top_n); advanced_layout.addRow(self.simplify); advanced_layout.addRow("Cutoff de similaridade",self.simplify_cutoff)
        self.explanation=QLabel("Frequência mostra quantas proteínas possuem cada anotação. Enriquecimento testa se uma anotação aparece mais do que o esperado em relação ao background. FDR reduz falsos positivos entre muitos testes."); self.explanation.setWordWrap(True)
        self.readiness=QLabel("Abra um projeto."); self.run_button=QPushButton("6. Executar GO"); self.run_button.setEnabled(False); self.progress=QProgressBar(); self.progress.setRange(0,0); self.progress.hide()
        config=QWidget(); config_layout=QVBoxLayout(config); config_layout.addLayout(form); config_layout.addWidget(advanced); config_layout.addWidget(self.explanation); config_layout.addWidget(self.readiness); config_layout.addWidget(self.run_button); config_layout.addWidget(self.progress)
        self.source=QLabel("Fonte GO desta execução: —"); self.summary=QLabel("Nenhum resultado disponível."); self.summary.setWordWrap(True)
        self.table_picker=QComboBox(); self.table=QTableWidget(); self.detail=QPlainTextEdit(); self.detail.setReadOnly(True)
        self.export_table=QPushButton("Exportar tabela atual"); self.export_workbook=QPushButton("Exportar workbook"); self.open_folder=QPushButton("Abrir pasta GO")
        actions=QHBoxLayout(); [actions.addWidget(x) for x in (self.export_table,self.export_workbook,self.open_folder)]
        results=QWidget(); rl=QVBoxLayout(results); rl.addWidget(self.source); rl.addWidget(self.summary); rl.addWidget(self.table_picker); rl.addWidget(self.table,1); rl.addWidget(self.detail); rl.addLayout(actions)
        self.graph_picker=QComboBox(); self.graph=QLabel("Nenhum gráfico."); self.export_graph=QPushButton("Exportar gráfico")
        graphs=QWidget(); gl=QVBoxLayout(graphs); gl.addWidget(self.graph_picker); gl.addWidget(self.graph,1); gl.addWidget(self.export_graph)
        self.history=QComboBox(); history=QWidget(); hl=QVBoxLayout(history); hl.addWidget(QLabel("Execuções preservadas em disco")); hl.addWidget(self.history); hl.addStretch()
        tabs=QTabWidget(); tabs.addTab(config,"Configuração"); tabs.addTab(results,"Resultados"); tabs.addTab(graphs,"Gráficos"); tabs.addTab(history,"Histórico")
        layout=QVBoxLayout(self); layout.addWidget(tabs)
        self.run_button.clicked.connect(self._emit_run); self.table_picker.currentIndexChanged.connect(self._load_table); self.table.currentCellChanged.connect(self._show_detail); self.graph_picker.currentIndexChanged.connect(self._show_graph)
        self.export_table.clicked.connect(lambda:self.export_table_requested.emit(str(self.table_picker.currentData() or ""))); self.export_workbook.clicked.connect(self.export_workbook_requested); self.export_graph.clicked.connect(lambda:self.export_graph_requested.emit(str(self.graph_picker.currentData() or ""))); self.open_folder.clicked.connect(self.open_folder_requested)

    def set_project(self,project:Project|None)->None:
        self.project=project; self.target.clear(); self.background.clear(); self.history.clear()
        state=go_readiness(project)
        if project and state.ready:
            for key,label in available_sets(project).items(): self.target.addItem(label,key)
            for key,label in available_sets(project).items():
                if key in {"all_experiment","mapped","manual"}: self.background.addItem(label,key)
            index=self.background.findData("mapped"); self.background.setCurrentIndex(max(0,index)); self.history.addItems(list_go_runs(project))
        self.readiness.setText(state.reason); self.run_button.setEnabled(state.ready)

    def _emit_run(self)->None:
        rows=[]
        try: rows=[int(x.strip()) for x in self.manual_rows.text().split(",") if x.strip()]
        except ValueError: self.readiness.setText("Seleção manual deve conter números de linha."); return
        self.run_requested.emit({"target_selection":self.target.currentData(),"background_selection":self.background.currentData(),"manual_rows":rows,"ontologies":self.ontology.currentData(),"evidence_filter":self.evidence.currentData(),"fdr_cutoff":self.fdr.value(),"p_cutoff":self.pvalue.value(),"min_count":self.min_count.value(),"top_n":self.top_n.value(),"simplify":self.simplify.isChecked(),"simplify_cutoff":self.simplify_cutoff.value()})

    def set_running(self,running:bool)->None:self.progress.setVisible(running);self.run_button.setEnabled(not running)

    def show_outputs(self,outputs:GOOutputs)->None:
        self.outputs=outputs;m=outputs.metadata;self.source.setText(f"Fonte GO desta execução: {m.get('organism')} | {m.get('orgdb')} {m.get('orgdb_version')} | GO.db {m.get('go_db_version')} | Evidências: {m.get('evidence_filter')}")
        significant=int(outputs.summary.get("significant_terms",pd.Series(dtype=int)).sum());self.summary.setText(f"Target: {m.get('target_count')} | Background: {m.get('background_count')} | Com GO: {m.get('annotated_count')} | Sem GO: {m.get('unannotated_count')} | Termos testados: {int(outputs.summary.get('terms_tested',pd.Series(dtype=int)).sum())} | Significativos: {significant}")
        self.table_picker.clear()
        for path in outputs.tables:self.table_picker.addItem(path.name,str(path))
        self.graph_picker.clear()
        for path in outputs.graphs:self.graph_picker.addItem(path.stem,str(path))
        self.history.clear(); self.history.addItems(list_go_runs(self.project))

    def _load_table(self)->None:
        path=self.table_picker.currentData()
        if not path:return
        frame=pd.read_csv(path).head(500);self.table.setRowCount(len(frame));self.table.setColumnCount(len(frame.columns));self.table.setHorizontalHeaderLabels(list(frame.columns))
        for r,row in enumerate(frame.itertuples(index=False,name=None)):
            for c,value in enumerate(row):self.table.setItem(r,c,QTableWidgetItem("" if pd.isna(value) else str(value)))

    def _show_detail(self,row:int,*_args)->None:
        if row<0:return
        headers=[self.table.horizontalHeaderItem(i).text() for i in range(self.table.columnCount())];wanted=("GO_ID","Description","GO_term","Ontology","ontology","p.adjust","Count","genes","input_id","gene_symbol","uniprot_accession")
        lines=[f"{x}: {self.table.item(row,headers.index(x)).text()}" for x in wanted if x in headers and self.table.item(row,headers.index(x))]
        if self.outputs and "input_id" in headers and self.table.item(row,headers.index("input_id")):
            protein=self.table.item(row,headers.index("input_id")).text();annotations=self.outputs.annotations
            related=annotations[annotations["input_id"].astype(str)==protein]
            if len(related): lines.append("Termos da proteína:\n"+"\n".join(f"{item.ontology}: {item.GO_ID} — {item.GO_term}" for item in related.itertuples()))
        self.detail.setPlainText("\n\n".join(lines))

    def _show_graph(self)->None:
        path=self.graph_picker.currentData()
        if path:self.graph.setPixmap(QPixmap(path).scaled(900,600))
