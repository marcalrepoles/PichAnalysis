"""Offline, auditable navigation among already persisted analysis results."""
from __future__ import annotations

import csv

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QCheckBox, QFileDialog, QHBoxLayout, QLabel,
    QLineEdit, QMessageBox, QPlainTextEdit, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget)

from ..core.cross_module_integration import (open_cross_module_index,
    rebuild_cross_module_index)
from ..core.entity_identity import MatchStatus


class CrossModuleExplorer(QWidget):
    open_target_requested = Signal(str, str, str, str, str)
    biological_context_requested = Signal(object)

    def __init__(self):
        super().__init__()
        self.project = None
        self.index = None
        self.context = None
        self.matches = []
        self.source_records = []
        layout = QVBoxLayout(self)
        self.help = QLabel("Cross-module navigation links existing analysis results. It does not recompute analyses or reinterpret an experimental feature as multiple quantitative measurements.")
        self.help.setWordWrap(True); layout.addWidget(self.help)
        search_row = QHBoxLayout()
        self.search = QLineEdit(); self.search.setPlaceholderText("Feature ID, source row, original ID, UniProt, Gene Symbol or NCBI Gene")
        self.search_button = QPushButton("Search local index")
        self.rebuild_button = QPushButton("Rebuild index")
        for widget in (self.search, self.search_button, self.rebuild_button): search_row.addWidget(widget)
        layout.addLayout(search_row)
        self.source_table = QTableWidget(0, 6)
        self.source_table.setHorizontalHeaderLabels(["Module", "Run", "Type", "Feature", "Source row", "Identifier"])
        layout.addWidget(QLabel("Select a source feature or entity")); layout.addWidget(self.source_table, 1)
        self.identity = QPlainTextEdit(); self.identity.setReadOnly(True)
        layout.addWidget(self.identity)
        self.cross_lineage = QCheckBox("Include matches from other input lineages")
        layout.addWidget(self.cross_lineage)
        caveat = QLabel("A biological identifier match from another input lineage does not represent the same experimental feature.")
        caveat.setWordWrap(True); layout.addWidget(caveat)
        self.matches_table = QTableWidget(0, 8)
        self.matches_table.setHorizontalHeaderLabels(["Module", "Run", "Run date", "Database snapshot", "Match", "Match basis", "Target", "Details"])
        layout.addWidget(self.matches_table, 2)
        self.details = QPlainTextEdit(); self.details.setReadOnly(True)
        layout.addWidget(self.details)
        actions = QHBoxLayout()
        self.open_button = QPushButton("Open in analysis")
        self.open_button.setEnabled(False)
        self.export_button = QPushButton("Export cross-module summary...")
        actions.addWidget(self.open_button); actions.addWidget(self.export_button)
        self.biological_context_button = QPushButton("Biological context...")
        self.biological_context_button.setEnabled(False)
        actions.addWidget(self.biological_context_button)
        layout.addLayout(actions)
        self.status = QLabel("Open a project to explore existing runs.")
        layout.addWidget(self.status)
        self.search.returnPressed.connect(self._search)
        self.search_button.clicked.connect(self._search)
        self.rebuild_button.clicked.connect(self._rebuild)
        self.source_table.itemSelectionChanged.connect(self._source_selected)
        self.matches_table.itemSelectionChanged.connect(self._match_selected)
        self.cross_lineage.toggled.connect(self._resolve)
        self.open_button.clicked.connect(self._open_target)
        self.export_button.clicked.connect(self._export)
        self.biological_context_button.clicked.connect(self._open_biological_context)

    def set_project(self, project):
        self.project = project
        self.index = None; self.context = None; self.matches = []
        self.biological_context_button.setEnabled(False)
        self.source_table.setRowCount(0); self.matches_table.setRowCount(0)
        self.identity.clear(); self.details.clear()
        self.status.setText("Search the local index or select a feature in another analysis."
            if project else "Open a project to explore existing runs.")

    def _ensure_index(self):
        if not self.project:
            return False
        try:
            self.index = open_cross_module_index(self.project)
            return True
        except Exception as error:
            self.status.setText(f"Could not open cross-module index: {error}")
            return False

    def _rebuild(self):
        if not self.project: return
        try:
            self.index = rebuild_cross_module_index(self.project)
            self.status.setText("Local cross-module index rebuilt from persisted runs.")
            if self.search.text().strip(): self._search()
        except Exception as error:
            self.status.setText(f"Index rebuild failed: {error}")

    def _search(self):
        if not self._ensure_index(): return
        self.source_records = self.index.search(self.search.text())
        self.source_table.setRowCount(len(self.source_records))
        for row, record in enumerate(self.source_records):
            for column, value in enumerate((record["module"], record["run_id"], record["record_type"],
                record["feature_id"], record["source_row"], record["original_identifier"])):
                self.source_table.setItem(row, column, QTableWidgetItem(str(value or "")))
        self.status.setText(f"{len(self.source_records)} source records found. Select one; ambiguous results are never chosen automatically.")
        self.context = None; self.matches = []
        self.matches_table.setRowCount(0); self.identity.clear(); self.details.clear()
        self.open_button.setEnabled(False)

    def explore_feature(self, module_id, run_id, feature_id):
        self.search.setText(str(feature_id))
        self._search()
        candidates = [row for row, record in enumerate(self.source_records)
            if record["module"] == module_id and record["run_id"] == run_id
            and (feature_id in (record["feature_id"], record["original_identifier"], record["target"])
            or str(record["source_row"]) == str(feature_id))]
        preferred_type = {"proteomics_qc": "feature_metadata", "differential": "all_results",
            "presence_absence": "classification", "mapping": "catalog", "go": "annotations",
            "kegg": "mapping", "reactome": "mapping", "mitocarta": "mapping",
            "domains": "mapping", "string": "mapping", "complexes": "mapping",
            "mtdna_evidence": "mapping"}.get(module_id)
        preferred = [row for row in candidates
            if self.source_records[row]["record_type"] == preferred_type]
        if len(preferred) == 1:
            candidates = preferred
        if len(candidates) == 1:
            self.source_table.selectRow(candidates[0])
            self._source_selected()
            return True
        if len(candidates) > 1:
            self.status.setText("Multiple frozen source records match. Select the intended row explicitly.")
            return False
        if module_id == "mapping" and self.index:
            status = self.index.run_status(module_id, run_id)
            if status and status[0] != "Ready":
                self.status.setText(status[1])
                return False
        self.status.setText("No frozen source record was found for the selected feature.")
        return False

    def _source_selected(self):
        row = self.source_table.currentRow()
        if row < 0 or row >= len(self.source_records) or not self.index: return
        record = self.source_records[row]
        self.context = self.index.context(record["id"])
        context = self.context
        self.biological_context_button.setEnabled(True)
        lines = [f"Source module: {context.source_module}", f"Source run: {context.source_run_id}",
            f"Source result type: {context.source_result_type}",
            "Experimental identity:", f"  Feature ID: {context.feature_id or 'Not available'}",
            f"  Source row: {context.source_row or 'Not available'}",
            f"  Original identifier: {context.original_identifier or 'Not available'}",
            "Protein identity:", f"  UniProt: {', '.join(context.uniprot_accessions) or 'Not available'}",
            f"  Protein group members: {', '.join(context.protein_group_members) or 'Not available'}",
            "Gene identity:", f"  Gene Symbols: {', '.join(context.gene_symbols) or 'Not available'}",
            f"  NCBI Gene: {', '.join(context.ncbi_gene_ids) or 'Not available'}",
            "Database-specific identity:", f"  IDs: {context.module_specific_ids or 'Not available'}",
            f"  Snapshot: {context.database_snapshot_provenance or 'Not available'}",
            f"Input lineage: {context.input_lineage.frozen_input_hash if context.input_lineage else 'Not available'}"]
        self.identity.setPlainText("\n".join(lines))
        self._resolve()

    def _resolve(self):
        if not self.context or not self.index: return
        self.matches = self.index.matches(self.context,
            include_other_lineages=self.cross_lineage.isChecked())
        self.matches_table.setRowCount(len(self.matches))
        for row, match in enumerate(self.matches):
            for column, value in enumerate((match.target_module, match.target_run_id,
                match.run_date, match.database_snapshot, match.status.value, match.basis,
                match.target_entity, match.artifact_source)):
                self.matches_table.setItem(row, column, QTableWidgetItem(str(value or "")))
        self.status.setText(f"{len(self.matches)} matching records in persisted runs." if self.matches else
            "No compatible analysis run or unambiguous cross-module identifier is available.")
        self.details.clear()
        self.open_button.setEnabled(False)

    def _match_selected(self):
        row = self.matches_table.currentRow()
        if row < 0 or row >= len(self.matches): return
        match = self.matches[row]
        self.open_button.setEnabled(match.status in (MatchStatus.EXACT, MatchStatus.MAPPED, MatchStatus.AMBIGUOUS))
        lines = [f"Source: {self.context.source_module} / {self.context.source_run_id}",
            f"Target: {match.target_module} / {match.target_run_id} / {match.target_entity}",
            f"Match: {match.status.value}", f"Basis: {match.basis}",
            f"Source identifier: {match.source_identifier}",
            f"Target identifier: {match.target_identifier}",
            f"Input lineage compatible: {match.lineage_compatible}",
            f"Database snapshot: {match.database_snapshot or 'Not available'}",
            f"Artifact source: {match.artifact_source}"]
        if match.status == MatchStatus.AMBIGUOUS:
            lines.append("Multiple compatible target entities were found. Select one to continue.")
        self.details.setPlainText("\n".join(lines))

    def _open_target(self):
        row = self.matches_table.currentRow()
        if row < 0 or row >= len(self.matches): return
        match = self.matches[row]
        if match.status in (MatchStatus.INCOMPLETE_RUN, MatchStatus.NO_COMPATIBLE_RUN,
            MatchStatus.NOT_FOUND): return
        self.open_target_requested.emit(match.target_module, match.target_run_id,
            match.target_entity, match.details.get("record_type", ""),
            match.details.get("source_row", ""))

    def _open_biological_context(self):
        if self.context:
            self.biological_context_requested.emit(self.context)
    def _export(self):
        if not self.context: return
        path, _ = QFileDialog.getSaveFileName(self, "Export cross-module summary",
            "cross_module_summary.csv", "CSV files (*.csv)")
        if not path: return
        with open(path, "w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.writer(stream)
            writer.writerow(("Source module", "Source run", "Source feature", "Source row",
                "Target module", "Target run", "Target entity", "Match status", "Match basis",
                "Input lineage compatible", "Snapshot", "Artifact"))
            for match in self.matches:
                writer.writerow((self.context.source_module, self.context.source_run_id,
                    self.context.feature_id or "", self.context.source_row or "",
                    match.target_module, match.target_run_id, match.target_entity,
                    match.status.value, match.basis, match.lineage_compatible,
                    match.database_snapshot, match.artifact_source))
