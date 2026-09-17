from __future__ import annotations

from collections import Counter
from typing import Any

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QHeaderView, QLabel, QLineEdit,
    QPushButton, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from ..core.column_mapping import (
    ColumnRole, MappingValidation, QuantificationType, QUANTIFICATION_LABELS, ROLE_LABELS,
)
from ..core.identifier_detection import IDENTIFIER_LABELS, IdentifierType
from ..core.importer import ImportResult


class DataPage(QWidget):
    PREVIEW_ROWS = 200
    save_mapping_requested = Signal(dict)

    def __init__(self) -> None:
        super().__init__()
        self._detections: dict[str, dict[str, Any]] = {}
        self.import_button = QPushButton("Import XLSX, CSV, or TSV")
        self.summary = QLabel("Abra ou crie um projeto para importar dados.")
        self.summary.setWordWrap(True)
        self.preview = QTableWidget()
        self.preview.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.preview.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.columns = QTableWidget(0, 4)
        self.columns.setHorizontalHeaderLabels(["Name", "Type", "Present", "Missing"])
        self.columns.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.columns.horizontalHeader().setStretchLastSection(True)

        self.mapping = QTableWidget(0, 8)
        self.mapping.setHorizontalHeaderLabels([
            "Column", "Role", "Primary", "Identifier type", "Condition",
            "Replicate", "Quantification type", "Detection",
        ])
        self.mapping.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.mapping.verticalHeader().setVisible(False)
        self.mapping.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.mapping.horizontalHeader().setSectionResizeMode(7, QHeaderView.ResizeMode.Stretch)
        self.apply_button = QPushButton("Aplicar sugestões automáticas")
        self.save_button = QPushButton("Save and validate configuration")
        self.status = QLabel("Configuration: Not configured")
        self.status.setWordWrap(True)
        self.design_summary = QLabel()
        self.design_summary.setWordWrap(True)

        tabs = QTabWidget()
        preview_page = QWidget()
        preview_layout = QVBoxLayout(preview_page)
        preview_layout.addWidget(QLabel("Pré-visualização (até 200 linhas)"))
        preview_layout.addWidget(self.preview, 3)
        preview_layout.addWidget(QLabel("Detected columns"))
        preview_layout.addWidget(self.columns, 2)
        mapping_page = QWidget()
        mapping_layout = QVBoxLayout(mapping_page)
        mapping_layout.addWidget(QLabel("Column configuration — suggestions can be freely corrected."))
        mapping_layout.addWidget(self.mapping)
        mapping_layout.addWidget(self.apply_button)
        mapping_layout.addWidget(self.save_button)
        mapping_layout.addWidget(self.status)
        mapping_layout.addWidget(self.design_summary)
        tabs.addTab(preview_page, "Tabela e colunas")
        tabs.addTab(mapping_page, "Column configuration")

        layout = QVBoxLayout(self)
        layout.addWidget(self.import_button)
        layout.addWidget(self.summary)
        layout.addWidget(tabs, 1)
        self.apply_button.clicked.connect(self._apply_suggestions)
        self.save_button.clicked.connect(lambda: self.save_mapping_requested.emit(self.mapping_values()))

    @staticmethod
    def _combo(labels: dict[Any, str], selected: str) -> QComboBox:
        combo = QComboBox()
        for key, label in labels.items():
            combo.addItem(label, key.value)
        index = combo.findData(selected)
        combo.setCurrentIndex(max(index, 0))
        return combo

    def show_result(self, result: ImportResult, mapping: dict[str, dict[str, Any]] | None = None) -> None:
        frame = result.dataframe
        sheet = f" | Planilha: {result.sheet}" if result.sheet else ""
        self.summary.setText(
            f"Arquivo: {result.original_path.name} | Formato: {result.source_format.upper()}"
            f"{sheet} | Rows: {len(frame)} | Columns: {len(frame.columns)}"
        )
        preview = frame.head(self.PREVIEW_ROWS)
        self.preview.setRowCount(len(preview))
        self.preview.setColumnCount(len(preview.columns))
        self.preview.setHorizontalHeaderLabels([str(value) for value in preview.columns])
        for row_index, row in enumerate(preview.itertuples(index=False, name=None)):
            for column_index, value in enumerate(row):
                self.preview.setItem(row_index, column_index, QTableWidgetItem("" if value is None else str(value)))
        metadata = result.column_metadata
        self.columns.setRowCount(len(metadata))
        for row_index, item in enumerate(metadata):
            values = (item["name"], item["dtype"], item["non_null"], item["missing"])
            for column_index, value in enumerate(values):
                self.columns.setItem(row_index, column_index, QTableWidgetItem(str(value)))
        self.set_mapping(mapping or {})

    def set_mapping(self, mapping: dict[str, dict[str, Any]]) -> None:
        self.mapping.setRowCount(len(mapping))
        self._detections = {}
        for row, (name, config) in enumerate(mapping.items()):
            self.mapping.setItem(row, 0, QTableWidgetItem(name))
            role = self._combo(ROLE_LABELS, config.get("role", ColumnRole.OTHER.value))
            primary = QCheckBox()
            primary.setChecked(bool(config.get("primary_identifier")))
            identifier = self._combo(IDENTIFIER_LABELS, config.get("identifier_type", IdentifierType.UNKNOWN.value))
            condition = QLineEdit(str(config.get("condition", "")))
            replicate = QLineEdit(str(config.get("replicate", "")))
            replicate.setMaximumWidth(100)
            quantification = self._combo(QUANTIFICATION_LABELS, config.get("quantification_type", QuantificationType.UNKNOWN.value))
            detection = config.get("detection") or {}
            self._detections[name] = detection
            confidence = round(float(detection.get("confidence", 0)) * 100)
            detection_text = f"{confidence}% — {detection.get('reason', 'Sem sugestão')}"
            if detection.get("contains_multiple"):
                detection_text += " Algumas células contêm múltiplos identificadores."
            self.mapping.setCellWidget(row, 1, role)
            self.mapping.setCellWidget(row, 2, primary)
            self.mapping.setCellWidget(row, 3, identifier)
            self.mapping.setCellWidget(row, 4, condition)
            self.mapping.setCellWidget(row, 5, replicate)
            self.mapping.setCellWidget(row, 6, quantification)
            self.mapping.setItem(row, 7, QTableWidgetItem(detection_text))
            role.currentIndexChanged.connect(lambda _index, current=row: self._update_row(current))
            primary.toggled.connect(lambda checked, current=row: self._exclusive_primary(current, checked))
            self._update_row(row)

    def _update_row(self, row: int) -> None:
        role = self.mapping.cellWidget(row, 1).currentData()
        is_identifier = role == ColumnRole.IDENTIFIER.value
        is_quantification = role == ColumnRole.QUANTIFICATION.value
        self.mapping.cellWidget(row, 2).setEnabled(is_identifier)
        self.mapping.cellWidget(row, 3).setEnabled(is_identifier)
        for column in (4, 5, 6):
            self.mapping.cellWidget(row, column).setEnabled(is_quantification)
        if not is_identifier:
            self.mapping.cellWidget(row, 2).setChecked(False)

    def _exclusive_primary(self, selected_row: int, checked: bool) -> None:
        if checked:
            for row in range(self.mapping.rowCount()):
                if row != selected_row:
                    self.mapping.cellWidget(row, 2).setChecked(False)

    def mapping_values(self) -> dict[str, dict[str, Any]]:
        values: dict[str, dict[str, Any]] = {}
        for row in range(self.mapping.rowCount()):
            name = self.mapping.item(row, 0).text()
            values[name] = {
                "role": self.mapping.cellWidget(row, 1).currentData(),
                "identifier_type": self.mapping.cellWidget(row, 3).currentData(),
                "primary_identifier": self.mapping.cellWidget(row, 2).isChecked(),
                "condition": self.mapping.cellWidget(row, 4).text().strip(),
                "replicate": self.mapping.cellWidget(row, 5).text().strip(),
                "quantification_type": self.mapping.cellWidget(row, 6).currentData(),
                "detection": self._detections.get(name),
            }
        return values

    def _apply_suggestions(self) -> None:
        for row in range(self.mapping.rowCount()):
            name = self.mapping.item(row, 0).text()
            detection = self._detections.get(name, {})
            detected_type = detection.get("detected_type", IdentifierType.UNKNOWN.value)
            if detected_type != IdentifierType.UNKNOWN.value and float(detection.get("confidence", 0)) >= 0.55:
                role = self.mapping.cellWidget(row, 1)
                role.setCurrentIndex(role.findData(ColumnRole.IDENTIFIER.value))
                identifier = self.mapping.cellWidget(row, 3)
                identifier.setCurrentIndex(identifier.findData(detected_type))

    def show_validation(self, validation: MappingValidation) -> None:
        if validation.valid:
            self.status.setText("Configuration: Valid")
            values = self.mapping_values()
            primary = next(((name, item) for name, item in values.items() if item["primary_identifier"]), None)
            counts = Counter(item["condition"] for item in values.values() if item["role"] == ColumnRole.QUANTIFICATION.value and item["condition"])
            lines: list[str] = []
            if primary:
                identifier_label = IDENTIFIER_LABELS[IdentifierType(primary[1]["identifier_type"])]
                lines.extend(["Identificador principal", f"{primary[0]} — {identifier_label}"])
            if counts:
                lines.append("\nCondições")
                lines.extend(f"{condition}: {count} coluna(s) quantitativa(s)" for condition, count in sorted(counts.items()))
            self.design_summary.setText("\n".join(lines))
        else:
            self.status.setText("Configuration: Invalid\n" + "\n".join(f"• {item}" for item in validation.errors))
            self.design_summary.clear()
        if validation.warnings:
            self.status.setText(self.status.text() + "\nAvisos:\n" + "\n".join(f"• {item}" for item in validation.warnings))
