from __future__ import annotations

from PySide6.QtWidgets import (
    QAbstractItemView, QHeaderView, QLabel, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

from ..core.importer import ImportResult


class DataPage(QWidget):
    PREVIEW_ROWS = 200

    def __init__(self) -> None:
        super().__init__()
        self.import_button = QPushButton("Importar XLSX, CSV ou TSV")
        self.summary = QLabel("Abra ou crie um projeto para importar dados.")
        self.summary.setWordWrap(True)
        self.preview = QTableWidget()
        self.preview.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.preview.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.columns = QTableWidget(0, 4)
        self.columns.setHorizontalHeaderLabels(["Nome", "Tipo", "Preenchidos", "Ausentes"])
        self.columns.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.columns.horizontalHeader().setStretchLastSection(True)
        layout = QVBoxLayout(self)
        layout.addWidget(self.import_button)
        layout.addWidget(self.summary)
        layout.addWidget(QLabel("Pré-visualização (até 200 linhas)"))
        layout.addWidget(self.preview, 3)
        layout.addWidget(QLabel("Colunas"))
        layout.addWidget(self.columns, 2)

    def show_result(self, result: ImportResult) -> None:
        frame = result.dataframe
        sheet = f" | Planilha: {result.sheet}" if result.sheet else ""
        self.summary.setText(
            f"Arquivo: {result.original_path.name} | Formato: {result.source_format.upper()}"
            f"{sheet} | Linhas: {len(frame)} | Colunas: {len(frame.columns)}"
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

