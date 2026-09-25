"""Shared, selection-driven Explore action for persisted analysis pages."""
from __future__ import annotations

from PySide6.QtWidgets import QLabel, QPushButton, QTableWidget


IDENTIFIER_COLUMNS = (
    "feature_id", "Feature ID", "Feature Id", "source_row", "Source row", "entity_key", "reactome_entity_key", "canonical_uniprot",
    "uniprot_accession", "ncbi_gene_id", "kegg_gene_id", "string_protein_id",
    "complex_id", "GO_ID", "go_id", "Reactome_ID", "pathway_id",
    "original_id", "original_identifier", "input_id", "gene_symbol",
    "UniProt", "Gene ID",
)


def selected_identifier(page):
    """Return a selected persisted identifier; never infer one from current input."""
    candidates = [table for table in page.findChildren(QTableWidget)
                  if table.isVisible() and table.currentRow() >= 0
                  and table.selectedItems()]
    if len(candidates) != 1:
        return None
    table = candidates[0]
    headers = {table.horizontalHeaderItem(column).text(): column
               for column in range(table.columnCount())
               if table.horizontalHeaderItem(column) is not None}
    for name in IDENTIFIER_COLUMNS:
        column = headers.get(name)
        if column is not None:
            cell = table.item(table.currentRow(), column)
            if cell and cell.text().strip():
                return cell.text().strip()
    return None


def attach_explore_action(page, module_id, callback):
    """Attach the same action to a page, using only its loaded run and selection."""
    button = QPushButton("Explore across analyses...")
    hint = QLabel("")
    hint.setWordWrap(True)
    page.layout().addWidget(button)
    page.layout().addWidget(hint)

    def explore():
        outputs = getattr(page, "statistics", None) if module_id == "differential" else getattr(page, "outputs", None)
        metadata = outputs.get("metadata", {}) if isinstance(outputs, dict) else getattr(outputs, "metadata", {})
        run_id = str(metadata.get("run_id", ""))
        identifier = selected_identifier(page)
        if not run_id or not identifier:
            hint.setText("Select one result row in the loaded run before exploring across analyses.")
            return
        hint.clear()
        callback(module_id, run_id, identifier)

    button.clicked.connect(explore)
    page.cross_module_explore_button = button
    page.cross_module_explore_hint = hint
    return button, hint

def attach_biological_context_action(page, module_id, callback):
    """One shared action over the selected row of a loaded, persisted run."""
    button = QPushButton("Biological context...")
    button.setEnabled(False)
    hint = QLabel("")
    hint.setWordWrap(True)
    page.layout().addWidget(button)
    page.layout().addWidget(hint)

    def selection():
        identifier = selected_identifier(page)
        rows = [table for table in page.findChildren(QTableWidget)
                if table.isVisible() and table.currentRow() >= 0 and table.selectedItems()]
        if len(rows) == 1:
            table = rows[0]
            headers = {table.horizontalHeaderItem(column).text(): column
                for column in range(table.columnCount()) if table.horizontalHeaderItem(column)}
            for name in ("feature_id", "uniprot_accession", "canonical_uniprot",
                         "reactome_entity_key", "string_protein_id", "kegg_gene_id",
                         "gene_symbol", "original_id", "original_identifier"):
                if name in headers and table.item(table.currentRow(), headers[name]):
                    value = table.item(table.currentRow(), headers[name]).text().strip()
                    if value:
                        identifier = value
                        break
        if not identifier and module_id in {"kegg", "reactome"}:
            protein = getattr(page, "protein", None)
            identifier = str(protein.currentData() or "") if protein else ""
        outputs = getattr(page, "outputs", None)
        metadata = outputs.get("metadata", {}) if isinstance(outputs, dict) else getattr(outputs, "metadata", {})
        return str(metadata.get("run_id", "")), identifier

    def update():
        run_id, identifier = selection()
        button.setEnabled(bool(run_id and identifier))

    def explore():
        run_id, identifier = selection()
        if not run_id or not identifier:
            hint.setText("Select one entity in a loaded run first.")
            return
        rows = [table for table in page.findChildren(QTableWidget)
                if table.isVisible() and table.currentRow() >= 0 and table.selectedItems()]
        source_row = None
        if len(rows) == 1:
            table = rows[0]
            for column in range(table.columnCount()):
                header = table.horizontalHeaderItem(column)
                if header and header.text() in {"source_row", "Source row", "Source Row"}:
                    item = table.item(table.currentRow(), column)
                    if item and item.text().isdigit():
                        source_row = int(item.text())
                    break
        hint.clear()
        callback(module_id, run_id, identifier, source_row)

    for table in page.findChildren(QTableWidget):
        table.itemSelectionChanged.connect(update)
    protein = getattr(page, "protein", None)
    if protein is not None and module_id in {"kegg", "reactome"}:
        protein.currentIndexChanged.connect(update)
    button.clicked.connect(explore)
    page.biological_context_button = button
    page.biological_context_hint = hint
    return button
