"""Open only the requested persisted run and retain its target identity."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from PySide6.QtWidgets import QLabel, QLineEdit, QTableWidget, QTabWidget

from ..core.go_analysis import GOOutputs
from ..core.presence_analysis import PresenceOutputs
from .cross_module_actions import IDENTIFIER_COLUMNS


def _historical_view(module_id, outputs):
    if module_id not in {"go", "presence_absence"}:
        return outputs
    root = Path(outputs["run_root"])
    metadata = outputs["metadata"]
    def frame(relative):
        path = root / relative
        return pd.read_csv(path) if path.is_file() else pd.DataFrame()
    if module_id == "presence_absence":
        return PresenceOutputs(frame("tables/presence_matrix.csv"),
            frame("tables/condition_detection.csv"), frame("tables/classification.csv"),
            frame("tables/condition_summary.csv"), metadata,
            tuple(sorted((root / "graphs").glob("*.png"))))
    return GOOutputs(frame("annotation/go_annotations.csv"),
        frame("annotation/unannotated_proteins.csv"), frame("summary.csv"), metadata,
        tuple(sorted(root.glob("**/*.csv"))), tuple(sorted((root / "graphs").glob("*.png"))))


def _target_hint(page):
    hint = getattr(page, "cross_module_target_hint", None)
    if hint is None:
        hint = QLabel()
        hint.setWordWrap(True)
        layout = page.layout()
        if hasattr(layout, "insertWidget"):
            layout.insertWidget(0, hint)
        else:
            layout.addWidget(hint)
        page.cross_module_target_hint = hint
    return hint


def open_persisted_target(project, adapter, page, run_id, target_identity, *, source_row=""):
    """Load the exact run, then focus a safe target or display an explicit fallback."""
    if run_id not in adapter.list_runs(project):
        raise FileNotFoundError(f"Run {run_id} is not available for {adapter.module_id}.")
    outputs = adapter.load_run(project, run_id)
    metadata = outputs["metadata"] if isinstance(outputs, dict) else outputs.metadata
    if str(metadata.get("run_id")) != str(run_id):
        raise RuntimeError("Loaded artifacts do not belong to the requested run.")
    view = _historical_view(adapter.module_id, outputs)
    display = getattr(page, "_show_statistics", None) if adapter.module_id == "differential" else getattr(page, "show_outputs", None)
    if display is None:
        raise RuntimeError("Analysis page cannot display persisted outputs.")
    display(view)
    history = getattr(page, "history", None) or getattr(page, "stat_history", None)
    if history is not None:
        index = history.findData(run_id)
        if index < 0:
            index = history.findText(run_id)
        if index >= 0:
            history.blockSignals(True)
            history.setCurrentIndex(index)
            history.blockSignals(False)
    loaded = getattr(page, "statistics", None) if adapter.module_id == "differential" else getattr(page, "outputs", None)
    loaded_metadata = loaded.get("metadata", {}) if isinstance(loaded, dict) else getattr(loaded, "metadata", {})
    if str(loaded_metadata.get("run_id")) != str(run_id):
        raise RuntimeError("Analysis page did not retain the requested run.")
    target = str(target_identity).strip()
    if not target:
        raise ValueError("Target identity is empty.")
    if adapter.module_id == "differential":
        page.stage_tabs.setCurrentIndex(2)
        page.result_tabs.setCurrentIndex(1)
    elif adapter.module_id == "proteomics_qc":
        page.tabs.setCurrentIndex(2)
        page.table_tabs["Detection"].setCurrentIndex(1)
    elif adapter.module_id in {"presence_absence", "go"}:
        tabs = page.findChildren(QTabWidget)
        if tabs:
            tabs[0].setCurrentIndex(1)
    focused = False
    if source_row:
        for table in page.findChildren(QTableWidget):
            if not table.isVisible():
                continue
            columns = [column for column in range(table.columnCount())
                if table.horizontalHeaderItem(column) is not None
                and table.horizontalHeaderItem(column).text() in {"source_row", "Source Row", "Source row"}]
            for column in columns:
                matches = [row for row in range(table.rowCount())
                    if table.item(row, column) and table.item(row, column).text() == str(source_row)]
                if len(matches) == 1:
                    table.selectRow(matches[0])
                    table.scrollToItem(table.item(matches[0], column))
                    focused = True
                    break
            if focused:
                break
    for name in (() if focused else ("pathway", "protein", "gene_choice", "entity_choice", "complex_choice")):
        combo = getattr(page, name, None)
        if combo is not None and hasattr(combo, "findData"):
            index = combo.findData(target)
            if index < 0:
                index = combo.findText(target)
            if index >= 0:
                combo.setCurrentIndex(index)
                focused = True
                break
    if not focused:
        for table in page.findChildren(QTableWidget):
            if not table.isVisible():
                continue
            id_columns = [column for column in range(table.columnCount())
                if table.horizontalHeaderItem(column) is not None
                and table.horizontalHeaderItem(column).text() in IDENTIFIER_COLUMNS]
            for row in range(table.rowCount()):
                if any(table.item(row, column) and table.item(row, column).text() == target
                       for column in id_columns):
                    table.selectRow(row)
                    table.scrollToItem(table.item(row, 0))
                    focused = True
                    break
            if focused:
                break
    if not focused:
        search = getattr(page, "search", None)
        if isinstance(search, QLineEdit):
            search.setText(target)
    _target_hint(page).setText(
        f"Loaded run: {run_id} | Target: {target}"
        + (" | Selected in results." if focused else " | Locate in persisted results."))
    return focused
