"""Simple readable, on-demand context figures; official images remain untouched."""
from __future__ import annotations

import math

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen

from .biological_context import GREEN, RED, GRAY, frozen_direction
from . import differential_analysis

MAX_NODES = 40


def render_context(project, source, item, differential_run=None):
    key = {"kegg": "kegg_gene_id", "reactome": "reactome_entity_key",
           "string": "string_protein_id"}[item.module]
    rows = item.members.fillna("").astype(str).to_dict("records")
    unique = {}
    for row in rows:
        identity = str(row.get(key, "")).strip()
        if identity:
            unique.setdefault(identity, row)
    selected = unique.pop(item.selected_entity, None)
    ordered = ([(item.selected_entity, selected)] if selected is not None else []) + sorted(unique.items())
    hidden = max(0, len(ordered) - MAX_NODES)
    ordered = ordered[:MAX_NODES]
    columns = min(5, max(1, math.ceil(math.sqrt(len(ordered)))))
    lines = max(1, math.ceil(len(ordered) / columns))
    width = 1250
    height = max(520, 255 + lines * 126)
    image = QImage(width, height, QImage.Format_ARGB32)
    image.fill(QColor("white"))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.setPen(QColor("#222222"))
    painter.setFont(QFont("Arial", 19, QFont.Bold))
    painter.drawText(QRectF(35, 22, width - 70, 38), Qt.AlignLeft, item.title)
    painter.setFont(QFont("Arial", 11))
    painter.drawText(QRectF(35, 64, width - 70, 28), Qt.AlignLeft,
        f"{item.module.upper()} {item.item_id}  |  Run {item.run_id}  |  Snapshot {item.snapshot or 'not recorded'}")
    painter.drawText(QRectF(35, 94, width - 70, 28), Qt.AlignLeft,
        f"Directional source: {differential_run or 'None'} (Condition A - Condition B)")
    painter.drawText(QRectF(35, 122, width - 70, 28), Qt.AlignLeft,
        "Status of detected members in this sample/run; not pathway activity.")
    if hidden:
        painter.drawText(QRectF(35, 150, width - 70, 28), Qt.AlignLeft,
            f"Showing {len(ordered)} of {len(ordered) + hidden} detected members for readability.")
    differential_results = (differential_analysis.load_run(project, differential_run)["tables"]["all_results"]
        if differential_run else None)
    positions = {}
    cell_width = (width - 80) / columns
    for index, (identity, _) in enumerate(ordered):
        col, row = index % columns, index // columns
        positions[identity] = (40 + cell_width * (col + .5), 245 + row * 126)
    if item.module == "string":
        painter.setPen(QPen(QColor("#a7adb1"), 1.6))
        for edge in item.edges.to_dict("records"):
            left, right = str(edge.get("protein_a", "")), str(edge.get("protein_b", ""))
            if left in positions and right in positions:
                painter.drawLine(*map(int, positions[left] + positions[right]))
    for identity, row in ordered:
        x, y = positions[identity]
        direction = frozen_direction(project, differential_run, source, row, differential_results) if differential_run else "none"
        color = QColor({"up": GREEN, "down": RED}.get(direction, GRAY))
        painter.setBrush(color)
        painter.setPen(QPen(QColor("#202020"), 4 if identity == item.selected_entity else 1))
        painter.drawEllipse(QRectF(x - 23, y - 23, 46, 46))
        label = next((str(row.get(name, "")).strip() for name in
            ("Gene_symbol", "gene_symbol", "preferred_name", "UniProt", "uniprot_accession")
            if str(row.get(name, "")).strip()), identity)
        label = label[:22]
        painter.setPen(QColor("#202020"))
        painter.setFont(QFont("Arial", 10, QFont.Bold if identity == item.selected_entity else QFont.Normal))
        painter.drawText(QRectF(x - cell_width / 2 + 5, y + 27, cell_width - 10, 42),
            Qt.AlignHCenter | Qt.TextWordWrap, label)
    painter.setFont(QFont("Arial", 10))
    for offset, (label, color) in enumerate((("Up", GREEN), ("Down", RED),
                                              ("No change / no direction", GRAY))):
        x = 40 + offset * 200
        painter.setBrush(QColor(color))
        painter.setPen(QColor("#333333"))
        painter.drawEllipse(QRectF(x, height - 45, 14, 14))
        painter.drawText(QRectF(x + 22, height - 50, 190, 28), Qt.AlignLeft, label)
    painter.drawText(QRectF(width - 340, height - 50, 300, 28), Qt.AlignRight,
        "Thick outline = selected entity")
    painter.end()
    return image