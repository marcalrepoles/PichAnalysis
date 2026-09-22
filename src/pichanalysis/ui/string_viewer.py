"""Qt presentation of persisted STRING run nodes and edges.

Geometry, dragging and display filters never modify scientific CSV artifacts.
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QPainter, QPen, QPixmap
from PySide6.QtSvg import QSvgGenerator
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QGraphicsEllipseItem, QGraphicsItem,
    QGraphicsLineItem, QGraphicsScene, QGraphicsTextItem, QGraphicsView,
    QHBoxLayout, QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget,
)

VIEWER_NODE_LIMIT = 1500
HOP_COLORS = {0: QColor("#2563eb"), 1: QColor("#16a34a"), 2: QColor("#f59e0b")}


def _value(row, key, default=""):
    value = row.get(key, default)
    return default if pd.isna(value) else value


class NodeItem(QGraphicsEllipseItem):
    def __init__(self, node, click, moved):
        hop = int(_value(node, "hop_level", 0))
        radius = 13 if hop == 0 else 10
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self.node = node; self.hop_level = hop; self.click = click; self.moved = moved
        self.setBrush(QBrush(HOP_COLORS.get(hop, QColor("#64748b"))))
        self.setPen(QPen(QColor("#111827"), 3 if bool(_value(node, "is_hub", False)) else 1.5))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, True)
        self.setZValue(2)
        self.label = QGraphicsTextItem(str(_value(node, "gene_symbol") or _value(node, "preferred_name") or _value(node, "string_protein_id")))
        self.label.setDefaultTextColor(QColor("#111827")); self.label.setParentItem(self)
        self.label.setPos(radius + 2, -radius)

    def mousePressEvent(self, event):
        self.click(str(_value(self.node, "string_protein_id")))
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        super().mouseMoveEvent(event)
        self.moved()


class EdgeItem(QGraphicsLineItem):
    def __init__(self, edge, a: NodeItem, b: NodeItem, click):
        super().__init__(); self.edge = edge; self.a = a; self.b = b; self.click = click
        score = int(_value(edge, "combined_score", 0))
        self.setPen(QPen(QColor(100, 116, 139, 180), 1 + score / 400))
        self.setZValue(0); self.update_position()

    def update_position(self): self.setLine(self.a.pos().x(), self.a.pos().y(), self.b.pos().x(), self.b.pos().y())
    def mousePressEvent(self, event):
        self.click(str(_value(self.edge, "protein_a")), str(_value(self.edge, "protein_b")))
        super().mousePressEvent(event)


class NetworkGraphicsView(QGraphicsView):
    def __init__(self, scene):
        super().__init__(scene)
        self.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)

    def wheelEvent(self, event):
        self.scale(1.15 if event.angleDelta().y() > 0 else 1 / 1.15,
                   1.15 if event.angleDelta().y() > 0 else 1 / 1.15)


class StringNetworkViewer(QWidget):
    node_selected = Signal(str)
    edge_selected = Signal(str, str)

    def __init__(self):
        super().__init__(); self.outputs = None; self.node_items = {}; self.edge_items = []
        self.initial_positions = {}; self.visible_ids = set(); self.highlighted = set()
        self.title = QLabel("No STRING run loaded."); self.legend = QLabel()
        self.message = QLabel(); self.message.setWordWrap(True)
        self.layers = QComboBox()
        for label, data in (("All nodes", "all"), ("Seeds", "seeds"), ("Degree 1", "degree1"), ("Degree 2", "degree2"), ("Hubs only", "hubs")):
            self.layers.addItem(label, data)
        self.component = QComboBox(); self.component.addItem("All components", None)
        self.search = QLineEdit(); self.search.setPlaceholderText("Search Gene, UniProt, STRING ID or preferred name")
        self.show_labels = QCheckBox("Show labels"); self.show_labels.setChecked(True)
        self.highlight_common = QPushButton("Highlight common direct neighbors")
        self.seed_focus = QComboBox(); self.seed_focus.addItem("All seeds", None)
        self.fit_button = QPushButton("Fit to Window"); self.zoom_in = QPushButton("Zoom In")
        self.zoom_out = QPushButton("Zoom Out"); self.reset_view = QPushButton("Reset View")
        self.reset_layout_button = QPushButton("Reset Layout")
        self.export_button = QPushButton("Export network image...")
        self.scene = QGraphicsScene(self); self.view = NetworkGraphicsView(self.scene)
        controls = QHBoxLayout()
        for widget in (self.layers, self.component, self.seed_focus, self.search, self.show_labels): controls.addWidget(widget)
        actions = QHBoxLayout()
        for button in (self.highlight_common, self.fit_button, self.zoom_in, self.zoom_out, self.reset_view, self.reset_layout_button, self.export_button): actions.addWidget(button)
        layout = QVBoxLayout(self); layout.addWidget(self.title); layout.addWidget(self.legend); layout.addLayout(controls)
        layout.addWidget(self.message); layout.addWidget(self.view, 1); layout.addLayout(actions)
        self.layers.currentIndexChanged.connect(self.render); self.component.currentIndexChanged.connect(self.render)
        self.seed_focus.currentIndexChanged.connect(self.render); self.search.returnPressed.connect(self.search_node)
        self.show_labels.toggled.connect(self._toggle_labels); self.highlight_common.clicked.connect(self.highlight_common_nodes)
        self.fit_button.clicked.connect(self.fit_to_window); self.zoom_in.clicked.connect(lambda:self.view.scale(1.2,1.2))
        self.zoom_out.clicked.connect(lambda:self.view.scale(1/1.2,1/1.2))
        self.reset_view.clicked.connect(self.reset_view_transform); self.reset_layout_button.clicked.connect(self.reset_layout)
        self.export_button.clicked.connect(self.export_image_dialog)

    def load_outputs(self, outputs):
        self.outputs = outputs; m = outputs["metadata"]
        network = "Functional association" if m.get("network_type") == "functional" else "Physical"
        mode = " | Direct-neighbor mode: Strict common" if m.get("degree1_selection_mode") == "strict_common" else ""
        self.title.setText(f"STRING {network} network | Threshold ≥ {m.get('combined_score_threshold')} | Maximum hop: {m.get('max_hop')}{mode}")
        self.legend.setText(f"Seed: blue large node | Degree 1: green | Degree 2: amber | Hub: bold border | Edge: STRING association | Network type: {network} | Score threshold: {m.get('combined_score_threshold')}")
        self.component.blockSignals(True); self.component.clear(); self.component.addItem("All components", None)
        metrics = outputs["tables"]["node_metrics"]
        for component in sorted(metrics.component_id.dropna().unique()): self.component.addItem(f"Component {int(component)}", int(component))
        self.component.blockSignals(False)
        self.seed_focus.blockSignals(True); self.seed_focus.clear(); self.seed_focus.addItem("All seeds",None)
        for row in outputs["tables"]["seeds"].to_dict("records"):
            self.seed_focus.addItem(str(_value(row,"gene_symbol") or _value(row,"preferred_name") or _value(row,"string_protein_id")),str(row["string_protein_id"]))
        self.seed_focus.blockSignals(False); self.layers.setCurrentIndex(0); self.render()

    def _filtered_ids(self):
        if not self.outputs:return set()
        nodes=self.outputs["tables"]["expanded_nodes"];metrics=self.outputs["tables"]["node_metrics"]
        ids=set(nodes.string_protein_id.astype(str)); mode=self.layers.currentData()
        if mode in {"seeds","degree1","degree2"}:
            hop={"seeds":0,"degree1":1,"degree2":2}[mode]
            ids &= set(nodes.loc[nodes.hop_level==hop,"string_protein_id"].astype(str))
        if mode=="hubs":ids &= set(self.outputs["tables"]["hubs"].string_protein_id.astype(str))
        component=self.component.currentData()
        if component is not None:ids &= set(metrics.loc[metrics.component_id==component,"string_protein_id"].astype(str))
        focus=self.seed_focus.currentData()
        if focus:
            neighbors={focus}
            for edge in self.outputs["tables"]["expanded_edges"].to_dict("records"):
                if focus==str(edge["protein_a"]):neighbors.add(str(edge["protein_b"]))
                if focus==str(edge["protein_b"]):neighbors.add(str(edge["protein_a"]))
            ids &= neighbors
        return ids

    def render(self):
        self.scene.clear(); self.node_items={};self.edge_items=[];self.initial_positions={}
        ids=self._filtered_ids();self.visible_ids=ids
        if not self.outputs:return
        if len(ids)>VIEWER_NODE_LIMIT:
            self.message.setText(f"This network contains {len(ids)} nodes. Apply display filters before rendering the interactive network.")
            return
        nodes=self.outputs["tables"]["expanded_nodes"]
        if not ids:self.message.setText("No nodes match the current display filters.");return
        edges=self.outputs["tables"]["expanded_edges"]
        self.message.setText("No qualifying edges were found at this threshold." if edges.empty else "Display filters affect only this viewer; scientific tables and metrics remain unchanged.")
        metrics=self.outputs["tables"]["node_metrics"].set_index("string_protein_id")
        hubs=set(self.outputs["tables"]["hubs"].string_protein_id.astype(str))
        visible=nodes[nodes.string_protein_id.astype(str).isin(ids)].copy()
        for component_index,(component,group) in enumerate(visible.groupby(visible.string_protein_id.map(lambda x:metrics.loc[x,"component_id"] if x in metrics.index else 0),sort=True)):
            center_x=(component_index%4)*520;center_y=(component_index//4)*520
            for hop,part in group.groupby("hop_level",sort=True):
                part=part.sort_values("string_protein_id");count=len(part);radius=60+int(hop)*115
                for i,row in enumerate(part.to_dict("records")):
                    node_id=str(row["string_protein_id"]);row["is_hub"]=node_id in hubs
                    angle=2*math.pi*i/max(count,1)
                    position=QPointF(center_x+radius*math.cos(angle),center_y+radius*math.sin(angle))
                    item=NodeItem(row,self.node_selected.emit,self.refresh_edges);item.setPos(position);item.label.setVisible(self.show_labels.isChecked())
                    self.scene.addItem(item);self.node_items[node_id]=item;self.initial_positions[node_id]=position
        for edge in edges.to_dict("records"):
            a,b=str(edge["protein_a"]),str(edge["protein_b"])
            if a in self.node_items and b in self.node_items:
                item=EdgeItem(edge,self.node_items[a],self.node_items[b],self.edge_selected.emit)
                self.scene.addItem(item);self.edge_items.append(item)
        self.fit_to_window()

    def refresh_edges(self):
        for item in self.edge_items:item.update_position()
    def _toggle_labels(self,shown):
        for item in self.node_items.values():item.label.setVisible(shown)
    def highlight_common_nodes(self):
        if not self.outputs:return
        common=set(self.outputs["tables"]["common_direct_neighbors"].string_protein_id.astype(str))
        self.highlighted=common
        for node_id,item in self.node_items.items():
            item.setPen(QPen(QColor("#a21caf") if node_id in common else QColor("#111827"),4 if node_id in common else 3 if item.node.get("is_hub") else 1.5))
    def search_node(self):
        query=self.search.text().strip().casefold()
        if not query:return None
        for node_id,item in self.node_items.items():
            row=item.node
            if any(query in str(_value(row,key)).casefold() for key in ("gene_symbol","uniprot_accession","string_protein_id","preferred_name")):
                self.view.centerOn(item);item.setSelected(True);self.node_selected.emit(node_id);return node_id
        self.message.setText("No matching node is currently displayed. Adjust display filters to search another node.")
        return None
    def fit_to_window(self):
        if self.scene.items():self.view.fitInView(self.scene.itemsBoundingRect().adjusted(-40,-40,40,40),Qt.AspectRatioMode.KeepAspectRatio)
    def reset_view_transform(self):self.view.resetTransform();self.fit_to_window()
    def reset_layout(self):
        for node_id,pos in self.initial_positions.items():self.node_items[node_id].setPos(pos)
        self.refresh_edges();self.fit_to_window()
    def export_image(self,path):
        if not self.scene.items():raise ValueError("No network is rendered for export.")
        path=Path(path);bounds=self.scene.itemsBoundingRect().adjusted(-20,-20,20,20)
        if path.suffix.lower()==".svg":
            generator=QSvgGenerator();generator.setFileName(str(path));generator.setSize(bounds.size().toSize());generator.setViewBox(bounds.toRect())
            painter=QPainter(generator);self.scene.render(painter,QRectF(generator.viewBox()),bounds);painter.end()
        else:
            pixmap=QPixmap(bounds.size().toSize());pixmap.fill(Qt.GlobalColor.white)
            painter=QPainter(pixmap);self.scene.render(painter,QRectF(pixmap.rect()),bounds);painter.end();pixmap.save(str(path),"PNG")
    def export_image_dialog(self):
        path,_=QFileDialog.getSaveFileName(self,"Export network image","string_network.png","PNG (*.png);;SVG (*.svg)")
        if path:self.export_image(path)
