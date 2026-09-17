from __future__ import annotations

from PySide6.QtWidgets import QLabel, QPushButton, QVBoxLayout, QWidget

from ..core.project import Project


class ProjectPage(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.new_button = QPushButton("Novo projeto")
        self.open_button = QPushButton("Abrir projeto existente")
        self.folder_button = QPushButton("Abrir pasta do projeto")
        self.folder_button.setEnabled(False)
        self.details = QLabel("Nenhum projeto aberto.")
        self.details.setWordWrap(True)
        layout = QVBoxLayout(self)
        layout.addWidget(self.new_button)
        layout.addWidget(self.open_button)
        layout.addWidget(self.folder_button)
        layout.addWidget(self.details)
        layout.addStretch()

    def show_project(self, project: Project) -> None:
        imported = project.config.get("input", {}).get("original_file") or "Nenhum"
        self.details.setText(
            f"Projeto: {project.name}\nPasta: {project.root}\nArquivo importado: {imported}"
        )
        self.folder_button.setEnabled(True)

