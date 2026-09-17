from __future__ import annotations

from PySide6.QtWidgets import QLabel, QPushButton, QVBoxLayout, QWidget

from ..core.project import Project


class ProjectPage(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.new_button = QPushButton("Novo projeto")
        self.open_button = QPushButton("Open existing project")
        self.folder_button = QPushButton("Open project folder")
        self.folder_button.setEnabled(False)
        self.details = QLabel("No project is open.")
        self.details.setWordWrap(True)
        layout = QVBoxLayout(self)
        layout.addWidget(self.new_button)
        layout.addWidget(self.open_button)
        layout.addWidget(self.folder_button)
        layout.addWidget(self.details)
        layout.addStretch()

    def show_project(self, project: Project) -> None:
        imported = project.config.get("input", {}).get("original_file") or "None"
        self.details.setText(
            f"Project: {project.name}\nFolder: {project.root}\nImported file: {imported}"
        )
        self.folder_button.setEnabled(True)
