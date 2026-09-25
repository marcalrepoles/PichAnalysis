from __future__ import annotations

import json
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class ProjectError(Exception):
    """Raised when a project cannot be created, opened, or updated."""


PROJECT_DIRECTORIES = (
    "input/original",
    "input/processed",
    "mapping/tables",
    "mapping/raw",
    "annotations",
    "analyses",
    "reports",
    "logs",
    "scripts/runs",
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class Project:
    root: Path
    config: dict[str, Any]

    @property
    def name(self) -> str:
        return str(self.config["project_name"])

    @property
    def config_path(self) -> Path:
        return self.root / "project.json"

    def save(self) -> None:
        self.config["updated_at"] = _now()
        temporary = self.config_path.with_suffix(".json.tmp")
        try:
            temporary.write_text(
                json.dumps(self.config, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
            temporary.replace(self.config_path)
        except OSError as error:
            temporary.unlink(missing_ok=True)
            raise ProjectError(f"Could not save project.json: {error}") from error

    def relative(self, path: Path) -> str:
        return path.resolve().relative_to(self.root.resolve()).as_posix()

    def reserve_original(self, source: Path) -> Path:
        destination = self.root / "input" / "original" / source.name
        counter = 2
        while destination.exists():
            destination = destination.with_name(
                f"{source.stem}_{counter}{source.suffix}"
            )
            counter += 1
        return destination

    def copy_original(self, source: Path) -> Path:
        source = Path(source)
        if not source.is_file():
            raise ProjectError("The selected file does not exist.")
        destination = self.reserve_original(source)
        try:
            shutil.copy2(source, destination)
        except OSError as error:
            raise ProjectError(f"Could not copy the original file: {error}") from error
        return destination


def create_project(parent: Path, name: str) -> Project:
    parent = Path(parent).expanduser()
    clean_name = name.strip()
    if not clean_name or clean_name in {".", ".."}:
        raise ProjectError("Enter a valid project name.")
    if re.search(r'[<>:"/\\|?*]', clean_name):
        raise ProjectError("The project name contains forbidden characters.")
    root = parent / clean_name
    if root.exists():
        raise ProjectError("A file or folder with this name already exists.")
    created_at = _now()
    try:
        root.mkdir(parents=True)
        for directory in PROJECT_DIRECTORIES:
            (root / directory).mkdir(parents=True)
    except OSError as error:
        raise ProjectError(f"Could not create the project directory structure: {error}") from error
    project = Project(
        root=root.resolve(),
        config={
            "schema_version": 1,
            "application": "PichAnalysis",
            "project_name": clean_name,
            "created_at": created_at,
            "updated_at": created_at,
            "project_root": ".",
            "organism_name": None,
            "organism_tax_id": None,
            "input": {
                "original_file": None,
                "processed_file": None,
                "format": None,
                "sheet": None,
                "rows": None,
                "columns": None,
                "column_metadata": [],
            },
            "columns": {},
            "column_configuration": {
                "status": "not_configured",
                "errors": ["No primary identifier selected."],
                "warnings": [],
            },
        },
    )
    project.save()
    return project


def open_project(root: Path) -> Project:
    root = Path(root).expanduser().resolve()
    config_path = root / "project.json"
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise ProjectError("The folder does not contain project.json.") from error
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ProjectError(f"project.json is invalid or unreadable: {error}") from error
    required = {"schema_version", "application", "project_name", "input"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise ProjectError("project.json is missing required fields.")
    if config["application"] != "PichAnalysis" or config["schema_version"] != 1:
        raise ProjectError("The folder is not a compatible PichAnalysis project.")
    if not all((root / directory).is_dir() for directory in PROJECT_DIRECTORIES):
        raise ProjectError("The project directory structure is incomplete.")
    return Project(root=root, config=config)
