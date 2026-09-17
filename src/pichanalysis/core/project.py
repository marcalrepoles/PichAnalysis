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
            raise ProjectError(f"Não foi possível salvar project.json: {error}") from error

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
            raise ProjectError("O arquivo selecionado não existe.")
        destination = self.reserve_original(source)
        try:
            shutil.copy2(source, destination)
        except OSError as error:
            raise ProjectError(f"Não foi possível copiar o arquivo original: {error}") from error
        return destination


def create_project(parent: Path, name: str) -> Project:
    parent = Path(parent).expanduser()
    clean_name = name.strip()
    if not clean_name or clean_name in {".", ".."}:
        raise ProjectError("Informe um nome de projeto válido.")
    if re.search(r'[<>:"/\\|?*]', clean_name):
        raise ProjectError("O nome do projeto contém caracteres não permitidos.")
    root = parent / clean_name
    if root.exists():
        raise ProjectError("Já existe um arquivo ou pasta com esse nome.")
    created_at = _now()
    try:
        root.mkdir(parents=True)
        for directory in PROJECT_DIRECTORIES:
            (root / directory).mkdir(parents=True)
    except OSError as error:
        raise ProjectError(f"Não foi possível criar a estrutura do projeto: {error}") from error
    project = Project(
        root=root.resolve(),
        config={
            "schema_version": 1,
            "application": "PichAnalysis",
            "project_name": clean_name,
            "created_at": created_at,
            "updated_at": created_at,
            "project_root": ".",
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
                "errors": ["Nenhum identificador principal selecionado."],
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
        raise ProjectError("A pasta não contém project.json.") from error
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ProjectError(f"project.json inválido ou ilegível: {error}") from error
    required = {"schema_version", "application", "project_name", "input"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise ProjectError("project.json não contém os campos obrigatórios.")
    if config["application"] != "PichAnalysis" or config["schema_version"] != 1:
        raise ProjectError("A pasta não é um projeto PichAnalysis compatível.")
    if not all((root / directory).is_dir() for directory in PROJECT_DIRECTORIES):
        raise ProjectError("A estrutura de pastas do projeto está incompleta.")
    return Project(root=root, config=config)
