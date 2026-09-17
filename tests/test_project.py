import json

import pytest

from pichanalysis.core.project import PROJECT_DIRECTORIES, ProjectError, create_project, open_project


def test_create_project_structure_and_json(tmp_path):
    project = create_project(tmp_path, "Experimento")
    assert all((project.root / path).is_dir() for path in PROJECT_DIRECTORIES)
    config = json.loads((project.root / "project.json").read_text(encoding="utf-8"))
    assert config["application"] == "PichAnalysis"
    assert config["schema_version"] == 1
    assert config["project_root"] == "."


def test_open_valid_project(tmp_path):
    created = create_project(tmp_path, "Valido")
    opened = open_project(created.root)
    assert opened.name == "Valido"


def test_reject_non_project_directory(tmp_path):
    with pytest.raises(ProjectError, match="project.json"):
        open_project(tmp_path)


def test_original_copy_never_overwrites(tmp_path):
    project = create_project(tmp_path, "Copias")
    source = tmp_path / "dados.csv"
    source.write_text("a\n1\n", encoding="utf-8")
    first = project.copy_original(source)
    source.write_text("a\n2\n", encoding="utf-8")
    second = project.copy_original(source)
    assert first.name == "dados.csv"
    assert second.name == "dados_2.csv"
    assert first.read_text(encoding="utf-8") == "a\n1\n"

