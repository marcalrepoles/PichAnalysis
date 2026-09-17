import pandas as pd

from pichanalysis.core.importer import import_into_project, read_table
from pichanalysis.core.project import create_project


def test_read_csv(tmp_path):
    path = tmp_path / "sample.csv"
    path.write_text("protein,value\nP1,1\nP2,\n", encoding="utf-8")
    frame, sheet = read_table(path)
    assert list(frame.columns) == ["protein", "value"]
    assert len(frame) == 2
    assert sheet is None


def test_read_tsv(tmp_path):
    path = tmp_path / "sample.tsv"
    path.write_text("protein\tvalue\nP1\t1\n", encoding="utf-8")
    frame, _ = read_table(path)
    assert frame.loc[0, "protein"] == "P1"


def test_read_simple_xlsx(tmp_path):
    path = tmp_path / "sample.xlsx"
    pd.DataFrame({"protein": ["P1"], "value": [1]}).to_excel(path, index=False)
    frame, sheet = read_table(path)
    assert sheet == "Sheet1"
    assert frame.shape == (1, 2)


def test_import_records_relative_files_and_metadata(tmp_path):
    project = create_project(tmp_path, "Import")
    source = tmp_path / "sample.csv"
    source.write_text("protein,value\nP1,1\n", encoding="utf-8")
    result = import_into_project(project, source)
    assert result.original_path.exists()
    assert result.processed_path.exists()
    assert project.config["input"]["original_file"] == "input/original/sample.csv"
    assert project.config["input"]["column_metadata"][0]["role"] is None
