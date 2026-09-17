from pathlib import Path

from pichanalysis.core.r_runtime import RRuntime, discover_rscript


def test_discovers_latest_windows_program_files(tmp_path):
    older = tmp_path / "R" / "R-4.4.3" / "bin" / "Rscript.exe"
    newer = tmp_path / "R" / "R-4.6.0" / "bin" / "Rscript.exe"
    older.parent.mkdir(parents=True); newer.parent.mkdir(parents=True)
    older.touch(); newer.touch()
    found = discover_rscript(platform="win32", program_files=str(tmp_path), which=lambda _: None)
    assert found == str(newer)


def test_runtime_builds_arguments_without_shell_string(tmp_path):
    executable = tmp_path / "Rscript"
    executable.touch()
    runtime = RRuntime(str(executable))
    command = runtime.command(Path("folder with spaces/script.R"), "--input", "data with spaces.csv")
    assert command == (str(executable), str(Path("folder with spaces/script.R")), "--input", "data with spaces.csv")
