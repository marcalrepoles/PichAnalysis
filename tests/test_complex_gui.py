"""Offline Complex Portal GUI regression tests."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from pichanalysis.core.complex_analysis import ComplexParameters, MissingComplexOutputError, read_complex_outputs, run_complex_analysis
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.ui.complexes_page import COVERAGE_LABELS, ComplexPage
from smoke_complex_analysis import build_fixture


def test_complex_gui_python_r_history(tmp_path):
    app = QApplication.instance() or QApplication([])
    project, manager = build_fixture(tmp_path)
    page = ComplexPage(manager)
    page.set_project(project)
    assert "Ready" in page.database.text()
    assert page.run_button.isEnabled()
    assert page.fdr.value() == .05 and page.minimum_overlap.value() == 3 and page.top_n.value() == 20
    assert "P12345-2" in page.manual.item(1,1).text()
    assert page.manual.item(1,2).text()=="P12345"
    assert "Complex Portal lookup: P12345" in page.manual.item(1,1).toolTip()
    page.target.setCurrentIndex(page.target.findData("Manual selection"))
    page.manual.selectRow(0)
    page.manual.selectRow(2)
    page.background.setCurrentIndex(page.background.findData("Manual selection"))
    page.background_manual.item(0).setSelected(True)
    selected=page.parameters()
    assert selected["manual_rows"]==(1,3)
    assert selected["background_manual_rows"]==(1,)
    assert selected["target_selection"]==selected["background_selection"]=="Manual selection"
    page.target.setCurrentIndex(0)
    page.background.setCurrentIndex(0)
    runtime = RRuntime()
    if not runtime.available:
        import pytest
        pytest.skip("Rscript unavailable")
    a = run_complex_analysis(project,manager,runtime,run_id="gui_a",
        parameters=ComplexParameters(target_selection="Manual selection",manual_rows=(1,2,3,4,5,6)))
    page.show_outputs(a)
    page.show_complex_details("CPX-1")
    assert page.coverage_progress.format() == "Protein-component coverage: 3 / 3 groups"
    assert page.detail_tables["component_groups"].rowCount() == 3
    assert page.detail_tables["nonprotein_participants"].rowCount() == 2
    assert page.detail_tables["nested_complexes"].rowCount() == 1
    groups = page._frames[page.detail_tables["component_groups"]]
    alt = groups[groups.component_type == "alternative"].iloc[0]
    assert alt.detected_uniprot_accessions == "Q11111;Q22222"
    assert a["tables"]["complex_frequency"].set_index("complex_id").loc["CPX-1","target_member_count"] == 4
    assert page.detail_tables["direct_participants"].rowCount() > 0
    assert page.detail_tables["expanded_protein_components"].rowCount() > 0
    memberships=page._frames[page.tables["protein_to_complexes"]]
    protein=memberships.canonical_uniprot.iloc[0]
    row=int(memberships.index[memberships.canonical_uniprot==protein][0])
    page.tables["protein_to_complexes"].selectRow(row)
    assert page.protein_memberships.rowCount()==int((memberships.canonical_uniprot==protein).sum())
    assert page.graph_choice.count() == 6
    for code,label in COVERAGE_LABELS.items():
        assert label in [page.coverage_filter.itemText(i) for i in range(page.coverage_filter.count())]
    b_snapshot = manager.complex_portal.install_from_file(tmp_path/"9606.tsv")
    b = run_complex_analysis(project,manager,runtime,run_id="gui_b",
        parameters=ComplexParameters(target_selection="Manual selection",manual_rows=(1,3)))
    page.show_outputs(b)
    assert page.outputs["metadata"]["snapshot_id"] == b_snapshot.name
    page._load_history(page.history.findData("gui_a"))
    assert page.outputs["metadata"]["snapshot_id"] == a["metadata"]["snapshot_id"]
    assert page.tables["summary"].rowCount() >= 20
    assert page.tables["complex_coverage"].rowCount() == len(a["tables"]["complex_coverage"])
    assert page.tables["component_groups"].rowCount() == len(a["tables"]["component_groups"])
    assert page.tables["complex_frequency"].rowCount() > 0
    assert page.tables["complex_enrichment_all"].rowCount() > 0
    assert page.tables["protein_to_complexes"].rowCount() > 0
    assert page.tables["complex_to_proteins"].rowCount() > 0
    memberships=page._frames[page.tables["protein_to_complexes"]]
    protein=memberships.canonical_uniprot.iloc[0]
    row=int(memberships.index[memberships.canonical_uniprot==protein][0])
    page.tables["protein_to_complexes"].selectRow(row)
    assert page.protein_memberships.rowCount()==int((memberships.canonical_uniprot==protein).sum())
    assert page.graph_choice.count() == 6
    assert read_complex_outputs(project,"gui_a")["metadata"]["snapshot_id"] != b["metadata"]["snapshot_id"]
    (a["run_root"]/"plots/coverage_classes.pdf").unlink()
    a["workbook"].unlink()
    optional = read_complex_outputs(project,"gui_a")
    assert optional["workbook"] is None and len(optional["plots"]) == 11
    page.show_outputs(optional)
    (a["run_root"]/"summary.csv").unlink()
    import pytest
    with pytest.raises(MissingComplexOutputError):
        read_complex_outputs(project,"gui_a")
    page._load_history(page.history.findData("gui_a"))
    assert page.outputs is None and "incomplete" in page.history_status.text()
    page.close()
    app.processEvents()
