"""Multi-module persisted-identity bridges, without scientific recomputation."""
import pandas as pd
from PySide6.QtWidgets import QApplication

from pichanalysis.core.cross_module_integration import (
    CrossModuleAdapter, CrossModuleAdapterRegistry, default_registry,
    rebuild_cross_module_index,
)
from pichanalysis.core.project import create_project
from pichanalysis.ui.cross_module_explorer import CrossModuleExplorer


MODULES = ("mapping", "presence_absence", "go", "kegg", "reactome",
    "mitocarta", "domains", "string", "complexes", "mtdna_evidence",
    "proteomics_qc", "differential")


def test_default_registry_has_all_explicit_adapters():
    assert {adapter.module_id for adapter in default_registry().values()} == set(MODULES)


def test_multi_module_offline_bridges_from_persisted_ids(tmp_path):
    project = create_project(tmp_path, "integration")
    frames = {
        "source": pd.DataFrame([{"feature_id": "ROW000001", "source_row": "1",
            "original_identifier": "P12345-2", "uniprot_accession": "P12345-2",
            "ncbi_gene_id": "100"}]),
        "go": pd.DataFrame([{"input_id": "P12345-2", "GO_ID": "GO:0001",
            "uniprot_accession": "P12345-2", "ncbi_gene_id": "100"}]),
        "kegg": pd.DataFrame([{"kegg_gene_id": "hsa:100", "ncbi_gene_id": "100"}]),
        "reactome": pd.DataFrame([{"reactome_entity_key": "R-HSA-1",
            "ncbi_gene_id": "100"}]),
        "mitocarta": pd.DataFrame([{"entity_key": "NCBI:100", "ncbi_gene_id": "100"}]),
        "domains": pd.DataFrame([{"uniprot_accession": "P12345-2",
            "interpro_id": "IPR000001"}]),
        "string": pd.DataFrame([{"uniprot_accession": "P12345-2",
            "string_protein_id": "9606.ENSP1"}]),
        "complexes": pd.DataFrame([{"canonical_uniprot": "P12345-2",
            "complex_id": "CPX-1"}]),
        "mtdna_evidence": pd.DataFrame([{"entity_key": "NCBI:100",
            "ncbi_gene_id": "100"}]),
    }
    registry = CrossModuleAdapterRegistry()
    for module_id, frame in frames.items():
        registry.register(CrossModuleAdapter(module_id,
            lambda _project, _module=module_id: [f"{_module}_A"],
            lambda _project, run_id, _module=module_id, _frame=frame: {
                "run_root": project.root / "analyses" / _module / "runs" / run_id,
                "metadata": {"run_id": run_id, "snapshot_id": f"snapshot_{_module}"},
                "tables": {"records": _frame}}, ("records",)))
    index = rebuild_cross_module_index(project, registry)
    source = next(row for row in index.search("ROW000001") if row["module"] == "source")
    context = index.context(source["id"])
    assert context.uniprot_accessions == ("P12345-2",)
    assert context.ncbi_gene_ids == ("100",)
    assert not index.matches(context)
    matches = index.matches(context, include_other_lineages=True)
    assert {match.target_module for match in matches} == set(frames) - {"source"}
    assert all(not match.lineage_compatible for match in matches)
    assert all("different or unverified input lineage" in match.basis for match in matches)
    assert not any(match.source_identifier == "P12345" for match in matches)
    app = QApplication.instance() or QApplication([])
    explorer = CrossModuleExplorer()
    explorer.set_project(project)
    explorer.index = index
    explorer._ensure_index = lambda: True
    explorer.search.setText("ROW000001")
    explorer._search()
    assert len(explorer.source_records) == 1
    explorer.source_table.selectRow(0)
    explorer.cross_lineage.setChecked(True)
    assert {match.target_module for match in explorer.matches} == set(frames) - {"source"}
    opened = []
    explorer.open_target_requested.connect(lambda *args: opened.append(args))
    assert not explorer.open_button.isEnabled()
    target_row = next(i for i, match in enumerate(explorer.matches)
        if match.target_module == "string")
    explorer.matches_table.selectRow(target_row)
    assert explorer.open_button.isEnabled()
    explorer.open_button.click()
    assert opened[0][:2] == ("string", "string_A")
    explorer.close()
