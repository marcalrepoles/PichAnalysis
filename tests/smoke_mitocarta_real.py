"""Manual official MitoCarta3.0 smoke; intentionally excluded from pytest collection."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.databases.mitocarta import FILES, MitoCartaProvider


def main() -> int:
    temporary=Path(tempfile.mkdtemp(prefix="pichanalysis-mitocarta-real-"));report={"temporary_path":str(temporary)}
    try:
        manager=DatabaseManager(temporary/"databases");snapshot=manager.mitocarta.create_staging_snapshot();manager.mitocarta.mark_downloading(snapshot);provider=MitoCartaProvider(retries=2,timeout=120);files=[]
        for name in FILES:files.append(provider.download(name,snapshot/"raw"/name))
        snapshot=manager.mitocarta.finalize_snapshot(snapshot,files);manifest=manager.mitocarta.manifest();genes=pd.read_csv(manager.mitocarta.genes_path());sub=pd.read_csv(manager.mitocarta.subcompartment_membership_path());members=pd.read_csv(manager.mitocarta.pathway_membership_path());hierarchy=pd.read_csv(manager.mitocarta.table_path("hierarchy"))
        report.update({"snapshot_path":str(snapshot),"snapshot_id":snapshot.name,"source_files":files,"workbook_sheets":manifest["workbook_sheets"],"workbook_sheet_used":manifest["workbook_sheet_used"],"source_columns":manifest["source_columns"],"total_workbook_rows":manifest["total_workbook_rows"],"mitocarta_genes_detected":int(genes.is_mitocarta.sum()),"mitopathways_detected":manifest["mitopathway_count"],"subcompartment_memberships":len(sub),"pathway_memberships":len(members),"hierarchy_rows":len(hierarchy),"expected_genes":1136,"expected_pathways":149,"ready":manager.mitocarta.is_ready()})
        assert report["mitocarta_genes_detected"]==1136 and report["mitopathways_detected"]==149 and report["ready"]
        print("MITOCARTA_REAL_SMOKE "+json.dumps(report,indent=2));return 0
    except Exception as error:
        report["error"]=repr(error);print("MITOCARTA_REAL_SMOKE_FAILED "+json.dumps(report,indent=2),file=sys.stderr);return 1


if __name__=="__main__":raise SystemExit(main())
