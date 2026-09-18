from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

from pichanalysis.core.interpro_database import InterProDatabase


def main():
    temporary=Path(tempfile.mkdtemp(prefix="pichanalysis_interpro_smoke_"));database=InterProDatabase(temporary/"global")
    metadata=database.download_metadata();accessions=["P04637","P38398"]
    root=database.build_annotation_set(temporary/"project",accessions,annotation_set_id="official_smoke")
    loaded=database.load_annotation_set(temporary/"project","official_smoke");hash_index=root/"raw_response_hashes.csv"
    database.provider=object()
    offline={accession:{"interpro":len(database.find_interpro_memberships(temporary/"project",accession,"official_smoke")),"pfam":len(database.find_pfam_memberships(temporary/"project",accession,"official_smoke")),"interpro_locations":len(database.find_interpro_locations(temporary/"project",accession,"official_smoke")),"pfam_locations":len(database.find_pfam_locations(temporary/"project",accession,"official_smoke"))} for accession in accessions}
    assert all(values["interpro"] and values["pfam"] for values in offline.values())
    report={"temporary_path":str(temporary),"annotation_set_id":loaded.manifest["annotation_set_id"],"requested_uniprot_accessions":accessions,"interpro_release":loaded.manifest["interpro_release"],"pfam_release":loaded.manifest["pfam_release"],"raw_response_files":len(list((root/"raw").glob("*/*.json"))),"raw_response_hash_index_sha256":hashlib.sha256(hash_index.read_bytes()).hexdigest(),"interpro_memberships":len(loaded.interpro_membership),"pfam_memberships":len(loaded.pfam_membership),"interpro_locations":len(loaded.interpro_locations),"pfam_locations":len(loaded.pfam_locations),"no_match_count":loaded.manifest["no_match_protein_count"],"failed_count":loaded.manifest["failed_protein_count"],"metadata_snapshot":metadata.name,"offline_lookup":offline}
    print("INTERPRO_REAL_SMOKE "+json.dumps(report,sort_keys=True))

if __name__=="__main__":main()
