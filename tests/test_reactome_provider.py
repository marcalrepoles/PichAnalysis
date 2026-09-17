import io,tarfile
import pytest
from pichanalysis.core.database_registry import DATABASE_REGISTRY
from pichanalysis.core.databases.reactome import safe_extract_tar
from pichanalysis.core.databases.reactome import CORE_FILES
from pichanalysis.core.reactome_database import ReactomeDatabase
from pichanalysis.core.database_registry import DatabaseState
def test_reactome_registry():
 spec=DATABASE_REGISTRY["reactome"];assert spec.taxonomy_id=="9606" and spec.organism_code=="R-HSA-"
def test_safe_tar_rejects_traversal(tmp_path):
 archive=tmp_path/"bad.tgz"
 with tarfile.open(archive,"w:gz") as tar:
  data=b"bad";info=tarfile.TarInfo("../escape.txt");info.size=len(data);tar.addfile(info,io.BytesIO(data))
 with pytest.raises(ValueError,match="Unsafe"):safe_extract_tar(archive,tmp_path/"out")
 assert not (tmp_path/"escape.txt").exists()

def core_fixture(root,value="Pathway"):
 root.mkdir();data={"ReactomePathways.txt":f"R-HSA-1\t{value}\tHomo sapiens\n","ReactomePathwaysRelation.txt":"R-HSA-1\tR-HSA-2\n","UniProt2Reactome.txt":"P04637\tR-HSA-1\tTP53\tEvidence\tHomo sapiens\n","NCBI2Reactome.txt":"7157\tR-HSA-1\tTP53\tEvidence\tHomo sapiens\n","humanPathwaysWithDiagrams.txt":"R-HSA-1\n","pathway2summation.txt":"R-HSA-1\tA valid summary\n"}
 for name,text in data.items():(root/name).write_text(text,encoding="utf-8")
 return root
def test_valid_snapshot_manifest_hashes_and_active(tmp_path):
 source=core_fixture(tmp_path/"source");db=ReactomeDatabase(tmp_path/"db");snap=db.install_from_directory(source);m=db.manifest()
 assert m["status"]==DatabaseState.READY and db.active_snapshot()==snap and len(m["files"])==6 and db.is_core_ready()
 import hashlib
 for name,path in db.core_paths().items():assert db.source_hashes()[name]==hashlib.sha256(path.read_bytes()).hexdigest()
def test_missing_file_is_not_ready(tmp_path):
 source=core_fixture(tmp_path/"source");(source/CORE_FILES[0]).unlink();db=ReactomeDatabase(tmp_path/"db")
 with pytest.raises(FileNotFoundError):db.install_from_directory(source)
 assert not db.is_core_ready()
def test_structurally_invalid_file_is_not_ready(tmp_path):
 source=core_fixture(tmp_path/"source");(source/"ReactomePathways.txt").write_text("invalid",encoding="utf-8");db=ReactomeDatabase(tmp_path/"db")
 with pytest.raises(ValueError,match="ReactomePathways"):db.install_from_directory(source)
 assert not db.is_core_ready()
def test_failed_update_preserves_active_snapshot(tmp_path):
 db=ReactomeDatabase(tmp_path/"db");a=db.install_from_directory(core_fixture(tmp_path/"a","A"));before=(a/"raw/ReactomePathways.txt").read_bytes();bad=core_fixture(tmp_path/"bad");(bad/"ReactomePathways.txt").write_text("bad",encoding="utf-8")
 with pytest.raises(ValueError):db.install_from_directory(bad)
 assert db.active_snapshot()==a and (a/"raw/ReactomePathways.txt").read_bytes()==before
def test_valid_update_activates_b_and_preserves_a(tmp_path):
 db=ReactomeDatabase(tmp_path/"db");a=db.install_from_directory(core_fixture(tmp_path/"a","A"));b=db.install_from_directory(core_fixture(tmp_path/"b","B"));assert db.active_snapshot()==b and a.is_dir() and b!=a
def test_staging_is_incomplete_and_does_not_activate(tmp_path):
 db=ReactomeDatabase(tmp_path/"db");snap=db.create_staging_snapshot();assert db.manifest(snap)["status"]==DatabaseState.INCOMPLETE and db.active_snapshot() is None
