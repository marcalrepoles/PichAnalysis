from pathlib import Path
def test_key_ui_labels_are_english():
    root=Path(__file__).parents[1]/"src/pichanalysis/ui";text="\n".join(path.read_text(encoding="utf-8") for path in root.glob("*.py"))
    for label in ("Database Manager","Open existing project","Import XLSX, CSV, or TSV","Presence / absence","KEGG academic-use notice"):assert label in text
    for label in ("Abrir projeto existente","Importar XLSX","Nenhum resultado disponível","Configuração das colunas"):assert label not in text
