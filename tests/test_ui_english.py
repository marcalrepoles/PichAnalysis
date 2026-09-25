from pathlib import Path


def test_reviewed_runtime_text_is_english():
    root = Path(__file__).parents[1] / "src" / "pichanalysis"
    ui = "\n".join(path.read_text(encoding="utf-8") for path in (root / "ui").glob("*.py"))
    core = "\n".join(path.read_text(encoding="utf-8") for path in (root / "core").glob("*.py"))
    for label in (
        "Database Manager", "Open existing project", "Import XLSX, CSV, or TSV",
        "Presence / absence", "KEGG academic-use notice", "Identification and Annotation",
        "Application R scripts", "All evidence", "Rows analyzed:",
    ):
        assert label in ui
    for label in (
        "Novo projeto", "Nome científico", "Identificação e anotação",
        "Aplicar sugestões automáticas", "Pré-visualização", "Tabela e colunas",
        "Linhas separadas por vírgula", "Todas as evidências", "Executar GO",
        "Número mínimo de réplicas", "Zero será tratado como ausência",
        "Testar R", "Scripts R do aplicativo", "Importação concluída",
    ):
        assert label not in ui
    for message in (
        "A coluna", "Coluna sem valores", "Não foi possível", "Abra um projeto",
        "Pronto para análise GO", "Todas as proteínas", "Pronto para executar",
        "Nenhum identificador principal", "Resultados GO inválidos",
    ):
        assert message not in core


def test_reviewed_r_messages_and_plot_labels_are_english():
    root = Path(__file__).parents[1] / "r_scripts"
    text = "\n".join(path.read_text(encoding="utf-8") for path in root.rglob("*.R") if "tests" not in path.parts)
    for label in ("Condition / replicate", "GO frequency", "Missing R packages:"):
        assert label in text
    for label in (
        "Pacotes R ausentes", "Argumentos devem", "Coluna identificadora ausente",
        "Proteínas detectadas", "Fração detectada", "Termo GO", "Detectada",
    ):
        assert label not in text
