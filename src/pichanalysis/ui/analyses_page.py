from PySide6.QtWidgets import QLabel, QListWidget, QListWidgetItem, QVBoxLayout, QWidget


class AnalysesPage(QWidget):
    def __init__(self) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        message = QLabel("Nenhuma análise científica implementada nesta versão.")
        message.setWordWrap(True)
        layout.addWidget(message)
        modules = QListWidget()
        for name in (
            "Identificação e anotação", "Presença/ausência", "Diferencial", "GO",
            "KEGG", "Reactome", "MitoCarta", "Domínios", "STRING / redes", "Complexos",
        ):
            item = QListWidgetItem(name)
            item.setFlags(item.flags() & ~item.flags().ItemIsEnabled)
            modules.addItem(item)
        layout.addWidget(modules)

