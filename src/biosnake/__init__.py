from .biosnake_native import Client, SearchResults, Database
from . import handlers

VERSION = "1.0.0"

Client.dataframe = property(handlers._client_get_databases_dataframe)
SearchResults.dataframe = property(handlers._search_results_get_dataframe)
SearchResults.to_fasta = handlers._search_results_to_fasta
Database.dataframe = property(handlers._database_get_dataframe)

class Biosnake:
    databases: Client

    def __init__(self, storage_directory: str = "biosnake_storage"):
        self.databases = Client(storage_directory, VERSION)


__all__ = [
    "Biosnake",
]
