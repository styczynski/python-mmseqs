from .mmseqs_native import Client, SearchResults, Database
import pandas as pd

VERSION = "1.0.0"

def client_get_databases_dataframe(self):
    data = dict(name=[], description=[], type=[])
    for db in self:
        for key in data.keys():
            data[key].append(getattr(db, key))
    return pd.DataFrame(data)

Client.dataframe = property(client_get_databases_dataframe)

def search_results_get_dataframe(self):
    data = dict()
    for header in self._headers:
        header_data = []
        for row in self._records:
            for aln in row:
                header_data.append(getattr(aln, header))
        data[header] = header_data
    return pd.DataFrame(data)

SearchResults.dataframe = property(search_results_get_dataframe)

def database_get_dataframe(self):
    data = self.columns_data
    return pd.DataFrame(dict(
        id=data[0],
        length=data[2],
        sequence=data[1],
    ))

Database.dataframe = property(database_get_dataframe)

class MMSeqs:
    def __init__(self, storage_directory: str = "mmseqs_storage"):
        self.databases = Client(storage_directory, VERSION)


__all__ = [
    "MMSeqs",
]
