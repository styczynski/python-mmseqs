from .mmseqs_native import Databases, SearchResults
import pandas as pd

VERSION = "1.0.0"

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


class MMSeqs:
    def __init__(self, storage_directory: str = "mmseqs_storage"):
        self.databases = Databases(storage_directory, VERSION)


__all__ = [
    "MMSeqs",
]
