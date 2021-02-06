from .mmseqs_native import Databases

VERSION = "1.0.0"


class MMSeqs:
    def __init__(self, storage_directory: str = "mmseqs_storage"):
        self.databases = Databases(storage_directory, VERSION)


__all__ = [
    "MMSeqs",
]
