import os
import pathlib
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

from .mmseqs_native import MMSeqsCallArgs, _call_mmseqs
from .databases import Databases
from .base import MetaDatabase, MMSeqsSettings


class MMSeqs:
    def __init__(self, storage_directory='mmseqs_storage'):
        self.settings = MMSeqsSettings(
            storage_directory=storage_directory,
            seq_storage_directory=os.path.join(storage_directory, 'databases'),
            meta_db=MetaDatabase(os.path.join(storage_directory, 'mmseqs_meta.sqlite')),
        )
        self.databases = Databases(self.execute, self.settings)
        for directory in [self.settings.storage_directory, self.settings.seq_storage_directory]:
            if not os.path.isdir(directory):
                pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
        with self.settings.meta_db.open() as meta_db:
            if not meta_db.has('databases'):
                meta_db.set('databases', [])

        pass


    def execute(self, cli_args):
        print(f'execute ==> {repr(cli_args)}')
        call_args = MMSeqsCallArgs()
        call_args.cli_args = ["mmseq2", *[str(arg) for arg in cli_args]]
        return _call_mmseqs(call_args)