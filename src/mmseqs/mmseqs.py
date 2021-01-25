import os
import pathlib
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

from .mmseqs_native import MMSeqsCallConfig, _call_mmseqs
from .databases import Databases
from .base import MetaDatabase, MMSeqsSettings


class MMSeqs:
    def __init__(self, storage_directory='mmseqs_storage'):
        self.settings = MMSeqsSettings(
            storage_directory=storage_directory,
            seq_storage_directory=os.path.join(storage_directory, 'databases'),
            tmp_directory=os.path.join(storage_directory, 'workdir', ''),
            seq_results_directory=os.path.join(storage_directory, 'results', ''),
            meta_db=MetaDatabase(os.path.join(storage_directory, 'mmseqs_meta.sqlite')),
        )
        self.databases = Databases(self.execute, self.settings)
        for directory in [
            self.settings.storage_directory,
            self.settings.seq_storage_directory,
            self.settings.tmp_directory,
            self.settings.seq_results_directory]:
            if not os.path.isdir(directory):
                pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
        with self.settings.meta_db.open() as meta_db:
            if not meta_db.has('databases'):
                meta_db.set('databases', [])

        pass

    def execute(self, command_name, args):
        cli_args = MMSeqsCallConfig()
        for key in args.keys():
            setattr(cli_args, key, args[key])
            if key in ["db1", "db2", "db3"]:
                setattr(cli_args, f"{key}Index", f'{args[key]}.index')
                setattr(cli_args, f"{key}dbtype", f'{args[key]}.dbtype')
                hdr_key = f'hdr{key.replace("db", "")}'
                setattr(cli_args, f"{hdr_key}", f'{args[key]}_h')
                setattr(cli_args, f"{hdr_key}Index", f'{args[key]}_h.index')
                setattr(cli_args, f"{hdr_key}dbtype", f'{args[key]}_h.dbtype')
        print(f'execute ==> {args}')
        cli_args.baseTmpPath = self.settings.tmp_directory
        return _call_mmseqs(command_name, cli_args)
