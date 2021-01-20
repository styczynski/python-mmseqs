from typing import Sequence, List, Optional
import os
from dataclasses import dataclass
from glob import glob
from datetime import date

from .utils import remove_paths, to_args
from .base import MMSeqsBase

from .mmseqs_native import MMSeqsMultiParamString

PARAM_DB_TYPE_MAPPING = dict(
    auto=0,
    amino_acid=1,
    nucleotides=2,
)

PARAM_CREATEDB_MODE_MAPPING = dict(
    copy=0,
    soft_link=1,
)

PARAM_DB_SEARCH_TYPE_MAPPING = dict(
    auto=0,
    protein=1,
    translated=2,
    nucleotides=3,
    translated_nucleotides_aligned=4,
)

@dataclass
class IndexStats:
    db_size: int
    index_entries: int
    avg_kmer_size: float

@dataclass
class Database:
    """Class for keeping track of an item in inventory."""
    name: str
    description: str
    input_files: List[str]
    database_type: str
    created_on: date

    def copy(self, name):
        with self._base.settings.meta_db.open() as meta_db:
            src_dbs, _ = meta_db.list_filter('databases', self._base, lambda db: db.name != self.name)
            src_db = src_dbs[0]
            db_path = os.path.join(self._base.settings.seq_storage_directory, src_db.name)
            db_copy_path = os.path.join(self._base.settings.seq_storage_directory, name)
            self._base._execute_cli("cpdb", dict(db1=db_path, db2=db_copy_path))
            new_db = Database(
                name=name,
                description=src_db.description,
                input_files=src_db.input_files,
                created_on=date.today(),
                database_type=src_db.database_type,
            )
            meta_db.list_append('databases', new_db)
        return new_db

    def remove(self):
        with self._base.settings.meta_db.open() as meta_db:
            removed_dbs, left_dbs = meta_db.list_filter('databases', self._base, lambda db: db.name != self.name)
            for db in removed_dbs:
                db_path = os.path.join(self._base.settings.seq_storage_directory, db.name)
                self._base._execute_cli("rmdb", dict(db1=db_path))
                self._base._execute_cli("rmdb", dict(db1=f'{db_path}_h'))
            meta_db.set('databases', left_dbs)

    def search(self, search_input: str):
        tmp_dir = 'tmp'
        results_path = os.path.join(self._base.settings.seq_results_directory, f'{self.name}.query_results.m8')
        seq_db_path = os.path.join(self._base.settings.seq_storage_directory, self.name)
        out = self._base._execute_cli('easy-search', dict(
            filenames=[search_input, seq_db_path, results_path, tmp_dir],
            shuffleDatabase=False,
            sensitivity=5.7,
            removeTmpFiles=True,
            writeLookup=False,
            alignmentMode=3,
        ))
        print(out.vars_str)

    def create_index(self, search_type: str = 'nucleotides') -> IndexStats:
        tmp_dir = 'tmp'
        seq_db_path = os.path.join(self._base.settings.seq_storage_directory, self.name)

        # par.orfStartMode = 1;
        # //    par.orfMinLength = 30;
        # //    par.orfMaxLength = 32734;
        # //    par.kmerScore = 0; // extract all k-mers
        # //    par.maskMode = 0;
        # //    par.spacedKmer = false;
        # //    // VTML has a slightly lower sensitivity in the regression test
        # //    par.seedScoringMatrixFile = MultiParam<char*>("blosum62.out", "nucleotide.out");
        # //
        seed_scoring_matrix_file = MMSeqsMultiParamString()
        seed_scoring_matrix_file.aminoacids = "blosum62.out"
        seed_scoring_matrix_file.nucleotides = "nucleotide.out"
        out = self._base._execute_cli('createindex', dict(
            db1=seq_db_path,
            db2=tmp_dir,
            filenames=[tmp_dir],
            searchType=PARAM_DB_SEARCH_TYPE_MAPPING[search_type],
            orfStartMode=1,
            orfMinLength=30,
            orfMaxLength=32734,
            kmerScore=0,
            maskMode=1,
            sensitivity=7.5,
            removeTmpFiles=True,
        ))
        print("WELP DONE!")
        return IndexStats(
            db_size=int(out.vars_str['INDEX_TABLE_DB_SIZE']),
            index_entries=int(out.vars_str['INDEX_TABLE_ENTRIES']),
            avg_kmer_size=float(out.vars_str['INDEX_AVG_KMER_SIZE']),
        )

def inject_base(objs: List[any]):
    for obj in objs:
        setattr(obj, '_base')


class Databases(MMSeqsBase):
    def __init__(self, *args):
        super().__init__(*args)

    def __getitem__(self, index):
        return self.list()[index]

    def list(self) -> List[Database]:
        with self.settings.meta_db.open() as meta_db:
            return meta_db.list_get('databases', self)

    def create(self,
               name: str,
               description: str,
               input_files: Sequence[str],
               mode: str = "copy",
               database_type: str = "auto",
               offset: int = 0,
               shuffle: bool = False,
    ) -> Database:
        """
        Create a sequence database from multiple FASTA file

        :param name: Name of the new database
        :param description: Description of the new database
        :param input_files: FASTA/Q input files
        :param mode: Database creation mode
            copy - Default mode
            soft_link - Soft link data and write new index (works only with single line fasta/q)
        :param database_type: Database type
            auto - Automatically infer the type
            amino_acid - Database with amino acids sequences
            nucleotides - Database with nucleotide sequences
        :param offset: Numeric ids in index file are offset by this value
        :param shuffle: Shuffle the input database
        :return:
        """
        input_files = list(input_files)
        self._execute_cli("createdb", dict(
            filenames=[*input_files, os.path.join(self.settings.seq_storage_directory, name)],
            identifierOffset=offset,
            dbType=PARAM_DB_TYPE_MAPPING[database_type],
            createdbMode=PARAM_CREATEDB_MODE_MAPPING[mode],
            shuffleDatabase=int(shuffle),
        ))
        with self.settings.meta_db.open() as meta_db:
            new_db = Database(
                name=name,
                description=description,
                input_files=input_files,
                created_on=date.today(),
                database_type=database_type,
            )
            meta_db.list_append('databases', new_db)
        return new_db

