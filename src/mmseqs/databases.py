from typing import Sequence, List, Optional
import os
from dataclasses import dataclass
from glob import glob
from datetime import date
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

from .utils import remove_paths, to_args
from .base import MMSeqsBase, ObjectWithBaseRef
from .convert import GenericSequencesFilePathConverter, GenericSequences

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
class SearchRecord:
    query_sequence_id: str
    target_sequence_id: str
    sequence_identity: float
    alignment_length: int
    number_of_mismatches: int
    number_of_gap_openings: int
    domain_start_index_query: int
    domain_end_index_query: int
    domain_start_index_target: int
    domain_end_index_target: int
    e_value: float
    bit_score: int

@dataclass
class IndexStats:
    db_size: int
    index_entries: int
    avg_kmer_size: float

@dataclass
class SearchResults:
    records: Sequence[Sequence[SearchRecord]]

    @property
    def dataframe(self):
        return pd.DataFrame([[record.query_sequence_id,
                            record.target_sequence_id,
                            record.sequence_identity,
                            record.alignment_length,
                            record.number_of_mismatches,
                            record.number_of_gap_openings,
                            record.domain_start_index_query,
                            record.domain_end_index_query,
                            record.domain_start_index_target,
                            record.domain_end_index_target,
                            record.e_value,
                            record.bit_score] for alignments in self.records for record in alignments], columns=[
            'query_sequence_id',
            'target_sequence_id',
            'sequence_identity',
            'alignment_length',
            'number_of_mismatches',
            'number_of_gap_openings',
            'domain_start_index_query',
            'domain_end_index_query',
            'domain_start_index_target',
            'domain_end_index_target',
            'e_value',
            'bit_score',
        ])

class ComputedDatabase(ObjectWithBaseRef):
    _operation: str
    _args: List[any]
    _computed_db: any

    def __init__(self, operation, args):
        super().__init__()
        self._operation = operation
        self._args = args
        self._computed_db = None

    def __add__(self, other):
        db = ComputedDatabase(operation='merge', args=[self, other])
        setattr(db, '_base', self._base)
        return db

    def __sub__(self, other):
        db = ComputedDatabase(operation='substract', args=[self, other])
        setattr(db, '_base', self._base)
        return db

    def __truediv__(self, chunks_no):
        db = ComputedDatabase(operation='split', args=[self, chunks_no])
        setattr(db, '_base', self._base)
        return db

    @property
    def _name(self):
        if self._operation == 'merge':
            return f'{self._args[0].name}_merge_{self._args[1].name}'
        elif self._operation == 'substract':
            return f'{self._args[0].name}_substract_{self._args[1].name}'
        elif self._operation == 'split':
            return f'{self._args[0].name}_split_{self._args[1]}'

    def _compute(self):
        if self._computed_db is not None:
            return self._computed_db
        new_db = None
        if self._operation == 'merge':
            # merge dbs
            with self._base.settings.meta_db.open() as meta_db:
                path_in_a = os.path.join(self._base.settings.seq_storage_directory, self._args[0].name)
                path_out = os.path.join(self._base.settings.seq_storage_directory, self._name)
                path_in_b = os.path.join(self._base.settings.seq_storage_directory, self._args[1].name)
                self._base._execute_cli("mergedbs", dict(
                    db1=path_in_a,
                    db2=path_out,
                    db3=path_in_b,
                    filenames=[path_in_a, path_out, path_in_b]))
                new_db = Database(
                    name=self._name,
                    description="Merged database",
                    input_files=[*self._args[0].input_files, *self._args[1].input_files],
                    created_on=date.today(),
                    database_type=self._args[0].database_type,
                )
                meta_db.list_append('databases', new_db)
        elif self._operation == 'substract':
            # substract dbs
            with self._base.settings.meta_db.open() as meta_db:
                path_in_a = os.path.join(self._base.settings.seq_storage_directory, self._args[0].name)
                path_out = os.path.join(self._base.settings.seq_storage_directory, self._name)
                path_in_b = os.path.join(self._base.settings.seq_storage_directory, self._args[1].name)
                self._base._execute_cli("subtractdbs", dict(
                    db1=path_in_a,
                    db2=path_in_b,
                    db3=path_out,
                    filenames=[path_in_a, path_in_b, path_out]))
                new_db = Database(
                    name=self._name,
                    description="Substracted databases",
                    input_files=[*self._args[0].input_files, *self._args[1].input_files],
                    created_on=date.today(),
                    database_type=self._args[0].database_type,
                )
                meta_db.list_append('databases', new_db)
        elif self._operation == 'split':
            # substract dbs
            with self._base.settings.meta_db.open() as meta_db:
                path_in_a = os.path.join(self._base.settings.seq_storage_directory, self._args[0].name)
                path_out = os.path.join(self._base.settings.seq_storage_directory, self._name)
                self._base._execute_cli("splitdb", dict(
                    db1=path_in_a,
                    db2=path_out,
                    filenames=[path_in_a, path_out],
                    split=int(self._args[1])))
                new_db = Database(
                    name=self._name,
                    description="Splitted database",
                    input_files=[self._args[0].input_files],
                    created_on=date.today(),
                    database_type=self._args[0].database_type,
                    chunks=[],
                )
                meta_db.list_append('databases', new_db)
        self._computed_db = new_db
        setattr(self._computed_db, '_base', self._base)
        return self._computed_db

    def __getattr__(self, name):
        return getattr(self._compute(), name)

class Database(ObjectWithBaseRef):
    """Class for keeping track of an item in inventory."""
    description: str
    input_files: List[str]
    database_type: str
    created_on: date
    _name: str

    def __init__(self, description="", input_files=None, database_type="", created_on=None, name=""):
        self.description = description
        self.input_files = input_files
        self.database_type = database_type
        self.created_on = created_on
        self._name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        original_name = self._name
        with self._base.settings.meta_db.open() as meta_db:
            path_in = os.path.join(self._base.settings.seq_storage_directory, original_name)
            path_out = os.path.join(self._base.settings.seq_storage_directory, name)
            self._base._execute_cli("mvdb", dict(
                db1=path_in,
                db2=path_out,
                filenames=[path_in, path_out]))
            self._name = name
            meta_db.list_replace('databases', self._base, self, lambda x: x.name == original_name, save_back=True)

    @property
    def records(self):
        db_path = os.path.join(self._base.settings.seq_storage_directory, self._name)
        with open(f'{db_path}', 'r') as db_file, open(f'{db_path}_h', 'r') as h_file:
            for seq_contents, header in zip(db_file, h_file):
                if len(seq_contents) > 1:
                    yield SeqRecord(Seq(seq_contents[:-1].lstrip('\x00')), id=header[:-1].lstrip('\x00'))

    def to_fasta(self, output_path):
        # convert2fasta
        db_path = os.path.join(self._base.settings.seq_storage_directory, self.name)
        self._base._execute_cli('convert2fasta', dict(
            db1=db_path,
            db2=output_path,
            filenames=[db_path, output_path],
        ))
        #with open(output_path, 'w') as output_handle:
        #    SeqIO.write(self.records, output_handle, "fasta")

    def __add__(self, other):
        db = ComputedDatabase(operation='merge', args=[self, other])
        setattr(db, '_base', self._base)
        return db

    def __sub__(self, other):
        db = ComputedDatabase(operation='substract', args=[self, other])
        setattr(db, '_base', self._base)
        return db

    def __truediv__(self, chunks_no):
        db = ComputedDatabase(operation='split', args=[self, chunks_no])
        setattr(db, '_base', self._base)
        return db

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
            removed_dbs, left_dbs = meta_db.list_filter('databases', self._base, lambda db: db.name != self._name)
            for db in removed_dbs:
                db_path = os.path.join(self._base.settings.seq_storage_directory, db.name)
                self._base._execute_cli("rmdb", dict(db1=db_path))
                self._base._execute_cli("rmdb", dict(db1=f'{db_path}_h'))
            meta_db.set('databases', left_dbs)

    def search(self,
       search_input: GenericSequences,
       search_type: str = 'auto',
    ) -> SearchResults:
        with GenericSequencesFilePathConverter(search_input) as input_file_paths:
            tmp_dir = 'tmp'
            results_path = os.path.join(self._base.settings.seq_results_directory, f'{self._name}.query_results.m8')
            seq_db_path = os.path.join(self._base.settings.seq_storage_directory, self._name)
            out = self._base._execute_cli('easy-search', dict(
                filenames=[*input_file_paths, seq_db_path, results_path, tmp_dir],
                shuffleDatabase=False,
                sensitivity=5.7,
                removeTmpFiles=False,
                writeLookup=False,
                outfmt="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
                searchType=PARAM_DB_SEARCH_TYPE_MAPPING[search_type],
                alignmentMode=3,
            ))
            return SearchResults(out.blast_tab_records)

    def create_index(self, search_type: str = 'nucleotides') -> IndexStats:
        tmp_dir = 'tmp'
        seq_db_path = os.path.join(self._base.settings.seq_storage_directory, self._name)
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
        return IndexStats(
            db_size=int(out.vars_str['INDEX_TABLE_DB_SIZE']),
            index_entries=int(out.vars_str['INDEX_TABLE_ENTRIES']),
            avg_kmer_size=float(out.vars_str['INDEX_AVG_KMER_SIZE']),
        )


class Databases(MMSeqsBase):
    def __init__(self, *args):
        super().__init__(*args)

    def __getitem__(self, index):
        return self.list()[index]

    def list(self) -> List[Database]:
        """
        List all available databases
        :return: List of databases
        """
        with self.settings.meta_db.open() as meta_db:
            return meta_db.list_get('databases', self)

    def create(self,
               name: str,
               description: str,
               source: GenericSequences,
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
        with GenericSequencesFilePathConverter(source) as input_files:
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
                setattr(new_db, '_base', self)
                meta_db.list_append('databases', new_db)
            return new_db

