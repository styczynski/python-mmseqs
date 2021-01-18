from typing import Sequence, List, Optional
import os
from dataclasses import dataclass
from glob import glob
from datetime import date

from .utils import remove_paths
from .base import MMSeqsBase

PARAM_DB_TYPE_MAPPING = dict(
    auto=0,
    amino_acid=1,
    nucleotides=2,
)

PARAM_CREATEDB_MODE_MAPPING = dict(
    copy=0,
    soft_link=1,
)

# INPUT="$1"
# if [ -n "$TRANSLATED" ]; then
#     # 1. extract orf
#     if notExists "$2/orfs_aa.dbtype"; then
#         # shellcheck disable=SC2086
#         "$MMSEQS" extractorfs "$INPUT" "$2/orfs_aa" ${ORF_PAR} \
#             || fail "extractorfs died"
#     fi
#
#     # shellcheck disable=SC2086
#     "$MMSEQS" $INDEXER "$2/orfs_aa" "$INPUT" ${INDEX_PAR} \
#         || fail "indexdb died"
#
#     if [ -n "$REMOVE_TMP" ]; then
#         # shellcheck disable=SC2086
#         "$MMSEQS" rmdb "$2/orfs_aa" ${VERBOSITY}
#         rm -f "$2/createindex.sh"
#     fi
# elif [ -n "$LIN_NUCL" ] || [ -n "$NUCL" ]; then
#       # 1. extract orf
#     if notExists "$2/nucl_split_seq.dbtype"; then
#         # shellcheck disable=SC2086
#         "$MMSEQS" splitsequence "$INPUT" "$2/nucl_split_seq" ${SPLIT_SEQ_PAR} \
#             || fail "splitsequence died"
#     fi
#
#     # shellcheck disable=SC2086
#     "$MMSEQS" $INDEXER "$2/nucl_split_seq" "$INPUT" ${INDEX_PAR} \
#         || fail "indexdb died"
#
#     if [ -n "$REMOVE_TMP" ]; then
#         # shellcheck disable=SC2086
#         "$MMSEQS" rmdb "$2/nucl_split_seq" ${VERBOSITY}
#         rm -f "$2/createindex.sh"
#     fi
# else
#     # shellcheck disable=SC2086
#     "$MMSEQS" $INDEXER "$INPUT" "$INPUT" ${INDEX_PAR} \
#         || fail "indexdb died"
# fi



@dataclass
class Database:
    """Class for keeping track of an item in inventory."""
    name: str
    input_files: List[str]
    database_type: str
    created_on: date

    def delete(self):
        with self._base.settings.meta_db.open() as meta_db:
            removed_dbs = meta_db.list_filter('databases', self._base, lambda db: db.name != self.name)
            for db in removed_dbs:
                remove_paths([
                    f'{db.name}',
                    f'{db.name}.*',
                    f'{db.name}_h',
                    f'{db.name}_h.*',
                ], base_path=self._base.settings.seq_storage_directory, is_glob=True)

    def create_index(self):
        tmp_dir = 'tmp'
        seq_db_path = os.path.join(self._base.settings.seq_storage_directory, self.name)
        if not os.path.isfile(f'{seq_db_path}.dbtype'):
            raise Exception('Invalid database state')
        if not os.path.isdir(tmp_dir):
            raise Exception('Missing temporary directory')


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
               input_files: Sequence[str],
               mode: str = "copy",
               database_type: str = "auto",
               offset: int = 0,
               shuffle: bool = False,
    ) -> Database:
        """
        Create a sequence database from multiple FASTA file

        :param name: Name of the new database
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
        self._execute_cli([
            'createdb', *input_files, os.path.join(self.settings.seq_storage_directory, name),
            '--shuffle', int(shuffle),
            '--createdb-mode', PARAM_CREATEDB_MODE_MAPPING[mode],
            '--dbtype', PARAM_DB_TYPE_MAPPING[database_type],
            '--id-offset', offset
        ])
        with self.settings.meta_db.open() as meta_db:
            new_db = Database(
                name=name,
                input_files=input_files,
                created_on=date.today(),
                database_type=database_type,
            )
            meta_db.list_append('databases', new_db)
        return new_db

