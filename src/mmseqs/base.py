from typing import Optional
from sqlitedict import SqliteDict
from dataclasses import dataclass
import copy

class MetaDatabase:
    def __init__(self, path):
        self.path = path

    def open(self):
        return MetaDatabaseConnection(connection=SqliteDict(self.path, autocommit=True))

@dataclass
class MMSeqsSettings:
    """Class for keeping track of an item in inventory."""
    storage_directory: str
    seq_storage_directory: str
    seq_results_directory: str
    tmp_directory: str
    meta_db: Optional[MetaDatabase]


class MMSeqsBase:
    _execute_cli: any
    settings: any

    def __init__(self, execute_cli, settings):
        self._execute_cli = execute_cli
        self.settings = settings

class ObjectWithBaseRef:
    _base: MMSeqsBase

@dataclass
class MetaDatabaseConnection:
    connection: SqliteDict

    def set(self, name, val):
        self.connection[name] = val

    def list_append(self, name, obj):
        list_vals = self.connection[name]
        obj = copy.copy(obj)
        if hasattr(obj, '_base'):
            delattr(obj, '_base')
        list_vals.append(obj)
        self.connection[name] = list_vals

    def list_replace(self, name, base, replacement, filter_fn, save_back=False):
        result_objs = []
        called_once = False
        for obj in self.connection[name]:
            if filter_fn(obj):
                if not called_once:
                    if hasattr(replacement, '_base'):
                        delattr(replacement, '_base')
                    result_objs.append(replacement)
                called_once = True
            else:
                if hasattr(obj, '_base'):
                    delattr(obj, '_base')
                result_objs.append(obj)
        if save_back:
            self.connection[name] = result_objs
        for obj in result_objs:
            setattr(obj, '_base', base)
        return result_objs

    def list_filter(self, name, base, filter_fn, save_back=False):
        removed_objs = []
        left_objs = []
        for obj in self.connection[name]:
            if filter_fn(obj):
                left_objs.append(obj)
            else:
                setattr(obj, '_base', base)
                removed_objs.append(obj)
        if save_back:
            self.connection[name] = left_objs
        return removed_objs, left_objs

    def list_get(self, name, base):
        results = []
        for obj in self.connection[name]:
            setattr(obj, '_base', base)
            results.append(obj)
        return results

    def has(self, name):
        return name in self.connection

    def __enter__(self):
        return MetaDatabaseConnection(self.connection.__enter__())

    def __exit__(self, *args):
        self.connection.__exit__(*args)

