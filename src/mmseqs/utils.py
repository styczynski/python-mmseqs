import os
from typing import List, Optional
from glob import glob

def to_args(args_str: str) -> List[str]:
   return [param.strip() for param in args_str.split(' ') if len(param.strip()) > 0]


def remove_paths(paths: List[str], base_path: Optional[str] = None, is_glob=False):
    for glob_pattern in paths:
        if base_path is not None:
            glob_pattern = os.path.join(base_path, glob_pattern)
        if is_glob:
            evaluated_paths = glob(glob_pattern)
        else:
            evaluated_paths = [glob_pattern]
        for file_path in evaluated_paths:
            os.unlink(file_path)
