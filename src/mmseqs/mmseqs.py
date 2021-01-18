import os
import pathlib
from typing import Optional, Sequence, Tuple

from .mmseqs_native import MMSeqsCallArgs, _call_mmseqs


class MMSeqs:
    def __init__(self):
        pass

    def execute(self, cli_args):
        call_args = MMSeqsCallArgs()
        call_args.cli_args = cli_args
        return _call_mmseqs(call_args)