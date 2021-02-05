import collections
from tempfile import NamedTemporaryFile
from typing import Iterator, Sequence, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Generic sequence that can be anything that is string or biopython sequence-like
GenericSequence = Union[Seq, SeqRecord, str]
# List of generic sequences or a single one
GenericSequences = Union[
    GenericSequence, Sequence[GenericSequence], Iterator[GenericSequence]
]


def normalize_single_sequence(
    sequence: GenericSequence,
    i: int,
) -> SeqRecord:
    """
    Take a single sequence and normalize it's form to Biopython SeqRecord object.
    If the sequence has no ID the function uses i parameter to assign numeric ID to the sequence record.
    Supported input sequence types are:
        - Seq object
        - SeqRecord object
        - string
    :param sequence: Input sequence of supported type
    :param i: Index of that sequence
    :return: SeqRecord object
    """
    if isinstance(sequence, SeqRecord):
        return sequence
    elif isinstance(sequence, Seq):
        return SeqRecord(sequence, str(i))
    elif isinstance(sequence, str):
        return SeqRecord(Seq(sequence), str(i))


def convert_to_biopython_sequences(
    sequences: GenericSequences,
) -> Sequence[SeqRecord]:
    """
    Take a bunch of sequences and normalize their form to list of Biopython SeqRecord objects.
    Supported input sequences types are:
        - List of Seq objects (*)
        - List of SeqRecord objects
        - List of strings (*)
        - Seq object (**)
        - SeqRecord object (**)
        - string (**)
    If the input sequences does not provide information about the sequence id (marked with asterisks *) then
    the list index (starting at 0) will be used for sequence ID.
    If the input sequences are in fact only a single object then the integer value of 0 will be used for sequence ID.
    (marked with double asterisks **)
    :param sequences: Some input of supported type
    :return: List of Biopython SeqRecord objects
    """
    if isinstance(sequences, SeqRecord):
        return [sequences]
    elif isinstance(sequences, Seq):
        return [SeqRecord(sequences, 0)]
    elif isinstance(sequences, str):
        return [SeqRecord(Seq(sequences), 0)]
    elif isinstance(sequences, collections.Iterable):
        return [
            normalize_single_sequence(seq, i)
            for i, seq in enumerate(sequences)
        ]
    raise Exception("Invalid input type for normalize_sequences() was given.")


class GenericSequencesFilePathConverter:
    def __init__(self, input: GenericSequences):
        self.input = input
        self.tmp_file = None

    def __enter__(self):
        if isinstance(self.input, str):
            return [self.input]
        else:
            self.tmp_file = NamedTemporaryFile(suffix=".fasta", mode="w")
            output_handle = self.tmp_file.__enter__()
            output_handle.seek(0)
            records = convert_to_biopython_sequences(self.input)
            SeqIO.write(records, output_handle, "fasta")
            output_handle.seek(0)
            return [self.tmp_file.name]

    def __exit__(self, type, value, traceback):
        if self.tmp_file is not None:
            self.tmp_file.__exit__(type, value, traceback)
