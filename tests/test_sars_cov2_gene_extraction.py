"""
Module providing CLI utility to extract reference genes from input sequences.
"""
from __future__ import with_statement
import os
import click
from collections import defaultdict
import subprocess
import shutil
import pandas as pd


import mmap
from typing import Optional


TEST_REFERENCE_GENES_FILE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "data", "test_sars_cov2_gene_extraction", "reference_genes.fasta")
)

TEST_GENOMES_FILE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "data", "test_sars_cov2_gene_extraction", "genomes.fasta")
)

HEADER_MAPPINGS = {
    "query_sequence_id": 'query',
    "target_sequence_id": 'target',
    "target_sequence_content": 'tseq',
    "query_sequence_length": 'qlen',
    "domain_start_index_query": 'qstart',
    "domain_end_index_query": 'qend',
    "domain_start_index_target": 'tstart',
    "domain_end_index_target": 'tend',
}


def get_file_lines_count(filename: str, text: Optional[str] = None) -> int:
    """
    Count lines in file.
    This is the fastest implementation possible in Python using mmap.
    Optionally function takes a text parameter.
    If that parameter is not None, then only lines containing this text are counted.
    :param filename: Path to the input file
    :param text: Optional text to match in counted lines
    :return: Number of counted lines
    """
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    if text is None:
        while readline():
            lines += 1
    else:
        text = text.encode()
        for line in buf:
            if text in line:
                lines += 1
    return lines


def get_fasta_size(filename: str) -> int:
    """
    Fastest way to determine number of sequences in FASTA file.
    :param filename: Path to the input file
    :return: Number of sequences in FASTA file
    """
    return get_file_lines_count(filename, ">")


def extract_genes(
    gene_fasta: str,
    input_file: str,
    use_native: bool = True,
    remove_temp_files: bool = True,
):
    """
    For given FASTA input file and set of reference genes extract them into separate files.
    :param gene_fasta: Path to input FASTA file with reference genes
    :param input_file: Path to input FASTA file with sequences
    :param use_native: If true will use native bindings. Spawn process otherwise.
    :param remove_temp_files: Remove temporary results?
    """
    sensitivity = 6
    max_sequence_length = 30000
    df = None
    if use_native:
        import biosnake
        client = biosnake.Biosnake()
        input_sequences = client.databases.create(
            "input_sequences", "Input sequences to cut", input_file
        )
        print("OK?")
        df = input_sequences.search_file(
            gene_fasta,
            sensitivity=sensitivity,
            max_sequence_length=max_sequence_length,
            max_results_count_per_query=input_sequences.size,
            search_type="nucleotides",
            headers=[
                "query_sequence_id",
                "target_sequence_id",
                "target_sequence_content",
                "query_sequence_length",
                "domain_start_index_query",
                "domain_end_index_query",
                "domain_start_index_target",
                "domain_end_index_target"
            ],
        ).dataframe
    else:
        headers = list(HEADER_MAPPINGS.keys())
        args = [
            "mmseqs",
            "easy-search",
            gene_fasta,
            input_file,
            "result.m8",
            "biosnake_storage",
            "--search-type",
            "3",
            "--format-output",
            ",".join([HEADER_MAPPINGS[header] for header in headers]),
            "-s",
            str(sensitivity),
            "--max-seq-len",
            str(max_sequence_length),
            "--max-seqs",
            str(get_fasta_size(input_file)),
            "--split-memory-limit",
            "45000",
        ]
        # Spawn mmseqs
        p = subprocess.Popen(" ".join(args), shell=True)
        p.communicate()
        if p.returncode != 0:
            raise Exception("Error running mmseq. Is it installed?")

        # Convert blast M8 format to a dataframe for easier use
        df = pd.read_table("result.m8", names=headers, header=None)

    # Remove temporary files
    if remove_temp_files:
        os.unlink("result.m8")
    shutil.rmtree("biosnake_storage")

    # Process the dataframe.
    # This loop resizes all matches so that all genes have the same size
    # as the reference ones
    df = df.assign(start=0, end=0)
    for index, row in df.iterrows():
        # If domain_start_index_query > 0 then we will take more nucleotides at the start
        sequence_start = (
                row["domain_start_index_target"] - row["domain_start_index_query"]
        )
        # We should take (query_sequence_length-1)-domain_end_index_query additional elements
        # at the match end (query and target are 0-indexed)
        # However Python slice notation list[:index] extracts all items except single one
        # under last index, so we discard -1 here.
        sequence_end = (
                row["domain_end_index_target"]
                + row["query_sequence_length"]
                - row["domain_end_index_query"]
        )
        # Extract the sequence content
        df.at[index, "start"] = sequence_start
        df.at[index, "end"] = sequence_end
        #df.at[index, "alignment"] = row["target_sequence_content"][
        #                            sequence_start:sequence_end
        #                            ]
    return df[["query_sequence_id", "target_sequence_id", "start", "end"]].rename(
        columns=dict(
            query_sequence_id="gene",
            target_sequence_id="id",
        )
    )


def assert_extracted_genes_are_the_same(our: pd.DataFrame, ref: pd.DataFrame):
    def clean_dataframe(input_df: pd.DataFrame) -> pd.DataFrame:
        df = input_df.copy()
        df.drop_duplicates(subset=['gene', 'id', 'start', 'end'], inplace=True)
        df["key"] = df["gene"] + df["id"]
        df.set_index("key", inplace=True)
        return df
    our, ref = clean_dataframe(our), clean_dataframe(ref)
    pd.testing.assert_frame_equal(our, ref, check_like=True)


def test_sars_cov2_gene_extraction():
    biosnake_result, mmseqs_result = tuple(extract_genes(
        TEST_REFERENCE_GENES_FILE,
        TEST_GENOMES_FILE,
        use_native=use_native,
    ) for use_native in [True, False])
    assert_extracted_genes_are_the_same(biosnake_result, mmseqs_result)


if __name__ == "__main__":
    test_sars_cov2_gene_extraction()
