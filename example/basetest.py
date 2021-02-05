from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import mmseqs

#
# Demonstration of basic mmseqs2 operations
#

# Create a client
client = mmseqs.MMSeqs()

# Create a database from fasta file
# Here we specify name of the database, description and input file
# (The input can also be a Seq/SeqRecord list/iterator/etc.)
client.databases.create("test", "Test database", "a.fasta")
print(client.databases[0].description)
import sys
sys.exit(0)

# Perform search on a database
# Note that the search queries can be a string with a patch to the FASTA file with queries
results = client.databases[0].search(
    [
        "ACTAGCTCAGTCAACTAGCTCAGTCCTCAGTCAACTAGCTCAGTCTATATATATACAAC",
        "ACTAGCTCAGTCAACTAGCTCAGTCCTCAGTCAACTAGCT",
        "ACTAGCTCAGTCAACTAGCT",
        "ACTAGCTCAGT",
    ],
    search_type="nucleotides",
)

# results.records is a list of lists. Each item contains alignments for each query.
# Each list of alignments consists of single result
print(results)

# You can also get a pandas dataframe with the following columns:
#  (Compatible with Blast M8 format - http://www.pangloss.com/wiki/Blast)
#     query_sequence_id: str
#     target_sequence_id: str
#     sequence_identity: float
#     alignment_length: int
#     number_of_mismatches: int
#     number_of_gap_openings: int
#     domain_start_index_query: int
#     domain_end_index_query: int
#     domain_start_index_target: int
#     domain_end_index_target: int
#     e_value: float
#     bit_score: int
