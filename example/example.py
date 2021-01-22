import mmseqs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#
# Demonstration of basic mmseqs2 operations
#

# Create a client
client = mmseqs.MMSeqs()

# Create a database from fasta file
# Here we specify name of the database, description and input file
# (The input can also be a Seq/SeqRecord list/iterator/etc.)
client.databases.create("test", "Test database", "DB.fasta")

# We can access stored records as usual, getting SeqRecords back
for record in client.databases[0].records:
    print(record)

# We can also dump a database back to FASTA
client.databases[0].to_fasta("dump.fasta")

# Demonstration of basic database manipulations:
# List all databases
print(client.databases)

# Clone a database
cloned_db = client.databases[0].copy('new_db')

# Remove the clone
cloned_db.remove()

# Create an index for faster searching
client.databases[0].create_index()

# Perform search on a database
# Note that the search queries can be a string with a patch to the FASTA file with queries
results = client.databases[0].search([
    SeqRecord(Seq("ACACACAAAACTACACACAAAACTACACACAAAAAAA"), id="q1"),
    SeqRecord(Seq("TACACACAAAACTAACAAAACTAACAAAACTAACAAAACTTTAGAA"), id="q2"),
    SeqRecord(Seq("CAAAACTAACAAAACTACACAAAACTACACACGGGAGAGCGA"), id="q3"),
    SeqRecord(Seq("CACGGGAGAGAAAACACGTAACTCTAGTCTAGTAACAAACAAAAACAAAACGTAACTCTAGTCTAGTA"), id="q4"),
], search_type='nucleotides')

# results.records is a list of lists. Each item contains alignments for each query.
# Each list of alignments consists of single result
print(results.records[0])

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
print(results.dataframe)

