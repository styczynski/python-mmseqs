import biosnake

#
# Demonstration of basic biosnake2 operations
#

# Create a client
client = biosnake.Biosnake()

# Create a database from fasta file
# Here we specify name of the database, description and input file
# (The input can also be a Seq/SeqRecord list/iterator/etc.)
client.databases.create("test", "Test database", "a.fasta")

# Get description of the database
print(client.databases[0].description)

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
    headers=["query_sequence_id",
      "target_sequence_id",
      "sequence_identity",
      "target_sequence_aligned"],
)

# Load queries from file:
# results = client.databases[0].search_file("input.fasta", search_type="nucleotides")

# You can pass list of headers to get:
#   query_sequence_id
#   target_sequence_id
#   query_sequence_content
#   target_sequence_content
#   sequence_identity
#   alignment_length
#   number_of_mismatches
#   number_of_gap_openings
#   domain_start_index_query
#   domain_end_index_query
#   domain_start_index_target
#   domain_end_index_target
#   e_value
#   bit_score
# For example:
# results2 = client.databases[0].search(
#     [
#         "ACTAGCTCAGTCAACTAGCTCAGTCCTCAGTCAACTAGCTCAGTCTATATATATACAAC",
#         "ACTAGCTCAGTCAACTAGCTCAGTCCTCAGTCAACTAGCT",
#         "ACTAGCTCAGTCAACTAGCT",
#         "ACTAGCTCAGT",
#     ],
#     search_type="nucleotides",
#     headers=["query_sequence_id", "target_sequence_id", "sequence_identity", "alignment_length", "number_of_mismatches"]
# )

# results.records is a list of lists. Each item contains alignments for each query.
# Each list of alignments consists of single result
# print(results.records)

# You can also get a pandas dataframe
#print(results.dataframe)
results.to_fasta('search_results_{query_sequence_id}.fasta', "{query_sequence_id}", "{target_sequence_aligned}")
