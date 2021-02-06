# MMseqs2 bindings for Python

This project provides bidings for mmseqs.
```python
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
print(results.records)
```