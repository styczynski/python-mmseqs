import mmseqs

#
# Demonstration of basic mmseqs2 operations
#

# Create a client
client = mmseqs.MMSeqs()
a = client.databases.create("testa", "Test database A", "a.fasta")
print(a.name)
