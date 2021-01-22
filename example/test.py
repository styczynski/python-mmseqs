import mmseqs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#
# Demonstration of basic mmseqs2 operations
#

# Create a client
client = mmseqs.MMSeqs()
a = client.databases.create("testa", "Test database A", "a.fasta")
b = client.databases.create("testb", "Test database B", "b.fasta")
print([db.name for db in client.databases])
a.name = 'renamed'
print([db.name for db in client.databases])

# a = client.databases.create("testa", "Test database A", "a.fasta")
# b = client.databases.create("testb", "Test database B", "b.fasta")
#
# a.name = 'testarenamed'
