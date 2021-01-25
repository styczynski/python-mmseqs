import mmseqs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#
# Demonstration of basic mmseqs2 operations
#

# Create a client
client = mmseqs.MMSeqs()
a = client.databases.create("testa", "Test database A", "a.fasta")
print(a.name)