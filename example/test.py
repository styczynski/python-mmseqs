import biosnake

#
# Demonstration of basic biosnake2 operations
#

# Create a client
client = biosnake.Biosnake()
a = client.databases.create("testa", "Test database A", "a.fasta")
print(a.size)
