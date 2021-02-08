import mmseqs


client = mmseqs.MMSeqs()
client.databases.create("test", "Test database", "example/a.fasta")

for db in client.databases:
    print(db.name)