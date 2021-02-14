import biosnake


client = biosnake.Biosnake()
client.databases.create("test", "Test database", "example/a.fasta")

for db in client.databases:
    print(db.name)

v = client.databases[0].test()
print(repr(v))