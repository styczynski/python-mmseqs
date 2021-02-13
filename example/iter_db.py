import biosnake

client = biosnake.Biosnake()
client.databases.create("test", "Test database", "example/a.fasta")

print(client.databases.dataframe)
    # for record in db:
    #     print(record.seq.tostring())