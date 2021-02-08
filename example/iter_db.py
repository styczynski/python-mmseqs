import mmseqs

client = mmseqs.MMSeqs()
client.databases.create("test", "Test database", "example/a.fasta")

print(client.databases.dataframe)
    # for record in db:
    #     print(record.seq.tostring())