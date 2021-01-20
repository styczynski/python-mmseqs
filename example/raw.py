import mmseqs

client = mmseqs.MMSeqs()
# client.databases.create("test", "Test database", ["input.fasta"])
# client.databases[0].create_index()
client.databases[0].search('query.fasta')


# print(client.databases.list())
# client.databases.create("test", "Test database", ["input.fasta"])
# print(client.databases.list())
# new_db = client.databases[0].copy('new_db')
# print(new_db)
# client.databases[0].remove()
# print(client.databases.list())
#client.databases.create("test", ["input.fasta"])