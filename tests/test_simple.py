import mmseqs

def test_empty():
    client = mmseqs.MMSeqs()
    client.databases.create("test", "Test database", "./example/a.fasta")
