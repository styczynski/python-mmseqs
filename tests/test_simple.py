import biosnake

def test_empty():
    client = biosnake.Biosnake()
    client.databases.create("test", "Test database", "./example/a.fasta")
