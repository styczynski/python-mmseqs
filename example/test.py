import mmseqs

# Utility modules to manipulate DBs
#   view                  Print DB entries given in --id-list to stdout
#   apply                 Execute given program on each DB entry
#   filterdb              DB filtering by given conditions
#   swapdb                Transpose DB with integer values in first column
#   prefixid              For each entry in a DB prepend the entry key to the entry itself
#   suffixid              For each entry in a DB append the entry key to the entry itself
#   renamedbkeys          Create a new DB with original keys renamed


# Input database creation
#   databases             List and download databases
#   createdb              Convert FASTA/Q file(s) to a sequence DB
#   createindex           Store precomputed index on disk to reduce search overhead
#   createlinindex        Create linsearch index
#   convertmsa            Convert Stockholm/PFAM MSA file to a MSA DB
#   tsv2db                Convert a TSV file to any DB
#   tar2db                Convert content of tar archives to any DB
#   msa2profile           Convert a MSA DB to a profile DB

if True:
    mmseqs.MMSeqs().execute(["createindex", "mmseqs_storage/databases/test", "tmp", "--search-type", "3"])
    # mmseqs.MMSeqs().execute(["databases", "-h"])
else:
    client = mmseqs.MMSeqs()
    print(client.databases.list())
    client.databases.create("test", ["input.fasta"])
    print(client.databases.list())
    #client.databases[0].delete()
    #print(client.databases.list())
