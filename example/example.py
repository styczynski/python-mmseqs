import mmseqs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

client = mmseqs.MMSeqs()
client.databases.create("test", "Test database", "DB.fasta")

for record in client.databases[0].records:
    pass
client.databases[0].to_fasta("dump.fasta")

res = client.databases[0].search([
    SeqRecord(Seq("ACACACAAAACTACACACAAAACTACACACAAAAAAA"), id="q1"),
    SeqRecord(Seq("TACACACAAAACTAACAAAACTAACAAAACTAACAAAACTTTAGAA"), id="q2"),
    SeqRecord(Seq("CAAAACTAACAAAACTACACAAAACTACACACGGGAGAGCGA"), id="q3"),
    SeqRecord(Seq("CACGGGAGAGAAAACACGTAACTCTAGTCTAGTAACAAACAAAAACAAAACGTAACTCTAGTCTAGTA"), id="q4"),
], search_type='nucleotides')

print(res.records[0])

#
# print(client.databases[0].search([
#   SeqRecord(Seq("ACACACAAAACTACACACAAAACTACACACAAAAAAA"), id="q1"),
#   SeqRecord(Seq("TACACACAAAACTAACAAAACTAACAAAACTAACAAAACTTTAGAA"), id="q2"),
#   SeqRecord(Seq("CAAAACTAACAAAACTACACAAAACTACACACGGGAGAGCGA"), id="q3"),
#   SeqRecord(Seq("CACGGGAGAGAAAACACGTAACTCTAGTCTAGTAACAAACAAAAACAAAACGTAACTCTAGTCTAGTA"), id="q4"),
# ], search_type='nucleotides'))


# print(client.databases.list())
# client.databases.create("test", "Test database", ["input.fasta"])
# print(client.databases.list())
# new_db = client.databases[0].copy('new_db')
# print(new_db)
# client.databases[0].remove()
# print(client.databases.list())
# client.databases.create("test", ["input.fasta"])
