import biosnake
import os

TEST_INPUT_FILE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "data", "test_basic_operations", "input.fasta")
)

def test_basic_search():
    client = biosnake.Biosnake()
    db = client.databases.create("test", "Test database", TEST_INPUT_FILE)
    assert db.description == "Test database"
    assert db.name == "test"
    search_results = db.search(
        [
            "ACTAGCTCAGTCAACTAGCTCAGTCCTCAGTCAACTAGCTCAGTCTATATATATACAAC",
            "ACTAGCTCAGTCAACTAGCTCAGTCCTCAGTCAACTAGCT",
            "ACTAGCTCAGTCAACTAGCT",
            "ACTAGCTCAGT",
        ],
        search_type="nucleotides",
        headers=["query_sequence_id",
          "target_sequence_id",
          "sequence_identity",
          "target_sequence_aligned"],
    )
    search_results.dataframe.to_csv('hujemuje.csv')

if __name__ == '__main__':
    test_basic_search()