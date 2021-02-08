#ifndef API_DATABASE_H
#define API_DATABASE_H

#include <mmseqs/api/client.h>
#include <mmseqs/api/search_results.h>
#include <mmseqs/api/parameters.h>
#include <mmseqs/output.h>

class Database {
public:
    Database();
    explicit Database(StateDatabase state, Client* parent);

    std::string to_fasta(std::string output_path = "");

    Database copy(std::string db_copy_name = "");

    std::string remove();

    SearchResults search(
        std::vector<std::string> sequences,
        std::string search_type = DEFAULT_SEARCH_TYPE,
        std::vector<std::string> headers = {},
        float sensitivity = DEFAULT_SENSITIVITY,
        int max_sequence_length = DEFAULT_MAX_SEQUENCE_LENGTH
    );

    SearchResults search_file(
        std::string search_input_fasta,
        std::string search_type = DEFAULT_SEARCH_TYPE,
        std::vector<std::string> headers = {},
        float sensitivity = DEFAULT_SENSITIVITY,
        int max_sequence_length = DEFAULT_MAX_SEQUENCE_LENGTH
    );

    mmseqs_output create_index(
        std::string search_type = DEFAULT_SEARCH_TYPE,
        float sensitivity = DEFAULT_SENSITIVITY,
        int max_sequence_length = DEFAULT_MAX_SEQUENCE_LENGTH,
        int max_orf_length = DEFAULT_MAX_ORF_LENGTH,
        int min_orf_length = DEFAULT_MIN_ORF_LENGTH,
        int orf_start_mode = DEFAULT_ORF_START_MODE
    );

    void setName(const std::string& new_name);

    const std::string& getName() const;

    void setDescription(const std::string& new_description);

    const std::string& getDescription() const;

    const std::string& getType() const;

private:
    StateDatabase _state;
    Client* _parent;
};

#endif
