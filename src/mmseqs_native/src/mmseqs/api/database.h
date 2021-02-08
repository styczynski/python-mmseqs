#ifndef API_DATABASE_H
#define API_DATABASE_H

#include <mmseqs/api/client.h>
#include <mmseqs/api/search_results.h>
#include <mmseqs/api/parameters.h>
#include <mmseqs/api/python_utils.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/dBReader.h>
#include <memory>

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

    struct Record
    {
        Record(const char* sequence, const size_t sequence_length, const char* header, const size_t header_length):
            _sequence(py::array(py::cast(std::string(sequence, sequence_length)))),
            _header(header, header_length) {};
        Record() {}

        py::array _sequence;
        std::string _header;
    };

    struct DBIterator
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Record;
        using pointer           = size_t;
        using reference         = Record;

        DBIterator(std::shared_ptr<DBReader<unsigned int>> db_reader, std::shared_ptr<DBReader<unsigned int>> header_reader, pointer ptr);

        reference operator*();
        pointer operator->();
        DBIterator& operator++();
        DBIterator operator++(int);
        friend bool operator== (const DBIterator& a, const DBIterator& b);
        friend bool operator!= (const DBIterator& a, const DBIterator& b);

    private:
        std::shared_ptr<DBReader<unsigned int>> _db_reader;
        std::shared_ptr<DBReader<unsigned int>> _header_reader;
        pointer m_ptr;
    };

    std::shared_ptr<DBReader<unsigned int>> db_reader;
    std::shared_ptr<DBReader<unsigned int>> header_reader;

    DBIterator begin();
    DBIterator end();
    void _init_readers();

    const std::tuple<py::array, py::array, py::array> getColumnData();

private:
    StateDatabase _state;
    Client* _parent;
    mmseqs_output _out;
};

#endif
