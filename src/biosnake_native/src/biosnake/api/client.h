#ifndef API_CLIENT_H
#define API_CLIENT_H

#include <biosnake/api/parameters.h>
#include <state/state.h>

#include <iterator>
#include <cstddef>

class Database;

class Client {
public:

    Client(const std::string& storage_path, const std::string& version);

    void prepare_to_execute_command() const;

    Database get(int index);

    std::vector<Database> list();
    int database_count() const;

    Database create(
        std::string name,
        std::string description,
        std::string input_fasta,
        std::string mode = "copy",
        std::string database_type = DEFAULT_SEARCH_TYPE,
        int offset = 0,
        bool shuffle = false
    );

    ~Client();

    std::string get_storage_path(std::string path) const;
    std::string get_workdir_path() const;
    std::string get_results_path(std::string path) const;

    State state;

    struct DatabasesIterator
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Database;
        using pointer           = int;
        using reference         = Database;

        DatabasesIterator(Client* client, pointer ptr);

        reference operator*();
        pointer operator->();
        DatabasesIterator& operator++();
        DatabasesIterator operator++(int);
        friend bool operator== (const DatabasesIterator& a, const DatabasesIterator& b);
        friend bool operator!= (const DatabasesIterator& a, const DatabasesIterator& b);

    private:
        Client* _client;
        pointer m_ptr;
    };

    DatabasesIterator begin();
    DatabasesIterator end();

private:
    std::string _storage_path, _storage_path_databases, _storage_path_results, _storage_path_workdir;
    std::string _config_path;
};

#endif