
#include <mmseqs/commons/application.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/output.h>
#include <mmseqs/api/client.h>
#include <mmseqs/api/database.h>
#include <mmseqs/commons/fileUtil.h>
#include <state/state.h>

Client::Client(const std::string& storage_path, const std::string& version):
   _storage_path(storage_path),
   _storage_path_results(storage_path + "/results"),
   _storage_path_databases(storage_path + "/databases"),
   _storage_path_workdir(storage_path + "/workdir"),
   _config_path(storage_path + "/config.toml"),
   state(version)
{
    mmseqs_output* out;
    if (!FileUtil::fileExists(out, _config_path.c_str())) {
        prepare_to_execute_command();
        state.dump(std::ofstream(_config_path, std::ios_base::binary));
    }
    std::string config_path = get_storage_path(_config_path);
    state = State::load(std::ifstream(_config_path, std::ios_base::binary), config_path);
}

void Client::prepare_to_execute_command() const {
    mmseqs_output* out;
    if (!FileUtil::directoryExists(out, _storage_path.c_str())) {
        FileUtil::makeDir(out, _storage_path.c_str());
    }
    if (!FileUtil::directoryExists(out, _storage_path_databases.c_str())) {
        FileUtil::makeDir(out, _storage_path_databases.c_str());
    }
    if (!FileUtil::directoryExists(out, _storage_path_workdir.c_str())) {
        FileUtil::makeDir(out, _storage_path_workdir.c_str());
    }
    if (!FileUtil::directoryExists(out, _storage_path_results.c_str())) {
        FileUtil::makeDir(out, _storage_path_results.c_str());
    }
}

Database Client::get(int index) {
    return Database(state.databases.find(state.databases_ids[index])->second, this);
}

int Client::database_count() const {
    return state.databases.size();
}

std::vector<Database> Client::list() {
    std::vector<Database> db_list;
    for (int i=0; i<state.databases_ids.size(); ++i) {
        db_list.push_back(get(i));
    }
    return db_list;
}

Database Client::create(
    std::string name,
    std::string description,
    std::string input_fasta,
    std::string mode,
    std::string database_type,
    int offset,
    bool shuffle
) {
    prepare_to_execute_command();
    Parameters args;
    args.baseTmpPath = get_workdir_path();
    args.filenames = {
        input_fasta,
        get_storage_path(name),
    };
    args.identifierOffset = offset;
    args.dbType = PARAM_DB_TYPE_MAPPING[database_type];
    args.createdbMode = PARAM_CREATEDB_MODE_MAPPING[mode];
    args.shuffleDatabase = (int)shuffle;
    auto out = call_mmseqs("createdb", args);

    auto db_state = state.create_database(name, description, database_type);
    return Database(db_state, this);
}

Client::~Client() {
    state.dump(std::ofstream(_config_path, std::ios_base::binary));
}

std::string Client::get_storage_path(std::string path) const {
    return _storage_path_databases + "/" + path;
}

std::string Client::get_workdir_path() const {
    return _storage_path_workdir + "/";
}

std::string Client::get_results_path(std::string path) const {
    return _storage_path_results + "/" + path;
}

Client::DatabasesIterator::DatabasesIterator(Client* client, pointer ptr) : m_ptr(ptr), _client(client) {}
Client::DatabasesIterator::reference Client::DatabasesIterator::operator*() { return _client->get(m_ptr); }
Client::DatabasesIterator::pointer Client::DatabasesIterator::operator->() { return m_ptr; }
Client::DatabasesIterator& Client::DatabasesIterator::operator++() { m_ptr++; return *this; }
Client::DatabasesIterator Client::DatabasesIterator::operator++(int) { Client::DatabasesIterator tmp = *this; ++(*this); return tmp; }
bool operator== (const Client::DatabasesIterator& a, const Client::DatabasesIterator& b) { return a.m_ptr == b.m_ptr; };
bool operator!= (const Client::DatabasesIterator& a, const Client::DatabasesIterator& b) { return a.m_ptr != b.m_ptr; };
Client::DatabasesIterator Client::begin() { return Client::DatabasesIterator(this, 0); }
Client::DatabasesIterator Client::end() { return Client::DatabasesIterator(this, database_count()); }
