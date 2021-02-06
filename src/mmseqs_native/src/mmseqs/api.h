#ifndef API_H
#define API_H

#include <state/state.h>
#include <mmseqs/rand.h>
#include <mmseqs/commons/fileUtil.h>

static std::map<std::string, int> PARAM_DB_TYPE_MAPPING = {
    {"auto", 0},
    {"amino_acid", 1},
    {"nucleotides", 2},
};

static std::map<std::string, int> PARAM_CREATEDB_MODE_MAPPING = {
    {"copy", 0},
    {"soft_link", 1},
};

static std::map<std::string, int> PARAM_DB_SEARCH_TYPE_MAPPING = {
    {"auto", 0},
    {"protein", 1},
    {"translated", 2},
    {"nucleotides", 3},
    {"translated_nucleotides_aligned", 4},
};

class Database;

class SearchResults {
public:

    SearchResults() {}
    explicit SearchResults(const std::vector<std::vector<mmseqs_blast_tab_record>>& records): _records(records) {}

    const std::vector<std::vector<mmseqs_blast_tab_record>> &getRecords() const {
        return _records;
    }

private:
    std::vector<std::vector<mmseqs_blast_tab_record>> _records;
};

class Databases {
public:

    Databases(const std::string& storage_path, const std::string& version);

    void prepare_to_execute_command();

    Database get(int index);

    std::vector<Database> list();

    Database create(
        std::string name,
        std::string description,
        std::string input_fasta,
        std::string mode = "copy",
        std::string database_type = "auto",
        int offset = 0,
        bool shuffle = false
    );

    ~Databases();

    std::string get_storage_path(std::string path);
    std::string get_workdir_path();
    std::string get_results_path(std::string path);

    State state;

private:
    std::string _storage_path, _storage_path_databases, _storage_path_results, _storage_path_workdir;
    std::string _config_path;
};

std::string write_temp_fasta(std::vector<std::string> sequences, std::string temp_path) {
    std::string file_path = temp_path + get_uuid() + "query_inputs.fasta";
    std::ofstream out(file_path, std::ios_base::binary);
    for (int i=0;i<sequences.size();++i) {
        out << ">query" + SSTR(i) + "\n" + sequences[i] + "\n";
    }
    return file_path;
}


class Database {
public:
    Database() {}
    explicit Database(StateDatabase state, Databases* parent): _state(state), _parent(parent) {}

    std::string to_fasta(std::string output_path = "") {
        _parent->prepare_to_execute_command();

        if (output_path.size() == 0) {
            output_path = _state.name + "_output.fasta";
        }

        std::string db_path = _parent->get_storage_path(_state.name);

        Parameters args;
        args.baseTmpPath = _parent->get_workdir_path();
        args.setDBFields(1, db_path),
        args.setDBFields(2, output_path);
        args.filenames = {db_path, output_path};
        call_mmseqs("convert2fasta", args);
        return output_path;
    }

    Database copy(std::string db_copy_name = "") {
        _parent->prepare_to_execute_command();

        if (db_copy_name.size() == 0) {
            int no = 0;
            while (no<5000) {
                mmseqs_output* out;
                if (!FileUtil::fileExists(out, (_parent->get_storage_path(_state.name + "_copy" + SSTR(no))).c_str())) {
                    db_copy_name = _state.name + "_copy" + SSTR(no);
                    break;
                }
                no++;
            }
        }

        std::string db_path = _parent->get_storage_path(_state.name);
        std::string db_copy_path = _parent->get_storage_path(db_copy_name);
        
        Parameters args;
        args.baseTmpPath = _parent->get_workdir_path();
        args.setDBFields(1, db_path);
        args.setDBFields(2, db_copy_path);
        call_mmseqs("cpdb", args);

        auto copy_state = _parent->state.create_database(db_copy_name, _state.description, _state.database_type);
        return Database(copy_state, _parent);
    }

    std::string remove() {
        _parent->prepare_to_execute_command();
        std::string seq_db_path = _parent->get_storage_path(_state.name);
        std::string seq_db_path_h = _parent->get_storage_path(_state.name + "_h");

        Parameters args;
        args.baseTmpPath = _parent->get_workdir_path();
        args.setDBFields(1, seq_db_path);
        call_mmseqs("rmdb", args);

        args.setDBFields(1, seq_db_path_h);
        call_mmseqs("rmdb", args);

        _parent->state.remove_database(_state.name);
        return _state.name;
    }

    SearchResults search(
        std::vector<std::string> sequences,
        std::string search_type = "auto"
    ) {
        auto path = write_temp_fasta(sequences, _parent->get_workdir_path());
        auto search_results = search_file(path, search_type);
        mmseqs_output* out;
        FileUtil::remove(out, path.c_str());
        return search_results;
    }

    SearchResults search_file(
        std::string search_input_fasta = "nucleotides",
        std::string search_type = "auto"
    ) {
        _parent->prepare_to_execute_command();
        std::string tmp_dir = "tmp_" + get_uuid();
        std::string results_path = _parent->get_results_path(_state.name + get_uuid() + ".query_results.m8");
        std::string seq_db_path = _parent->get_storage_path(_state.name);

        Parameters args;
        args.baseTmpPath = _parent->get_workdir_path();
        args.filenames = {
            search_input_fasta,
            seq_db_path,
            results_path,
            tmp_dir,
        };
        args.shuffleDatabase=false;
        args.sensitivity=5.7;
        args.removeTmpFiles=false;
        args.writeLookup=false;
        args.outfmt="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits";
        args.searchType=PARAM_DB_SEARCH_TYPE_MAPPING[search_type];
        args.alignmentMode=3;

        auto search_output = call_mmseqs("easy-search", args);

        mmseqs_output* out;
        if (FileUtil::fileExists(out, results_path.c_str())) {
            FileUtil::remove(out, results_path.c_str());
        }
        if (FileUtil::fileExists(out, tmp_dir.c_str())) {
            FileUtil::remove(out, tmp_dir.c_str());
        }

        return SearchResults(search_output.blast_tab_records);
    }

    mmseqs_output create_index(std::string search_type = "nucleotides") {
        _parent->prepare_to_execute_command();
        Parameters args;
        args.baseTmpPath = _parent->get_workdir_path();

        std::string seq_db_path = _parent->get_storage_path(_state.name);
        std::string tmp_dir = "tmp_dir";

        args.setSeedSubstitutionMatrices("blosum62.out", "nucleotide.out");
        args.setDBFields(1, seq_db_path);
        args.setDBFields(2, tmp_dir);
        args.filenames = {tmp_dir};
        args.searchType = PARAM_DB_SEARCH_TYPE_MAPPING[search_type];
        args.orfStartMode=1;
        args.orfMinLength=30;
        args.orfMaxLength=32734;
        args.kmerScore=0;
        args.maskMode=1;
        args.sensitivity=7.5;
        args.removeTmpFiles=true;

        auto out = call_mmseqs("createindex", args);
        return out;
    }

    void setName(const std::string &name_) {

    }

    const std::string &getName() const {
        return _state.name;
    }
    
    void setDescription(const std::string &name_) {

    }

    const std::string &getDescription() const {
        return _state.description;
    }

    const std::string &getType() const {
        return _state.database_type;
    }

private:
    StateDatabase _state;
    Databases* _parent;
};

Databases::Databases(const std::string& storage_path, const std::string& version):
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

void Databases::prepare_to_execute_command() {
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

Database Databases::get(int index) {
    return Database(state.databases[state.databases_ids[index]], this);
}

std::vector<Database> Databases::list() {
    std::vector<Database> db_list;
    for (int i=0; i<state.databases_ids.size(); ++i) {
        db_list.push_back(get(i));
    }
    return db_list;
}

Database Databases::create(
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

Databases::~Databases() {
    state.dump(std::ofstream(_config_path, std::ios_base::binary));
}

std::string Databases::get_storage_path(std::string path) {
    return _storage_path_databases + "/" + path;
}

std::string Databases::get_workdir_path() {
    return _storage_path_workdir + "/";
}

std::string Databases::get_results_path(std::string path) {
    return _storage_path_results + "/" + path;
}

#endif