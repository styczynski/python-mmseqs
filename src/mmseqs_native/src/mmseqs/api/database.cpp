
#include <mmseqs/commons/application.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/output.h>
#include <mmseqs/api/database.h>
#include <mmseqs/api/rand.h>
#include <mmseqs/api/utils.h>
#include <mmseqs/commons/fileUtil.h>

Database::Database() {}
Database::Database(StateDatabase state, Client* parent): _state(state), _parent(parent) {}

std::string Database::to_fasta(std::string output_path) {
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

Database Database::copy(std::string db_copy_name) {
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

std::string Database::remove() {
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

SearchResults Database::search(
    std::vector<std::string> sequences,
    std::string search_type,
    std::vector<std::string> headers,
    float sensitivity,
    int max_sequence_length
) {
    auto path = write_temp_fasta(sequences, _parent->get_workdir_path());
    auto search_results = search_file(path, search_type, headers, sensitivity, max_sequence_length);
    mmseqs_output* out;
    FileUtil::remove(out, path.c_str());
    return search_results;
}

SearchResults Database::search_file(
    std::string search_input_fasta,
    std::string search_type,
    std::vector<std::string> headers,
    float sensitivity,
    int max_sequence_length
) {
    _parent->prepare_to_execute_command();
    std::string tmp_dir = "tmp_" + get_uuid();
    std::string results_path = _parent->get_results_path(_state.name + get_uuid() + ".query_results.m8");
    std::string seq_db_path = _parent->get_storage_path(_state.name);

    if (headers.size() == 0) {
        headers = PARAM_SEARCH_COL_NAMES_DEFAULT;
    }

    Parameters args;
    args.baseTmpPath = _parent->get_workdir_path();
    args.filenames = {
        search_input_fasta,
        seq_db_path,
        results_path,
        tmp_dir,
    };
    args.shuffleDatabase=false;
    args.sensitivity=sensitivity;
    args.removeTmpFiles=false;
    args.writeLookup=false;
    args.maxSeqLen=max_sequence_length;
    args.outfmt="";
    for (auto& header_name: headers) {
        args.outfmt += PARAM_SEARCH_COL_NAMES_MAPPING[header_name] + ",";
    }
    args.outfmt = args.outfmt.substr(0, args.outfmt.size()-1);

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

    return SearchResults(search_output.blast_tab_records, headers);
}

mmseqs_output Database::create_index(
    std::string search_type,
    float sensitivity,
    int max_sequence_length,
    int max_orf_length,
    int min_orf_length,
    int orf_start_mode
) {
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
    args.orfStartMode=orf_start_mode;
    args.orfMinLength=min_orf_length;
    args.orfMaxLength=(max_sequence_length > max_orf_length) ? max_sequence_length : max_orf_length;
    args.maxSeqLen=max_sequence_length;
    args.kmerScore=0;
    args.maskMode=1;
    args.sensitivity=sensitivity;
    args.removeTmpFiles=true;

    auto out = call_mmseqs("createindex", args);
    return out;
}

void Database::setName(const std::string& new_name) {
    _parent->prepare_to_execute_command();

    std::string old_name = _state.name;
    std::string seq_db_path = _parent->get_storage_path(old_name);
    std::string new_seq_db_path = _parent->get_storage_path(new_name);

    Parameters args;
    args.baseTmpPath = _parent->get_workdir_path();
    args.filenames = {seq_db_path, new_seq_db_path};
    args.setDBFields(1, seq_db_path);
    args.setDBFields(2, new_seq_db_path);

    _state = _parent->state.rename_database(old_name, new_name);
}

const std::string& Database::getName() const {
    return _state.name;
}

void Database::setDescription(const std::string& new_description) {
    _state = _parent->state.set_database_description(_state.name, new_description);
}

const std::string& Database::getDescription() const {
    return _state.description;
}

const std::string& Database::getType() const {
    return _state.database_type;
}
