
#include <biosnake/commons/application.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/output.h>
#include <biosnake/api/database.h>
#include <biosnake/api/rand.h>
#include <biosnake/api/utils.h>
#include <biosnake/commons/fileUtil.h>


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
    call_biosnake("convert2fasta", args);
    return output_path;
}

Database Database::copy(std::string db_copy_name) {
    _parent->prepare_to_execute_command();

    if (db_copy_name.size() == 0) {
        int no = 0;
        while (no<5000) {
            biosnake_output* out;
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
    call_biosnake("cpdb", args);

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
    call_biosnake("rmdb", args);

    args.setDBFields(1, seq_db_path_h);
    call_biosnake("rmdb", args);

    _parent->state.remove_database(_state.name);
    return _state.name;
}

SearchResults Database::search(
    std::vector<std::string> sequences,
    std::string search_type,
    std::vector<std::string> headers,
    float sensitivity,
    int max_sequence_length,
    int max_results_count_per_query,
    int max_orf_length,
    int min_orf_length,
    int search_steps,
    float start_sensitivity
) {
    auto path = write_temp_fasta(sequences, _parent->get_workdir_path());
    auto search_results = search_file(path, search_type,
        headers, sensitivity,
        max_sequence_length, max_results_count_per_query,
        max_orf_length, min_orf_length,
        search_steps, start_sensitivity);
    biosnake_output* out;
    FileUtil::remove(out, path.c_str());
    return search_results;
}

SearchResults Database::search_file(
    std::string search_input_fasta,
    std::string search_type,
    std::vector<std::string> headers,
    float sensitivity,
    int max_sequence_length,
    int max_results_count_per_query,
    int max_orf_length,
    int min_orf_length,
    int search_steps,
    float start_sensitivity
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
    args.maxResListLen=max_results_count_per_query;
    args.orfMinLength=min_orf_length;
    args.orfMaxLength=max_orf_length;
    args.startSens=start_sensitivity;
    args.sensSteps=search_steps;
    args.splitMemoryLimit=5000000000;
    args.outfmt="";
    for (auto& header_name: headers) {
        args.outfmt += PARAM_SEARCH_COL_NAMES_MAPPING[header_name] + ",";
    }
    args.outfmt = args.outfmt.substr(0, args.outfmt.size()-1);

    args.searchType=PARAM_DB_SEARCH_TYPE_MAPPING[search_type];
    args.alignmentMode=3;

    auto search_output = call_biosnake("easy-search", args);

    biosnake_output* out;
    if (FileUtil::fileExists(out, results_path.c_str())) {
        FileUtil::remove(out, results_path.c_str());
    }
    if (FileUtil::fileExists(out, tmp_dir.c_str())) {
        FileUtil::remove(out, tmp_dir.c_str());
    }

    return SearchResults(search_output.blast_tab_records, headers);
}

biosnake_output Database::create_index(
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

    auto out = call_biosnake("createindex", args);
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

Database::DBIterator::DBIterator(std::shared_ptr<DBReader<unsigned int>> db_reader, std::shared_ptr<DBReader<unsigned int>> header_reader, Database::DBIterator::pointer ptr): _db_reader(db_reader), _header_reader(header_reader), m_ptr(ptr) {}

Database::DBIterator::reference Database::DBIterator::operator*() {
    unsigned int key = _db_reader->getDbKey(m_ptr);
    unsigned int headerKey = _header_reader->getId(key);
    const char* headerData = _header_reader->getData(headerKey, 0);
    const size_t headerLen = strlen(headerData);

    unsigned int bodyKey = _db_reader->getId(key);
    const char* bodyData = _db_reader->getData(bodyKey, 0);
    const size_t bodyLen = _db_reader->getEntryLen(bodyKey);

    return Record(bodyData, bodyLen-2, headerData, headerLen-1);
}
Database::DBIterator::pointer Database::DBIterator::operator->() { return m_ptr; }
Database::DBIterator& Database::DBIterator::operator++() { m_ptr++; return *this; }
Database::DBIterator Database::DBIterator::operator++(int) { Database::DBIterator tmp = *this; ++(*this); return tmp; }
bool operator== (const Database::DBIterator& a, const Database::DBIterator& b) { return a.m_ptr == b.m_ptr; };
bool operator!= (const Database::DBIterator& a, const Database::DBIterator& b) { return a.m_ptr != b.m_ptr; };
Database::DBIterator Database::begin() {
    _init_readers();
    return Database::DBIterator(db_reader, header_reader, 0);
}
Database::DBIterator Database::end() {
    _init_readers();
    return Database::DBIterator(db_reader, header_reader, db_reader->getSize());
}

void Database::_init_readers() {
    if (!db_reader) {
        std::string db_path = _parent->get_storage_path(_state.name);
        std::string db_index_path = _parent->get_storage_path(_state.name + ".index");
        db_reader = std::make_shared<DBReader<unsigned int>>(&_out, db_path.c_str(), db_index_path.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        db_reader->open(DBReader<unsigned int>::NOSORT);
    }
    if (!header_reader) {
        std::string header_path = _parent->get_storage_path(_state.name + "_h");
        std::string header_index_path = _parent->get_storage_path(_state.name + "_h.index");
        header_reader = std::make_shared<DBReader<unsigned int>>(&_out, header_path.c_str(), header_index_path.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        header_reader->open(DBReader<unsigned int>::NOSORT);
    }
}

const int Database::getSize() {
    _init_readers();
    return db_reader->getSize();
}

const std::tuple<py::array, py::array, py::array> Database::getColumnData() {
     _init_readers();

     const int count = db_reader->getSize();
     std::vector<std::string> headers;
     std::vector<std::string> sequences;
     std::vector<int> lengths;

     headers.reserve(count);
     sequences.reserve(count);
     lengths.reserve(count);

     for (size_t i=0; i<count; ++i) {
        unsigned int key = db_reader->getDbKey(i);
        unsigned int headerKey = header_reader->getId(key);
        const char* headerData = header_reader->getData(headerKey, 0);
        const size_t headerLen = strlen(headerData);

        unsigned int bodyKey = db_reader->getId(key);
        const char* bodyData = db_reader->getData(bodyKey, 0);
        const size_t bodyLen = db_reader->getEntryLen(bodyKey);

        headers.push_back(std::string(headerData, headerLen-1));
        sequences.push_back(std::string(bodyData, bodyLen-2));
        lengths.push_back(bodyLen-2);
     }
     return std::make_tuple(py::array(py::cast(headers)), py::array(py::cast(sequences)), py::array(py::cast(lengths)));
}