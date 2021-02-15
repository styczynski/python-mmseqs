#ifndef API_DEFAULTS_H
#define API_DEFAULTS_H

#include <string>
#include <vector>
#include <map>

static std::vector<std::string> PARAM_SEARCH_COL_NAMES_DEFAULT = {
  "query_sequence_id",
  "target_sequence_id",
  "sequence_identity",
  "alignment_length",
  "number_of_mismatches",
  "number_of_gap_openings",
  "domain_start_index_query",
  "domain_end_index_query",
  "domain_start_index_target",
  "domain_end_index_target",
  "e_value",
  "bit_score"
};

static std::map<std::string, std::string> PARAM_SEARCH_COL_NAMES_MAPPING = {
  {"query_sequence_id", "query"},
  {"target_sequence_id", "target"},
  {"query_sequence_content", "qseq"},
  {"target_sequence_content", "tseq"},
  {"query_sequence_aligned", "qaln"},
  {"target_sequence_aligned", "taln"},
  {"sequence_identity", "fident"},
  {"alignment_length", "alnlen"},
  {"number_of_mismatches", "mismatch"},
  {"number_of_gap_openings", "gapopen"},
  {"domain_start_index_query", "qstart"},
  {"domain_end_index_query", "qend"},
  {"domain_start_index_target", "tstart"},
  {"domain_end_index_target", "tend"},
  {"e_value", "evalue"},
  {"bit_score", "bits"}
};

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

static const std::string DEFAULT_SEARCH_TYPE = "auto";
static const float DEFAULT_SENSITIVITY = 6;
static const int DEFAULT_MAX_SEQUENCE_LENGTH = 30000;
static const int DEFAULT_MAX_ORF_LENGTH = 32734;
static const int DEFAULT_MIN_ORF_LENGTH = 30;
static const int DEFAULT_ORF_START_MODE = 1;
static const int DEFAULT_MAX_RESULTS_COUNT_PER_QUERY = 300;
static const int DEFAULT_SEARCH_STEPS = 1;
static const float DEFAULT_START_SENSITIVITY = 4.0;

#endif