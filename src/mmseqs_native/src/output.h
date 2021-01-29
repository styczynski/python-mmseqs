#ifndef OUTPUT_H
#define OUTPUT_H

#include "util.h"

#include <iostream>
#include <map>
#include <string>

struct mmseqs_blast_tab_record {
  std::string query_sequence_id;
  std::string target_sequence_id;
  float sequence_identity;
  int alignment_length;
  int number_of_mismatches;
  int number_of_gap_openings;
  int domain_start_index_query;
  int domain_end_index_query;
  int domain_start_index_target;
  int domain_end_index_target;
  float e_value;
  int bit_score;

  std::string toString() const {
    return std::string("MMSeqsSearchRecord(") + "query_sequence_id='" +
           SSTR(query_sequence_id) + "', " + "target_sequence_id='" +
           SSTR(target_sequence_id) + "', " +
           "sequence_identity=" + SSTR(sequence_identity) + ", " +
           "alignment_length=" + SSTR(alignment_length) + ", " +
           "number_of_mismatches=" + SSTR(number_of_mismatches) + ", " +
           "domain_start_index_query=" + SSTR(domain_start_index_query) + ", " +
           "domain_end_index_query=" + SSTR(domain_end_index_query) + ", " +
           "domain_start_index_target=" + SSTR(domain_start_index_target) +
           ", " + "domain_end_index_target=" + SSTR(domain_end_index_target) +
           ", " + "e_value=" + SSTR(e_value) + ", " +
           "bit_score=" + SSTR(bit_score) + ")";
  }
};

struct mmseqs_output {
  mmseqs_output() {
    vars_str = std::map<std::string, std::string>();
    blast_tab_records = std::vector<std::vector<mmseqs_blast_tab_record>>();
  }
  std::map<std::string, std::string> vars_str;
  std::vector<std::vector<mmseqs_blast_tab_record>> blast_tab_records;

  void output_string(std::string name, std::string value) {
    vars_str[name] = value;
  }

  void print() {
    std::cout << "Output{\n";
    for (auto it = vars_str.cbegin(); it != vars_str.cend(); ++it) {
      std::cout << "  " << it->first << ": [" << it->second << "]\n";
    }
    std::cout << "}\n";
  }
};

#endif
