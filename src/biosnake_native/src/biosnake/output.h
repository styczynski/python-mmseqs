#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <map>
#include <string>

struct biosnake_output;

template <typename T>
struct assert_false : std::false_type {};

template <typename T>
std::string SSTR(T) {
  static_assert(assert_false<T>::value, "Not implemented for requested type");
  return "";
}

template <>
std::string SSTR(const char *);
template <>
std::string SSTR(char *);
template <>
std::string SSTR(bool);
template <>
std::string SSTR(const char[]);
template <>
std::string SSTR(const std::string &);
template <>
std::string SSTR(std::string);
template <>
std::string SSTR(char);
template <>
std::string SSTR(short);
template <>
std::string SSTR(unsigned short);
template <>
std::string SSTR(int);
template <>
std::string SSTR(unsigned int);
template <>
std::string SSTR(long);
template <>
std::string SSTR(unsigned long);
template <>
std::string SSTR(long long);
template <>
std::string SSTR(unsigned long long);
template <>
std::string SSTR(double);
template <>
std::string SSTR(float);

#include <biosnake/commons/util.h>
#include <biosnake/commons/log.h>

struct biosnake_blast_tab_record {
  std::string query_sequence_id;
  std::string target_sequence_id;
  std::string query_sequence_content;
  std::string target_sequence_content;
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
    return std::string("BiosnakeSearchRecord(") + "query_sequence_id='" +
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

struct biosnake_output {
  biosnake_output(std::string file_log_output = ""): logger(file_log_output) {
    vars_str = std::map<std::string, std::string>();
    blast_tab_records = std::vector<std::vector<biosnake_blast_tab_record>>();
  }

  std::map<std::string, std::string> vars_str;
  std::vector<std::vector<biosnake_blast_tab_record>> blast_tab_records;
  Log logger;

  template<typename FormatString, typename... Args>
  void error(const FormatString &fmt, Args&&...args) {
      logger.error(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void warn(const FormatString &fmt, Args&&...args) {
      logger.warn(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void info(const FormatString &fmt, Args&&...args) {
      logger.info(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void debug(const FormatString &fmt, Args&&...args) {
      logger.debug(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void failure(const FormatString &fmt, Args&&...args) {
      logger.failure(fmt, std::forward<Args>(args)...);
  }

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
