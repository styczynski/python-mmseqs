#ifndef MMSEQS_COMMAND_H
#define MMSEQS_COMMAND_H

#include <vector>
#include <mmseqs/output.h>

struct MMseqsParameter;
struct DbType;

struct DbValidator {
  static std::vector<int> sequenceDb;
  static std::vector<int> nuclDb;
  static std::vector<int> aaDb;
  static std::vector<int> prefAlnResDb;
  static std::vector<int> taxSequenceDb;
  static std::vector<int> nuclAaDb;
  static std::vector<int> alignmentDb;
  static std::vector<int> prefilterDb;
  static std::vector<int> clusterDb;
  static std::vector<int> resultDb;
  static std::vector<int> ca3mDb;
  static std::vector<int> msaDb;
  static std::vector<int> genericDb;
  static std::vector<int> profileDb;
  static std::vector<int> csDb;
  static std::vector<int> indexDb;
  static std::vector<int> allDb;
  static std::vector<int> allDbAndFlat;
  static std::vector<int> taxResult;
  static std::vector<int> directory;
  static std::vector<int> flatfile;
  static std::vector<int> flatfileAndStdin;
  static std::vector<int> flatfileStdinAndGeneric;
  static std::vector<int> empty;
};

struct DbType {
  static const int ACCESS_MODE_INPUT = 1;
  static const int ACCESS_MODE_OUTPUT = 2;
  static const int NEED_DATA = 0;
  static const int NEED_HEADER = 1;
  static const int NEED_LOOKUP = 2;
  static const int NEED_TAXONOMY = 4;
  static const int VARIADIC = 8;
  static const int ZERO_OR_ALL = 16;

  const char *usageText;
  int accessMode;
  int specialType;
  std::vector<int> *validator;
};

#endif
