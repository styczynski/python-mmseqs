//
// Created by mad on 11/29/16.
//

#ifndef MMSEQS_LIBRARYREADER_H
#define MMSEQS_LIBRARYREADER_H

#include <mmseqs/output.h>
#include <sstream>
#include <vector>

class LibraryReader {
 public:
  bool StreamStartsWith(std::stringstream& in, const char* id);
  int ReadInt(mmseqs_output* out, const char* line, const char* label, const char* errmsg);
  double ReadDouble(mmseqs_output* out, const char* line, const char* label, const char* errmsg);
  std::string ReadString(mmseqs_output* out, const char* line, const char* label,
                         const char* errmsg);
  bool ReadBool(mmseqs_output* out, const char* line, const char* label, const char* errmsg);
  const char* strscn(const char* str);
  static std::vector<std::string> tokenize(const char* str, char sep);
  std::string getline(std::stringstream& in);
};

#endif  // MMSEQS_LIBRARYREADER_H
