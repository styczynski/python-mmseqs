//
// Created by mad on 11/29/16.
//

#ifndef BIOSNAKE_LIBRARYREADER_H
#define BIOSNAKE_LIBRARYREADER_H

#include <biosnake/output.h>
#include <sstream>
#include <vector>

class LibraryReader {
 public:
  bool StreamStartsWith(std::stringstream& in, const char* id);
  int ReadInt(biosnake_output* out, const char* line, const char* label, const char* errmsg);
  double ReadDouble(biosnake_output* out, const char* line, const char* label, const char* errmsg);
  std::string ReadString(biosnake_output* out, const char* line, const char* label,
                         const char* errmsg);
  bool ReadBool(biosnake_output* out, const char* line, const char* label, const char* errmsg);
  const char* strscn(const char* str);
  static std::vector<std::string> tokenize(const char* str, char sep);
  std::string getline(std::stringstream& in);
};

#endif  // BIOSNAKE_LIBRARYREADER_H
