#ifndef MMSEQS_COMMANDCALLER_H
#define MMSEQS_COMMANDCALLER_H

#include <cstddef>
#include <string>
#include <vector>

class CommandCaller {
 public:
  CommandCaller(mmseqs_output* output);

  void addVariable(const char* key, const char* value);
  void addVar(std::string key, std::string value);

  int callProgram(const char* program, size_t argc, const char** argv);

  static unsigned int getCallDepth();

  // Does not return on success
  void execProgram(const char* program, const std::vector<std::string>& argv);
  int callProgram(const char* program, const std::vector<std::string>& argv);

 private:
  mmseqs_output* out;
};

#endif  // MMSEQS_COMMANDCALLER_H
