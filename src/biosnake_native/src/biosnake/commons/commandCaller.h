#ifndef BIOSNAKE_COMMANDCALLER_H
#define BIOSNAKE_COMMANDCALLER_H

#include <biosnake/output.h>
#include <cstddef>
#include <string>
#include <vector>

class CommandCaller {
 public:
  CommandCaller(biosnake_output* output);

  void addVariable(const char* key, const char* value);
  void addVar(std::string key, std::string value);

  int callProgram(const char* program, size_t argc, const char** argv);

  static unsigned int getCallDepth(biosnake_output* out);

  // Does not return on success
  void execProgram(const char* program, const std::vector<std::string>& argv);
  int callProgram(const char* program, const std::vector<std::string>& argv);

 private:
  biosnake_output* out;
};

#endif  // BIOSNAKE_COMMANDCALLER_H
