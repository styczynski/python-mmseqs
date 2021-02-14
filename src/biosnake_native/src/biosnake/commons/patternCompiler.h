#ifndef BIOSNAKE_PATTERNCOMPILER_H
#define BIOSNAKE_PATTERNCOMPILER_H
#include <regex.h>
#include <biosnake/output.h>
#include <biosnake/commons/util.h>

class PatternCompiler {
 public:
  PatternCompiler(biosnake_output* output, const char *pattern): out(output) {
    if (regcomp(&regex, pattern, REG_EXTENDED | REG_NEWLINE) != 0) {
      out->failure("Error in regex {}", pattern);
    }
  }

  ~PatternCompiler() { regfree(&regex); }

  bool isMatch(const char *target) {
    return regexec(&regex, target, 0, NULL, 0) == 0;
  }

  bool isMatch(const char *target, size_t nmatch, regmatch_t *pmatch) {
    return regexec(&regex, target, nmatch, pmatch, 0) == 0;
  }

 private:
  biosnake_output* out;
  regex_t regex;
};

#endif  // BIOSNAKE_PATTERNCOMPILER_H
