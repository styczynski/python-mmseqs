#include <mmseqs/commons/parameters.h>
#include <mmseqs/output.h>

struct mmseqs_call_args {
  mmseqs_call_args() { cli_args = std::vector<std::string>(); }
  std::vector<std::string> cli_args;
};

mmseqs_output call_mmseqs(std::string command_name, Parameters args);
void subcall_mmseqs(mmseqs_output* out, std::string command_name,
                    Parameters args);
