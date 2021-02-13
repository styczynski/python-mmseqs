#include <biosnake/commons/parameters.h>
#include <biosnake/output.h>

struct biosnake_call_args {
  biosnake_call_args() { cli_args = std::vector<std::string>(); }
  std::vector<std::string> cli_args;
};

biosnake_output call_biosnake(std::string command_name, Parameters args);
void subcall_biosnake(biosnake_output* out, std::string command_name,
                    Parameters args);
