#include <biosnake/commons/application.h>
#include <biosnake/commons/command.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/timer.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#include <iomanip>

Command *getCommandByName(biosnake_output *out, const char *s) {
  for (size_t i = 0; i < commands.size(); i++) {
    Command &p = commands[i];
    if (!strcmp(s, p.cmd)) return &p;
  }
  out->failure("Cannot find a command called '{}'. Check if it has the correct spelling.", s);
  return NULL;
}

int runCommand(biosnake_output *out, Command *c, Parameters &par) {
  Timer timer;
  int status = c->commandFunction(out, par);
  out->info("Time for processing: {}", timer.lap());
  return status;
}

void subcall_biosnake(biosnake_output *out, std::string command_name,
                    Parameters args) {
  Parameters par(args);
  Command *c = getCommandByName(out, command_name.c_str());
  runCommand(out, c, par);
}

biosnake_output call_biosnake(std::string command_name, Parameters args) {
  biosnake_output out(args.logFilePath);
  subcall_biosnake(&out, command_name, args);
  return out;
}
