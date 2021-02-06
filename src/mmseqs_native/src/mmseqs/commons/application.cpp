#include <mmseqs/commons/application.h>
#include <mmseqs/commons/command.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/timer.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

#include <iomanip>

Command *getCommandByName(mmseqs_output *out, const char *s) {
  for (size_t i = 0; i < commands.size(); i++) {
    Command &p = commands[i];
    if (!strcmp(s, p.cmd)) return &p;
  }
  out->failure("Cannot find a command called '{}'. Check if it has the correct spelling.", s);
  return NULL;
}

int runCommand(mmseqs_output *out, Command *c, Parameters &par) {
  Timer timer;
  int status = c->commandFunction(out, par);
  out->info("Time for processing: {}", timer.lap());
  return status;
}

void subcall_mmseqs(mmseqs_output *out, std::string command_name,
                    Parameters args) {
  Parameters par(args);
  Command *c = getCommandByName(out, command_name.c_str());
  runCommand(out, c, par);
}

mmseqs_output call_mmseqs(std::string command_name, Parameters args) {
  mmseqs_output out(args.logFilePath);
  subcall_mmseqs(&out, command_name, args);
  return out;
}
