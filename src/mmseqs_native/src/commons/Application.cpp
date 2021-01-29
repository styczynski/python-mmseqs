#include "Application.h"
#include "Command.h"
#include "Debug.h"
#include "Parameters.h"
#include "Timer.h"
#include "Util.h"
#include "output.h"

#include <iomanip>

Command *getCommandByName(const char *s) {
  for (size_t i = 0; i < commands.size(); i++) {
    Command &p = commands[i];
    if (!strcmp(s, p.cmd)) return &p;
  }
  return NULL;
}

int runCommand(mmseqs_output *out, Command *c, Parameters &par) {
  Timer timer;
  int status = c->commandFunction(out, par);
  Debug(Debug::INFO) << "Time for processing: " << timer.lap() << "\n";
  return status;
}

void subcall_mmseqs(mmseqs_output *out, std::string command_name,
                    Parameters args) {
  Parameters par(args);
  Command *c = getCommandByName(command_name.c_str());
  runCommand(out, c, par);
}

mmseqs_output call_mmseqs(std::string command_name, Parameters args) {
  mmseqs_output out;
  subcall_mmseqs(&out, command_name, args);
  return out;
}