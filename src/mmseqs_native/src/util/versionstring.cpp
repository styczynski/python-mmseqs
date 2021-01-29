#include "command.h"
#include "debug.h"
#include "parameters.h"
#include "util.h"
#include "output.h"

int versionstring(mmseqs_output* out, Parameters& par) {
  Debug(Debug::INFO) << version << "\n";
  EXIT(EXIT_SUCCESS);
}
