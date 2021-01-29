#include <mmseqs/commons/command.h>
#include <mmseqs/commons/debug.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

int versionstring(mmseqs_output* out, Parameters& par) {
  Debug(Debug::INFO) << version << "\n";
  EXIT(EXIT_SUCCESS);
}
