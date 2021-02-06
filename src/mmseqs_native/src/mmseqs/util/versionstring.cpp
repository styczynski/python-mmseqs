#include <mmseqs/commons/command.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

int versionstring(mmseqs_output* out, Parameters& par) {
  out->info("Version: {}", version);
  return EXIT_SUCCESS;
}
