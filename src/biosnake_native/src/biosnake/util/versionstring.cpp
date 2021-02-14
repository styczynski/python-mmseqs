#include <biosnake/commons/command.h>
#include <biosnake/output.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

int versionstring(biosnake_output* out, Parameters& par) {
  out->info("Version: {}", version);
  return EXIT_SUCCESS;
}
