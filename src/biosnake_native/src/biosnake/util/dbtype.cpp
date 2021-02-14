#include <biosnake/commons/dBReader.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

int dbtype(biosnake_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, false, 0, 0);
  out->info("Database type: {}", Parameters::getDbTypeName(
      FileUtil::parseDbType(out, par.db1.c_str())));
  return EXIT_SUCCESS;
}
