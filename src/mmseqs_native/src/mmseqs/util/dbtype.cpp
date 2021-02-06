#include <mmseqs/commons/dBReader.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

int dbtype(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, false, 0, 0);
  out->info("Database type: {}", Parameters::getDbTypeName(
      FileUtil::parseDbType(out, par.db1.c_str())));
  return EXIT_SUCCESS;
}
