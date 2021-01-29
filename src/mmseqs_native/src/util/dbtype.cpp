#include "dBReader.h"
#include "debug.h"
#include "fileUtil.h"
#include "parameters.h"
#include "util.h"
#include "output.h"

int dbtype(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, false, 0, 0);
  Debug(Debug::INFO) << Parameters::getDbTypeName(
      FileUtil::parseDbType(par.db1.c_str()));
  EXIT(EXIT_SUCCESS);
}
