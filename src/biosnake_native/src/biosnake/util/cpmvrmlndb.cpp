#include <biosnake/commons/dBReader.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/output.h>

int rmdb(biosnake_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int>::removeDb(out, par.db1);
  return EXIT_SUCCESS;
}

int mvdb(biosnake_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int>::moveDb(out, par.db1.c_str(), par.db2.c_str());
  return EXIT_SUCCESS;
}

int cpdb(biosnake_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int>::copyDb(out, par.db1.c_str(), par.db2.c_str());
  return EXIT_SUCCESS;
}

int lndb(biosnake_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int>::softlinkDb(out, par.db1.c_str(), par.db2.c_str());
  return EXIT_SUCCESS;
}
