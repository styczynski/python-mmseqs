#include <mmseqs/commons/dBReader.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/output.h>

int rmdb(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int>::removeDb(par.db1);
  return EXIT_SUCCESS;
}

int mvdb(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int>::moveDb(par.db1.c_str(), par.db2.c_str());
  return EXIT_SUCCESS;
}

int cpdb(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int>::copyDb(par.db1.c_str(), par.db2.c_str());
  return EXIT_SUCCESS;
}

int lndb(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int>::softlinkDb(par.db1.c_str(), par.db2.c_str());
  return EXIT_SUCCESS;
}
