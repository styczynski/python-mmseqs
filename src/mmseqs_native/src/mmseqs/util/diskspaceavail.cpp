#include <mmseqs/commons/command.h>
#include <mmseqs/commons/debug.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

int diskspaceavail(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  size_t diskLimit = FileUtil::getFreeSpace(FileUtil::dirName(par.db1).c_str());
  Debug(Debug::INFO) << diskLimit << "\n";  // in bytes
  EXIT(EXIT_SUCCESS);
}
