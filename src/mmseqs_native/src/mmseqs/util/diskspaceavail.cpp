#include <mmseqs/commons/command.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

int diskspaceavail(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  size_t diskLimit = FileUtil::getFreeSpace(out, FileUtil::dirName(par.db1).c_str());
  out->info("Disk limit: {} bytes", diskLimit);
  return EXIT_SUCCESS;
}
