#include <biosnake/commons/command.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

int diskspaceavail(biosnake_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  size_t diskLimit = FileUtil::getFreeSpace(out, FileUtil::dirName(out, par.db1).c_str());
  out->info("Disk limit: {} bytes", diskLimit);
  return EXIT_SUCCESS;
}
