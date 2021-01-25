#include "Command.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "output.h"

int diskspaceavail(mmseqs_output* out, Parameters &par) {
//    Parameters &par = Parameters::getInstance();
    size_t diskLimit = FileUtil::getFreeSpace(FileUtil::dirName(par.db1).c_str());
    Debug(Debug::INFO) << diskLimit << "\n"; // in bytes
    EXIT(EXIT_SUCCESS);
}
