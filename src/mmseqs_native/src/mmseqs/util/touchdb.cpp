#include <mmseqs/commons/memoryMapped.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/prefiltering/prefilteringIndexReader.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

int touchdb(mmseqs_output* out, Parameters& par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  std::string db = par.db1;

  std::string indexDB = PrefilteringIndexReader::searchForIndex(db);
  if (indexDB.empty() == false) {
    db = indexDB;
  }

  MemoryMapped map(db, MemoryMapped::WholeFile,
                   MemoryMapped::CacheHint::SequentialScan);
  Util::touchMemory(reinterpret_cast<const char*>(map.getData()),
                    map.mappedSize());

  return EXIT_SUCCESS;
}