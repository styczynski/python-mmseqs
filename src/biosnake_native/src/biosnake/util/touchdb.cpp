#include <biosnake/commons/memoryMapped.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/prefiltering/prefilteringIndexReader.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

int touchdb(biosnake_output* out, Parameters& par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  std::string db = par.db1;

  std::string indexDB = PrefilteringIndexReader::searchForIndex(out, db);
  if (indexDB.empty() == false) {
    db = indexDB;
  }

  MemoryMapped map(out, db, MemoryMapped::WholeFile,
                   MemoryMapped::CacheHint::SequentialScan);
  Util::touchMemory(out, reinterpret_cast<const char*>(map.getData()),
                    map.mappedSize());

  return EXIT_SUCCESS;
}
