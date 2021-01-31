#include <mmseqs/commons/dBReader.h>
#include <mmseqs/commons/dBWriter.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/output.h>

#ifdef OPENMP
#include <omp.h>
#endif

int doCompression(mmseqs_output* out, Parameters& par, bool shouldCompress) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  DBReader<unsigned int> reader(
      par.db1.c_str(), par.db1Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  reader.open(DBReader<unsigned int>::NOSORT);
  if (shouldCompress == true && reader.isCompressed() == true) {
    out->info("Database is already compressed.");
    return EXIT_SUCCESS;
  }
  if (shouldCompress == false && reader.isCompressed() == false) {
    out->info("Database is already decompressed.");
    return EXIT_SUCCESS;
  }

  int dbtype = reader.getDbtype();
  dbtype = shouldCompress ? dbtype | (1 << 31) : dbtype & ~(1 << 31);
  DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                  shouldCompress, dbtype);
  writer.open();
  Debug::Progress progress(reader.getSize());

#pragma omp parallel
  {
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = (unsigned int)omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < reader.getSize(); ++i) {
      progress.updateProgress();
      writer.writeData(
          reader.getData(i, thread_idx),
          std::max(static_cast<unsigned int>(reader.getEntryLen(i)), 1u) - 1u,
          reader.getDbKey(i), thread_idx);
    }
  }
  writer.close();
  reader.close();

  return EXIT_SUCCESS;
}

int compress(mmseqs_output* out, Parameters& par) {
  return doCompression(out, par, true);
}

int decompress(mmseqs_output* out, Parameters& par) {
  return doCompression(out, par, false);
}
