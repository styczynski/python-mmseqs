#include <mmseqs/commons/dBReader.h>
#include <mmseqs/commons/dBWriter.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

#ifdef OPENMP
#include <omp.h>
#endif

int addid(mmseqs_output* out, const std::string &db1, const std::string &db1Index,
          const std::string &db2, const std::string &db2Index,
          const bool tsvOut, const std::string &mappingFile,
          const std::string &userStrToAdd, const bool isPrefix,
          const int threads, const int compressed) {
  DBReader<unsigned int> reader(
      out, db1.c_str(), db1Index.c_str(), threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

  const bool shouldCompress = tsvOut == false && compressed == true;
  // TODO: does generic db make more sense than copying db type here?
  const int dbType =
      tsvOut == true ? Parameters::DBTYPE_OMIT_FILE : reader.getDbtype();
  DBWriter writer(out, db2.c_str(), db2Index.c_str(), threads, shouldCompress,
                  dbType);
  writer.open();
  const bool shouldWriteNullByte = !tsvOut;

  size_t entries = reader.getSize();
  Log::Progress progress(entries);
  bool doMapping = false;
  DBReader<unsigned int> *lookupReader = NULL;
  if (mappingFile.size() > 0) {
    lookupReader =
        new DBReader<unsigned int>(out, mappingFile.c_str(), mappingFile.c_str(), 1,
                                   DBReader<unsigned int>::USE_LOOKUP);
    doMapping = true;
  }

#pragma omp parallel
  {
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

#pragma omp for schedule(dynamic, 100)
    for (size_t i = 0; i < entries; ++i) {
      progress.updateProgress();

      unsigned int key = reader.getDbKey(i);
      std::istringstream data(reader.getData(i, thread_idx));
      std::ostringstream ss;

      std::string line;
      while (std::getline(data, line)) {
        std::string strToAdd = "";
        if (userStrToAdd != "") {
          strToAdd = userStrToAdd;
        } else if (doMapping) {
          size_t lookupId = lookupReader->getLookupIdByKey(key);
          if (lookupId == SIZE_MAX) {
            out->failure("Could not find key {} in lookup", key);
          }
          strToAdd = lookupReader->getLookupEntryName(lookupId);
        } else {
          strToAdd = SSTR(key);
        }

        if (isPrefix) {
          ss << strToAdd << "\t" << line << "\n";
        } else {
          ss << line << "\t" << strToAdd << "\n";
        }
      }

      std::string result = ss.str();
      writer.writeData(result.c_str(), result.length(), key, thread_idx,
                       shouldWriteNullByte);
    }
  }
  writer.close(tsvOut);
  if (tsvOut) {
    FileUtil::remove(out, writer.getIndexFileName());
  }
  reader.close();
  if (doMapping) {
    delete lookupReader;
  }
  return EXIT_SUCCESS;
}

int prefixid(mmseqs_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  return (addid(out, par.db1, par.db1Index, par.db2, par.db2Index, par.tsvOut,
                par.mappingFile, par.prefix, true, par.threads,
                par.compressed));
}

int suffixid(mmseqs_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  return (addid(out, par.db1, par.db1Index, par.db2, par.db2Index, par.tsvOut,
                par.mappingFile, par.prefix, false, par.threads,
                par.compressed));
}
