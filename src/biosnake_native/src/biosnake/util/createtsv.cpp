#include <biosnake/commons/dBReader.h>
#include <biosnake/commons/dBWriter.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/indexReader.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#ifdef OPENMP
#include <omp.h>

#endif

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t)-1)
#endif

int createtsv(biosnake_output *out, Parameters &par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true,
  //    Parameters::PARSE_VARIADIC, 0);

  bool queryNucs = Parameters::isEqualDbtype(
      FileUtil::parseDbType(out, par.db1.c_str()), Parameters::DBTYPE_NUCLEOTIDES);
  bool targetNucs = Parameters::isEqualDbtype(
      FileUtil::parseDbType(out, par.db2.c_str()), Parameters::DBTYPE_NUCLEOTIDES);
  const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
  int queryHeaderType =
      (queryNucs) ? IndexReader::SRC_HEADERS : IndexReader::HEADERS;
  queryHeaderType = (par.idxSeqSrc == 0)
                        ? queryHeaderType
                        : (par.idxSeqSrc == 1) ? IndexReader::HEADERS
                                               : IndexReader::SRC_HEADERS;
  IndexReader qDbrHeader(
      out, par.db1, par.threads, queryHeaderType,
      (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
  IndexReader *tDbrHeader = NULL;
  DBReader<unsigned int> *queryDB = qDbrHeader.sequenceReader;
  DBReader<unsigned int> *targetDB = NULL;
  bool sameDB = (par.db2.compare(par.db1) == 0);
  const bool hasTargetDB = par.filenames.size() > 3;
  DBReader<unsigned int>::Index *qHeaderIndex =
      qDbrHeader.sequenceReader->getIndex();
  DBReader<unsigned int>::Index *tHeaderIndex = NULL;

  if (hasTargetDB) {
    if (sameDB) {
      tDbrHeader = &qDbrHeader;
      tHeaderIndex = qHeaderIndex;
      targetDB = queryDB;
    } else {
      int targetHeaderType =
          (targetNucs) ? IndexReader::SRC_HEADERS : IndexReader::HEADERS;
      targetHeaderType = (par.idxSeqSrc == 0)
                             ? targetHeaderType
                             : (par.idxSeqSrc == 1) ? IndexReader::HEADERS
                                                    : IndexReader::SRC_HEADERS;

      tDbrHeader =
          new IndexReader(out, par.db2, par.threads, targetHeaderType, touch);
      tHeaderIndex = tDbrHeader->sequenceReader->getIndex();
      targetDB = tDbrHeader->sequenceReader;
    }
  }

  DBReader<unsigned int> *reader;
  if (hasTargetDB) {
    reader = new DBReader<unsigned int>(
        out, par.db3.c_str(), par.db3Index.c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  } else {
    reader = new DBReader<unsigned int>(
        out, par.db2.c_str(), par.db2Index.c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  }
  reader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

  const std::string &dataFile = hasTargetDB ? par.db4 : par.db3;
  const std::string &indexFile = hasTargetDB ? par.db4Index : par.db3Index;
  const bool shouldCompress = par.dbOut == true && par.compressed == true;
  const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB
                                       : Parameters::DBTYPE_OMIT_FILE;
  DBWriter writer(out, dataFile.c_str(), indexFile.c_str(), par.threads,
                  shouldCompress, dbType);
  writer.open();

  const size_t targetColumn =
      (par.targetTsvColumn == 0) ? SIZE_T_MAX : par.targetTsvColumn - 1;
#pragma omp parallel
  {
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = (unsigned int)omp_get_thread_num();
#endif

    const char *columnPointer[255];
    char *dbKey = new char[par.maxSeqLen + 1];

    std::string outputBuffer;
    outputBuffer.reserve(10 * 1024);

#pragma omp for schedule(dynamic, 1000)
    for (size_t i = 0; i < reader->getSize(); ++i) {
      unsigned int queryKey = reader->getDbKey(i);
      size_t queryIndex = queryDB->getId(queryKey);

      char *headerData = queryDB->getData(queryIndex, thread_idx);
      if (headerData == NULL) {
        out->warn("Invalid header entry in query {}", queryKey);
        continue;
      }

      std::string queryHeader;
      if (par.fullHeader) {
        queryHeader = "\"";
        queryHeader.append(headerData, qHeaderIndex[queryIndex].length - 2);
        queryHeader.append("\"");
      } else {
        queryHeader = Util::parseFastaHeader(headerData);
      }

      size_t entryIndex = 0;

      char *data = reader->getData(i, thread_idx);
      while (*data != '\0') {
        if (targetColumn != SIZE_T_MAX) {
          size_t foundElements = Util::getWordsOfLine(data, columnPointer, 255);
          if (foundElements < targetColumn) {
            out->warn("Not enough columns!");
            continue;
          }
          Util::parseKey(columnPointer[targetColumn], dbKey);
        }
        std::string targetAccession;
        if (targetColumn == SIZE_T_MAX) {
          targetAccession = "";
        } else if (hasTargetDB) {
          unsigned int targetKey = (unsigned int)strtoul(dbKey, NULL, 10);
          size_t targetIndex = targetDB->getId(targetKey);
          char *targetData = targetDB->getData(targetIndex, thread_idx);
          if (targetData == NULL) {
            out->warn("Invalid header entry in query {} and target {}", queryKey, targetKey);
            continue;
          }
          if (par.fullHeader) {
            targetAccession = "\"";
            targetAccession.append(targetData,
                                   tHeaderIndex[targetIndex].length - 2);
            targetAccession.append("\"");
          } else {
            targetAccession = Util::parseFastaHeader(targetData);
          }
        } else {
          targetAccession = dbKey;
        }

        if (par.firstSeqRepr && !entryIndex) {
          queryHeader = targetAccession;
        }

        outputBuffer.append(queryHeader);
        outputBuffer.append("\t");
        outputBuffer.append(targetAccession);

        size_t offset = 0;
        if (targetColumn != 0) {
          outputBuffer.append("\t");
          offset = 0;
        } else {
          offset = strlen(dbKey);
        }

        char *nextLine = Util::skipLine(data);
        outputBuffer.append(data + offset, (nextLine - (data + offset)) - 1);
        outputBuffer.append("\n");
        data = nextLine;
        entryIndex++;
      }
      writer.writeData(outputBuffer.c_str(), outputBuffer.length(), queryKey,
                       thread_idx, par.dbOut);
      outputBuffer.clear();
    }
    delete[] dbKey;
  }
  writer.close(par.dbOut == false);

  if (par.dbOut == false) {
    if (hasTargetDB) {
      FileUtil::remove(out, par.db4Index.c_str());
    } else {
      FileUtil::remove(out, par.db3Index.c_str());
    }
  }

  reader->close();
  delete reader;
  if (hasTargetDB) {
    if (sameDB == false) {
      delete tDbrHeader;
    }
  }

  return EXIT_SUCCESS;
}
#undef SIZE_T_MAX
