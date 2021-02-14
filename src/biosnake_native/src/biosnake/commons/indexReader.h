#ifndef BIOSNAKE_INDEXREADER_H
#define BIOSNAKE_INDEXREADER_H

#include <biosnake/commons/dBReader.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/prefiltering/prefilteringIndexReader.h>

class IndexReader {
 public:
  const static int PRELOAD_NO = 0;
  const static int PRELOAD_DATA = 1;
  const static int PRELOAD_INDEX = 2;
  IndexReader(biosnake_output* output, const std::string &dataName, int threads,
              int databaseType = SEQUENCES | HEADERS, int preloadMode = false,
              int dataMode = (DBReader<unsigned int>::USE_INDEX |
                              DBReader<unsigned int>::USE_DATA))
      : out(output), sequenceReader(NULL), index(NULL) {
    int targetDbtype = FileUtil::parseDbType(out, dataName.c_str());
    if (Parameters::isEqualDbtype(targetDbtype, Parameters::DBTYPE_INDEX_DB)) {
      index = new DBReader<unsigned int>(
          out, dataName.c_str(), (dataName + ".index").c_str(), 1,
          DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
      index->open(DBReader<unsigned int>::NOSORT);
      if (PrefilteringIndexReader::checkIfIndexFile(index)) {
        PrefilteringIndexReader::printSummary(out, index);
        PrefilteringIndexData data =
            PrefilteringIndexReader::getMetadata(index);
        seqType = data.seqType;
        bool touchIndex = preloadMode & PRELOAD_INDEX;
        bool touchData = preloadMode & PRELOAD_DATA;
        if (databaseType & SRC_SEQUENCES) {
          sequenceReader = PrefilteringIndexReader::openNewReader(
              out,
              index, PrefilteringIndexReader::DBR2DATA,
              PrefilteringIndexReader::DBR2INDEX,
              dataMode & DBReader<unsigned int>::USE_DATA, threads, touchIndex,
              touchData);

        } else if (databaseType & SEQUENCES) {
          sequenceReader = PrefilteringIndexReader::openNewReader(
              out,
              index, PrefilteringIndexReader::DBR1DATA,
              PrefilteringIndexReader::DBR1INDEX,
              dataMode & DBReader<unsigned int>::USE_DATA, threads, touchIndex,
              touchData);
        } else if (databaseType & SRC_HEADERS) {
          sequenceReader = PrefilteringIndexReader::openNewHeaderReader(
              out, index, PrefilteringIndexReader::HDR2DATA,
              PrefilteringIndexReader::HDR2INDEX, threads, touchIndex,
              touchData);
        } else if (databaseType & HEADERS) {
          sequenceReader = PrefilteringIndexReader::openNewHeaderReader(
              out,
              index, PrefilteringIndexReader::HDR1DATA,
              PrefilteringIndexReader::HDR1INDEX, threads, touchIndex,
              touchData);
        }
        if (sequenceReader == NULL) {
          out->info("Index does not contain plain sequences. Using normal database instead.");
        }
        seqType = Parameters::DBTYPE_INDEX_DB;
      } else {
        out->warn("Outdated index version. Please recompute with 'createindex'");
        index->close();
        delete index;
        index = NULL;
      }
    }

    if (sequenceReader == NULL) {
      if (databaseType & (HEADERS | SRC_HEADERS)) {
        sequenceReader = new DBReader<unsigned int>(
            out, (dataName + "_h").c_str(), ((dataName + "_h") + ".index").c_str(),
            threads, dataMode);
      } else {
        sequenceReader = new DBReader<unsigned int>(
            out, dataName.c_str(), (dataName + ".index").c_str(), threads, dataMode);
      }
      sequenceReader->open(DBReader<unsigned int>::NOSORT);
      bool touchData = preloadMode & PRELOAD_DATA;
      if (touchData) {
        sequenceReader->readMmapedDataInMemory();
      }
      seqType = sequenceReader->getDbtype();
    }
  }

  static const int SEQUENCES = 1;
  static const int HEADERS = 2;
  static const int SRC_HEADERS = 4;
  static const int SRC_SEQUENCES = 8;

  int getDbtype() const { return seqType; }

  ~IndexReader() {
    if (sequenceReader != NULL) {
      sequenceReader->close();
      delete sequenceReader;
    }

    if (index != NULL) {
      index->close();
      delete index;
    }
  }

  DBReader<unsigned int> *sequenceReader;

 private:
  biosnake_output* out;
  DBReader<unsigned int> *index;
  int seqType;
};

#endif
