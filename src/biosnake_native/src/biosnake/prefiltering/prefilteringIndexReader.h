#ifndef PREFILTERINGINDEXREADER_H
#define PREFILTERINGINDEXREADER_H

#include <string>
#include <biosnake/commons/baseMatrix.h>
#include <biosnake/commons/dBReader.h>
#include <biosnake/prefiltering/indexTable.h>
#include <biosnake/output.h>

struct PrefilteringIndexData {
  int maxSeqLength;
  int kmerSize;
  int compBiasCorr;
  int alphabetSize;
  int mask;
  int spacedKmer;
  int kmerThr;
  int seqType;
  int srcSeqType;
  int headers1;
  int headers2;
  int splits;
};

class PrefilteringIndexReader {
 public:
  static const char *CURRENT_VERSION;
  static unsigned int VERSION;
  static unsigned int ENTRIES;
  static unsigned int ENTRIESOFFSETS;
  static unsigned int ENTRIESGRIDSIZE;
  static unsigned int SEQINDEXDATA;
  static unsigned int SEQINDEXDATASIZE;
  static unsigned int SEQINDEXSEQOFFSET;
  static unsigned int ENTRIESNUM;
  static unsigned int SEQCOUNT;
  static unsigned int META;
  static unsigned int SCOREMATRIXNAME;
  static unsigned int SEEDSCOREMATRIXNAME;
  static unsigned int SCOREMATRIX2MER;
  static unsigned int SCOREMATRIX3MER;
  static unsigned int DBR1INDEX;
  static unsigned int DBR1DATA;
  static unsigned int DBR2INDEX;
  static unsigned int DBR2DATA;
  static unsigned int HDR1INDEX;
  static unsigned int HDR1DATA;
  static unsigned int HDR2INDEX;
  static unsigned int HDR2DATA;
  static unsigned int GENERATOR;
  static unsigned int SPACEDPATTERN;

  static bool checkIfIndexFile(DBReader<unsigned int> *reader);
  static std::string indexName(const std::string &outDB);

  static void createIndexFile(
      biosnake_output *out, const std::string &outDb,
      DBReader<unsigned int> *dbr1, DBReader<unsigned int> *dbr2,
      DBReader<unsigned int> *hdbr1, DBReader<unsigned int> *hdbr2,
      BaseMatrix *seedSubMat, int maxSeqLen, bool spacedKmer,
      const std::string &spacedKmerPattern, bool compBiasCorrection,
      int alphabetSize, int kmerSize, int maskMode, int maskLowerCase,
      int kmerThr, int splits);

  static DBReader<unsigned int> *openNewHeaderReader(
      biosnake_output* out, DBReader<unsigned int> *dbr, unsigned int dataIdx, unsigned int indexIdx,
      int threads, bool touchIndex, bool touchData);

  static DBReader<unsigned int> *openNewReader(biosnake_output* out, DBReader<unsigned int> *dbr,
                                               unsigned int dataIdx,
                                               unsigned int indexIdx,
                                               bool includeData, int threads,
                                               bool touchIndex, bool touchData);

  static SequenceLookup *getSequenceLookup(biosnake_output* out, unsigned int split,
                                           DBReader<unsigned int> *dbr,
                                           int preloadMode);

  static IndexTable *getIndexTable(biosnake_output* out, unsigned int split,
                                   DBReader<unsigned int> *dbr,
                                   int preloadMode);

  static void printSummary(biosnake_output* out, DBReader<unsigned int> *dbr);

  static PrefilteringIndexData getMetadata(DBReader<unsigned int> *dbr);

  static std::string getSubstitutionMatrixName(DBReader<unsigned int> *dbr);

  static std::string getSubstitutionMatrix(DBReader<unsigned int> *dbr);

  static std::string getSpacedPattern(DBReader<unsigned int> *dbr);

  static ScoreMatrix get2MerScoreMatrix(DBReader<unsigned int> *dbr,
                                        int preloadMode);

  static ScoreMatrix get3MerScoreMatrix(DBReader<unsigned int> *dbr,
                                        int preloadMode);

  static std::string searchForIndex(biosnake_output* out, const std::string &pathToDB);

  static std::string dbPathWithoutIndex(std::string &dbname);

 private:
  static void printMeta(biosnake_output* out, int *meta);
};

#endif
