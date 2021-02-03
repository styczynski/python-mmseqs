#include <mmseqs/clustering/clustering.h>
#include <mmseqs/clustering/clusteringAlgorithms.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/timer.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/commons/itoa.h>

Clustering::Clustering(mmseqs_output* output, const std::string &seqDB, const std::string &seqDBIndex,
                       const std::string &alnDB, const std::string &alnDBIndex,
                       const std::string &outDB, const std::string &outDBIndex,
                       unsigned int maxIteration, int similarityScoreType,
                       int threads, int compressed)
    : out(output),
      maxIteration(maxIteration),
      similarityScoreType(similarityScoreType),
      threads(threads),
      compressed(compressed),
      outDB(outDB),
      outDBIndex(outDBIndex) {
  seqDbr =
      new DBReader<unsigned int>(out, seqDB.c_str(), seqDBIndex.c_str(), threads,
                                 DBReader<unsigned int>::USE_INDEX);
  seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);

  alnDbr = new DBReader<unsigned int>(
      out, alnDB.c_str(), alnDBIndex.c_str(), threads,
      DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
  alnDbr->open(DBReader<unsigned int>::NOSORT);
}

Clustering::~Clustering() {
  delete seqDbr;
  delete alnDbr;
}

void Clustering::run(int mode) {
  Timer timer;
  DBWriter *dbw = new DBWriter(out, outDB.c_str(), outDBIndex.c_str(), 1, compressed,
                               Parameters::DBTYPE_CLUSTER_RES);
  dbw->open();

  std::pair<unsigned int, unsigned int> *ret;
  ClusteringAlgorithms *algorithm = new ClusteringAlgorithms(
      out, seqDbr, alnDbr, threads, similarityScoreType, maxIteration);

  if (mode == Parameters::GREEDY) {
    out->info("Clustering mode: Greedy");
    ret = algorithm->execute(4);
  } else if (mode == Parameters::GREEDY_MEM) {
    out->info("Clustering mode: Greedy Low Mem");
    ret = algorithm->execute(4);
  } else if (mode == Parameters::SET_COVER) {
    out->info("Clustering mode: Set Cover");
    ret = algorithm->execute(1);
  } else if (mode == Parameters::CONNECTED_COMPONENT) {
    out->info("Clustering mode: Connected Component");
    ret = algorithm->execute(3);
  } else {
    out->failure("Wrong clustering mode");
  }

  Timer timerWrite;

  size_t dbSize = alnDbr->getSize();
  size_t seqDbSize = seqDbr->getSize();
  size_t cluNum = (dbSize > 0) ? 1 : 0;
  for (size_t i = 1; i < dbSize; i++) {
    cluNum += (ret[i].first != ret[i - 1].first);
  }
  out->info("Total time: {}\n", timer.lap());
  out->info("\nSize of the sequence database: {}\n", seqDbSize
                    );
  out->info("Size of the alignment database: {}\n", dbSize);
  out->info("Number of clusters: {}\n", cluNum);

  out->info("Writing results");
  writeData(dbw, ret, dbSize);
  out->info("Took {}", timerWrite.lap());
  delete[] ret;
  delete algorithm;

  dbw->close(false, false);
  seqDbr->close();
  alnDbr->close();
  delete dbw;
}

void Clustering::writeData(DBWriter *dbw,
                           const std::pair<unsigned int, unsigned int> *ret,
                           size_t dbSize) {
  std::string resultStr;
  resultStr.reserve(1024 * 1024 * 1024);
  char buffer[32];
  unsigned int prevRepresentativeKey = UINT_MAX;
  for (size_t i = 0; i < dbSize; i++) {
    unsigned int currRepresentativeKey = ret[i].first;
    // write query key first
    if (prevRepresentativeKey != currRepresentativeKey) {
      if (prevRepresentativeKey != UINT_MAX) {  // skip first
        dbw->writeData(resultStr.c_str(), resultStr.length(),
                       prevRepresentativeKey);
      }
      resultStr.clear();
      char *outpos = Itoa::u32toa_sse2(currRepresentativeKey, buffer);
      resultStr.append(buffer, (outpos - buffer - 1));
      resultStr.push_back('\n');
    }
    unsigned int memberKey = ret[i].second;
    if (memberKey != currRepresentativeKey) {
      char *outpos = Itoa::u32toa_sse2(memberKey, buffer);
      resultStr.append(buffer, (outpos - buffer - 1));
      resultStr.push_back('\n');
    }

    prevRepresentativeKey = currRepresentativeKey;
  }
  if (prevRepresentativeKey != UINT_MAX) {
    dbw->writeData(resultStr.c_str(), resultStr.length(),
                   prevRepresentativeKey);
  }
}
