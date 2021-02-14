#ifndef RESULT2PROFILE_H
#define RESULT2PROFILE_H

#include <biosnake/commons/dBReader.h>
#include <biosnake/commons/dBWriter.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/output.h>

#include <unordered_map>

class StatsComputer {
 public:
  StatsComputer(biosnake_output* output, const Parameters &par);
  ~StatsComputer();

  int run();

 private:
  biosnake_output* out;
  int stat;

  std::string queryDb;
  std::string queryDbIndex;

  std::string targetDb;
  std::string targetDbIndex;

  const bool tsvOut;

  DBReader<unsigned int> *resultReader;
  DBWriter *statWriter;

  int threads;

  template <typename T>
  struct PerSequence {
    typedef T (*type)(const char *);
  };

  template <typename T>
  int sequenceWise(biosnake_output* out, typename PerSequence<T>::type call,
                   bool onlyResultDb = false);

  int countNumberOfLines();
  int meanValue();
  int sumValue();
};

#endif
