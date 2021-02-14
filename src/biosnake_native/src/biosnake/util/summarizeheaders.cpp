#include <biosnake/commons/dBReader.h>
#include <biosnake/commons/dBWriter.h>
#include <biosnake/output.h>
#include <biosnake/commons/headerSummarizer.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#ifdef OPENMP
#include <omp.h>
#endif

int summarizeheaders(biosnake_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  DBReader<unsigned int> queryReader(
      out, par.db1.c_str(), par.db1Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  queryReader.open(DBReader<unsigned int>::NOSORT);

  DBReader<unsigned int> targetReader(
      out, par.db2.c_str(), par.db2Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  targetReader.open(DBReader<unsigned int>::NOSORT);

  DBReader<unsigned int> reader(
      out, par.db3.c_str(), par.db3Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  reader.open(DBReader<unsigned int>::NOSORT);

  DBWriter writer(out, par.db4.c_str(), par.db4Index.c_str(), par.threads,
                  par.compressed, Parameters::DBTYPE_GENERIC_DB);
  writer.open();

  HeaderSummarizer *summarizer;
  if (par.headerType == Parameters::HEADER_TYPE_METACLUST) {
    summarizer = new MetaclustHeaderSummarizer;
  } else if (par.headerType == Parameters::HEADER_TYPE_UNICLUST) {
    summarizer = new UniprotHeaderSummarizer;
  } else {
    out->error("Header type is not supported");
    return EXIT_FAILURE;
  }
  Log::Progress progress(reader.getSize());

#pragma omp parallel
  {
    int thread_idx = 0;
#ifdef OPENMP
    thread_idx = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 100)
    for (size_t i = 0; i < reader.getSize(); ++i) {
      progress.updateProgress();

      unsigned int id = reader.getDbKey(i);
      char *data = reader.getData(i, thread_idx);

      std::vector<std::string> headers;

      std::istringstream inStream(data);
      std::string line;
      size_t entry = 0;
      std::string representative;
      while (std::getline(inStream, line)) {
        char *header;
        if (entry == 0) {
          header = queryReader.getDataByDBKey(
              (unsigned int)strtoul(line.c_str(), NULL, 10), thread_idx);

          representative = line;
        } else {
          header = targetReader.getDataByDBKey(
              (unsigned int)strtoul(line.c_str(), NULL, 10), thread_idx);
        }
        headers.emplace_back(header);
        entry++;
      }

      std::ostringstream oss;
      oss << par.summaryPrefix << "-" << representative << "|"
          << summarizer->summarize(out, headers);

      std::string summary = oss.str();
      writer.writeData(summary.c_str(), summary.length(), id, thread_idx);
    }
  }
  writer.close();
  reader.close();
  targetReader.close();
  queryReader.close();
  delete summarizer;
  return EXIT_SUCCESS;
}
