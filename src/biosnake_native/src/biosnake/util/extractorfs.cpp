#include <biosnake/commons/dBReader.h>
#include <biosnake/commons/dBWriter.h>
#include <biosnake/output.h>
#include <biosnake/alignment/matcher.h>
#include <biosnake/commons/orf.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/translateNucl.h>
#include <biosnake/commons/util.h>
#include <biosnake/commons/itoa.h>
#include <biosnake/output.h>

#include <unistd.h>
#include <algorithm>
#include <climits>

#ifdef OPENMP
#include <omp.h>
#endif

int extractorfs(biosnake_output* out, Parameters& par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  DBReader<unsigned int> reader(
      out, par.db1.c_str(), par.db1Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  reader.open(DBReader<unsigned int>::NOSORT);

  DBReader<unsigned int> headerReader(
      out, par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  int outputDbtype = Parameters::DBTYPE_NUCLEOTIDES;
  headerReader.open(DBReader<unsigned int>::NOSORT);
  if (par.translate) {
    outputDbtype = Parameters::DBTYPE_AMINO_ACIDS;
  }
  DBWriter sequenceWriter(out, par.db2.c_str(), par.db2Index.c_str(), par.threads,
                          par.compressed, outputDbtype);
  sequenceWriter.open();

  DBWriter headerWriter(out, par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads,
                        false, Parameters::DBTYPE_GENERIC_DB);
  headerWriter.open();

  if ((par.orfStartMode == 1) && (par.contigStartMode < 2)) {
    out->failure("Parameter combination is illegal, orf-start-mode 1 can only go with contig-start-mode 2");
  }

  unsigned int forwardFrames = Orf::getFrames(par.forwardFrames);
  unsigned int reverseFrames = Orf::getFrames(par.reverseFrames);
  const char newline = '\n';
  Log::Progress progress(reader.getSize());
  TranslateNucl translateNucl(
      out, static_cast<TranslateNucl::GenCode>(par.translationTable));

#pragma omp parallel
  {
    Orf orf(out, par.translationTable, par.useAllTableStarts);
    int thread_idx = 0;
#ifdef OPENMP
    thread_idx = omp_get_thread_num();
#endif
    size_t querySize = 0;
    size_t queryFrom = 0;
    reader.decomposeDomainByAminoAcid(thread_idx, par.threads, &queryFrom,
                                      &querySize);
    if (querySize == 0) {
      queryFrom = 0;
    }
    char* aa = new char[par.maxSeqLen + 3 + 1];
    char buffer[1024];

    std::vector<Orf::SequenceLocation> res;
    res.reserve(1000);
    for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i) {
      progress.updateProgress();

      unsigned int key = reader.getDbKey(i);
      const char* data = reader.getData(i, thread_idx);
      size_t sequenceLength = reader.getSeqLen(i);
      if (!orf.setSequence(data, sequenceLength)) {
        out->warn("Invalid sequence with index {}!", i);
        continue;
      }

      const char* header = headerReader.getData(i, thread_idx);
      std::string headerAccession = Util::parseFastaHeader(header);
      orf.findAll(res, par.orfMinLength, par.orfMaxLength, par.orfMaxGaps,
                  forwardFrames, reverseFrames, par.orfStartMode);
      for (std::vector<Orf::SequenceLocation>::const_iterator it = res.begin();
           it != res.end(); ++it) {
        Orf::SequenceLocation loc = *it;

        if (par.contigStartMode < 2 &&
            (loc.hasIncompleteStart == par.contigStartMode)) {
          continue;
        }
        if (par.contigEndMode < 2 &&
            (loc.hasIncompleteEnd == par.contigEndMode)) {
          continue;
        }

        std::pair<const char*, size_t> sequence = orf.getSequence(loc);
        size_t fromPos = loc.from;
        size_t toPos = loc.to;
        if (loc.strand == Orf::STRAND_MINUS) {
          fromPos = (sequenceLength - 1) - loc.from;
          toPos = (sequenceLength - 1) - loc.to;
        }
        Orf::writeOrfHeader(buffer, key, fromPos, toPos, loc.hasIncompleteStart,
                            loc.hasIncompleteEnd);

        //                snprintf(buffer, LINE_MAX, "%.*s [Orf: %d, %zu, %zu,
        //                %d, %d]\n", (unsigned int)(headerAccession.size()),
        //                headerAccession.c_str(),
        //                          toPos, loc.hasIncompleteStart,
        //                          loc.hasIncompleteEnd);
        sequenceWriter.writeStart(thread_idx);
        if (par.translate) {
          if ((data[sequence.second] != '\n' && sequence.second % 3 != 0) &&
              (data[sequence.second - 1] == '\n' &&
               (sequence.second - 1) % 3 != 0)) {
            sequence.second = sequence.second - (sequence.second % 3);
          }

          if (sequence.second < 3) {
            continue;
          }

          if (sequence.second > (3 * par.maxSeqLen)) {
            sequence.second = (3 * par.maxSeqLen);
          }
          translateNucl.translate(aa, sequence.first, sequence.second);
          sequenceWriter.writeAdd(aa, (sequence.second / 3), thread_idx);

        } else {
          sequenceWriter.writeAdd(sequence.first, sequence.second, thread_idx);
        }
        sequenceWriter.writeAdd(&newline, 1, thread_idx);
        sequenceWriter.writeEnd(key, thread_idx);

        headerWriter.writeData(buffer, strlen(buffer), key, thread_idx);
      }
      res.clear();
    }
    delete[] aa;
  }
  headerWriter.close(true);
  sequenceWriter.close(true);
  headerReader.close();
  reader.close();
  // make identifiers stable
#pragma omp parallel
  {
#pragma omp single
    {
#pragma omp task
      { DBWriter::createRenumberedDB(out, par.hdr2, par.hdr2Index, "", ""); }

#pragma omp task
      {
        DBWriter::createRenumberedDB(out, par.db2, par.db2Index,
                                     par.createLookup ? par.db1 : "",
                                     par.createLookup ? par.db1Index : "");
      }
    }
  }
  DBReader<unsigned int>::softlinkDb(out, par.db1, par.db2, DBFiles::SOURCE);

  return EXIT_SUCCESS;
}
