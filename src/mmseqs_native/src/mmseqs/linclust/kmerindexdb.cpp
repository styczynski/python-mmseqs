#include <mmseqs/output.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/linclust/kmerIndex.h>
#include <mmseqs/linclust/linsearchIndexReader.h>
#include <mmseqs/commons/nucleotideMatrix.h>
#include <mmseqs/prefiltering/prefilteringIndexReader.h>
#include <mmseqs/prefiltering/reducedMatrix.h>
#include <mmseqs/commons/substitutionMatrix.h>
#include <mmseqs/commons/timer.h>
#include <mmseqs/linclust/kmersearch.h>
#include <mmseqs/output.h>

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t)-1)
#endif

int kmerindexdb(mmseqs_output *out, Parameters &par) {
  //    MMseqsMPI::init(argc, argv);
  //
  //    Parameters &par = Parameters::getInstance();
  //    setLinearFilterDefault(&par);
  //    par.parseParameters(argc, argv, command, true, 0,
  //    MMseqsParameter::COMMAND_CLUSTLINEAR);

  DBReader<unsigned int> seqDbr(
      out, par.db1.c_str(), par.db1Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  seqDbr.open(DBReader<unsigned int>::NOSORT);
  int querySeqType = seqDbr.getDbtype();

  setKmerLengthAndAlphabet(par, seqDbr.getAminoAcidDBSize(), querySeqType);
  // par.printParameters(command.cmd, argc, argv, *command.params);

  out->info("Database size: {} type: {}", seqDbr.getSize(), seqDbr.getDbTypeName());
  std::string indexDB = LinsearchIndexReader::indexName(out, par.db2);
  if (par.checkCompatible > 0 && FileUtil::fileExists(out, indexDB.c_str())) {
    out->info("Check index {}", indexDB);
    DBReader<unsigned int> index(
        out, indexDB.c_str(), (indexDB + ".index").c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    index.open(DBReader<unsigned int>::NOSORT);

    if (Parameters::isEqualDbtype(seqDbr.getDbtype(),
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        par.PARAM_ALPH_SIZE.wasSet) {
      out->warn("Alphabet size is not taken into account for compatibility check in nucleotide search");
    }

    std::string check;
    const bool compatible =
        LinsearchIndexReader::checkIfIndexFile(out, &index) &&
        (check = LinsearchIndexReader::findIncompatibleParameter(
             out, index, par, seqDbr.getDbtype())) == "";
    index.close();
    seqDbr.close();
    if (compatible) {
      out->info("Index is already up to date and compatible. Force recreation with --check-compatibility 0 parameter");
      return EXIT_SUCCESS;
    } else {
      if (par.checkCompatible == 2) {
        out->error("Index is incompatible. Incompatible parameter: {}", check);
        return EXIT_FAILURE;
      } else {
        out->warn("Index is incompatible and will be recreated. Incompatible parameter: {}", check);
      }
    }
  }

  BaseMatrix *subMat;
  if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
    subMat =
        new NucleotideMatrix(out, par.seedScoringMatrixFile.nucleotides, 1.0, 0.0);
  } else {
    if (par.alphabetSize.aminoacids == 21) {
      subMat = new SubstitutionMatrix(out, par.seedScoringMatrixFile.aminoacids, 2.0,
                                      0.0);
    } else {
      SubstitutionMatrix sMat(out, par.seedScoringMatrixFile.aminoacids, 2.0, 0.0);
      subMat = new ReducedMatrix(out, sMat.probMatrix, sMat.subMatrixPseudoCounts,
                                 sMat.aa2num, sMat.num2aa, sMat.alphabetSize,
                                 par.alphabetSize.aminoacids, 2.0);
    }
  }

  // seqDbr.readMmapedDataInMemory();

  // memoryLimit in bytes
  size_t memoryLimit = Util::computeMemory(out, par.splitMemoryLimit);

  out->info("\n");

  float kmersPerSequenceScale =
      (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES))
          ? par.kmersPerSequenceScale.nucleotides
          : par.kmersPerSequenceScale.aminoacids;
  size_t totalKmers = computeKmerCount(
      seqDbr, par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
  totalKmers *= par.pickNbest;
  size_t totalSizeNeeded = computeMemoryNeededLinearfilter<short>(totalKmers);
  // compute splits
  size_t splits = static_cast<size_t>(
      std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
  size_t totalKmersPerSplit =
      std::max(static_cast<size_t>(1024 + 1),
               static_cast<size_t>(std::min(totalSizeNeeded, memoryLimit) /
                                   sizeof(KmerPosition<short>)) +
                   1);
  std::vector<std::pair<size_t, size_t>> hashRanges =
      setupKmerSplits<short>(out, par, subMat, seqDbr, totalKmersPerSplit, splits);

  out->info("Process file into {} parts\n", hashRanges.size());
  std::vector<std::string> splitFiles;
  KmerPosition<short> *hashSeqPair = NULL;

  size_t writePos = 0;
  size_t mpiRank = 0;
  size_t adjustedKmerSize = par.kmerSize;
#ifdef HAVE_MPI
  splits = std::max(static_cast<size_t>(MMseqsMPI::numProc), splits);
  size_t fromSplit = 0;
  size_t splitCount = 1;
  mpiRank = MMseqsMPI::rank;
  // if split size is great than nodes than we have to
  // distribute all splits equally over all nodes
  unsigned int *splitCntPerProc = new unsigned int[MMseqsMPI::numProc];
  memset(splitCntPerProc, 0, sizeof(unsigned int) * MMseqsMPI::numProc);
  for (size_t i = 0; i < splits; i++) {
    splitCntPerProc[i % MMseqsMPI::numProc] += 1;
  }
  for (int i = 0; i < MMseqsMPI::rank; i++) {
    fromSplit += splitCntPerProc[i];
  }
  splitCount = splitCntPerProc[MMseqsMPI::rank];
  delete[] splitCntPerProc;

  for (size_t split = fromSplit; split < fromSplit + splitCount; split++) {
    std::string splitFileName = par.db2 + "_split_" + SSTR(split);
    size_t splitKmerCount =
        (splits > 1) ? static_cast<size_t>(
                           static_cast<double>(totalKmers / splits) * 1.2)
                     : totalKmers;
    int range =
        MathUtil::ceilIntDivision(USHRT_MAX + 1, static_cast<int>(splits));
    size_t rangeFrom = split * range;
    size_t rangeTo = (splits == 1) ? SIZE_T_MAX : splits * range + range;
    KmerSearch::ExtractKmerAndSortResult kmerRet =
        KmerSearch::extractKmerAndSort(splitKmerCount, rangeFrom, rangeTo,
                                       seqDbr, par, subMat);
    hashSeqPair = kmerRet.kmers;
    // assign rep. sequence to same kmer members
    // The longest sequence is the first since we sorted by kmer, seq.Len and id
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(),
                                  Parameters::DBTYPE_NUCLEOTIDES)) {
      writePos =
          LinsearchIndexReader::pickCenterKmer<Parameters::DBTYPE_NUCLEOTIDES>(
              hashSeqPair, splitKmerCount);
    } else {
      writePos =
          LinsearchIndexReader::pickCenterKmer<Parameters::DBTYPE_AMINO_ACIDS>(
              hashSeqPair, splitKmerCount);
    }

    LinsearchIndexReader::writeKmerIndexToDisk(splitFileName, hashSeqPair,
                                               writePos);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpiRank == 0) {
    for (size_t split = 0; split < splits; split++) {
      std::string splitFileName = par.db2 + "_split_" + SSTR(split);
      splitFiles.push_back(splitFileName);
    }
  }
#else
  for (size_t split = 0; split < hashRanges.size(); split++) {
    out->info("Generate k-mers list {}\n", split);

    std::string splitFileName = par.db2 + "_split_" + SSTR(split);

    KmerSearch::ExtractKmerAndSortResult kmerRet =
        KmerSearch::extractKmerAndSort(
            out,
            totalKmersPerSplit, hashRanges[split].first,
            hashRanges[split].second, seqDbr, par, subMat);
    hashSeqPair = kmerRet.kmers;
    adjustedKmerSize = std::max(adjustedKmerSize, kmerRet.adjustedKmer);
    // assign rep. sequence to same kmer members
    // The longest sequence is the first since we sorted by kmer, seq.Len and id
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(),
                                  Parameters::DBTYPE_NUCLEOTIDES)) {
      writePos =
          LinsearchIndexReader::pickCenterKmer<Parameters::DBTYPE_NUCLEOTIDES>(
              out, hashSeqPair, totalKmersPerSplit);
    } else {
      writePos =
          LinsearchIndexReader::pickCenterKmer<Parameters::DBTYPE_AMINO_ACIDS>(
              out, hashSeqPair, totalKmersPerSplit);
    }

    if (splits > 1) {
      LinsearchIndexReader::writeKmerIndexToDisk(out, splitFileName, hashSeqPair,
                                                 writePos);
      delete[] hashSeqPair;
      hashSeqPair = NULL;
    }

    splitFiles.push_back(splitFileName);
  }
#endif
  if (mpiRank == 0) {
    // write result
    DBWriter dbw(out, indexDB.c_str(), (indexDB + ".index").c_str(), 1,
                 par.compressed, Parameters::DBTYPE_INDEX_DB);
    dbw.open();

    out->info("Write VERSION ({})\n", PrefilteringIndexReader::VERSION
                      );
    dbw.writeData(
        (char *)PrefilteringIndexReader::CURRENT_VERSION,
        strlen(PrefilteringIndexReader::CURRENT_VERSION) * sizeof(char),
        PrefilteringIndexReader::VERSION, 0);
    dbw.alignToPageSize();

    out->info("Write META ({})\n", PrefilteringIndexReader::META
                      );
    const int mask = par.maskMode > 0;
    const int spacedKmer = (par.spacedKmer) ? 1 : 0;
    const bool sameDB = (par.db1 == par.db2);
    const int headers1 = 1;
    const int headers2 = (sameDB) ? 1 : 0;
    const int seqType = seqDbr.getDbtype();
    const int srcSeqType = FileUtil::parseDbType(out, par.db2.c_str());
    // Reuse the compBiasCorr field to store the adjustedKmerSize, It is not
    // needed in the linsearch
    int metadata[] = {static_cast<int>(par.maxSeqLen),
                      static_cast<int>(par.kmerSize),
                      static_cast<int>(adjustedKmerSize),
                      subMat->alphabetSize,
                      mask,
                      spacedKmer,
                      0,
                      seqType,
                      srcSeqType,
                      headers1,
                      headers2};
    char *metadataptr = (char *)&metadata;
    dbw.writeData(metadataptr, sizeof(metadata), PrefilteringIndexReader::META,
                  0);
    dbw.alignToPageSize();

    Timer timer;
    if (splits > 1) {
      seqDbr.unmapData();
      if (Parameters::isEqualDbtype(seqDbr.getDbtype(),
                                    Parameters::DBTYPE_NUCLEOTIDES)) {
        LinsearchIndexReader::mergeAndWriteIndex<
            Parameters::DBTYPE_NUCLEOTIDES>(
            out, dbw, splitFiles, subMat->alphabetSize, adjustedKmerSize);
      } else {
        LinsearchIndexReader::mergeAndWriteIndex<
            Parameters::DBTYPE_AMINO_ACIDS>(
            out, dbw, splitFiles, subMat->alphabetSize, adjustedKmerSize);
      }
    } else {
      if (Parameters::isEqualDbtype(seqDbr.getDbtype(),
                                    Parameters::DBTYPE_NUCLEOTIDES)) {
        LinsearchIndexReader::writeIndex<Parameters::DBTYPE_NUCLEOTIDES>(
            out, dbw, hashSeqPair, writePos, subMat->alphabetSize, adjustedKmerSize);
      } else {
        LinsearchIndexReader::writeIndex<Parameters::DBTYPE_AMINO_ACIDS>(
            out, dbw, hashSeqPair, writePos, subMat->alphabetSize, adjustedKmerSize);
      }
    }
    if (hashSeqPair) {
      delete[] hashSeqPair;
      hashSeqPair = NULL;
    }
    // SEQCOUNT
    out->info("Write SEQCOUNT ({})", PrefilteringIndexReader::SEQCOUNT);
    size_t tablesize = {seqDbr.getSize()};
    char *tablesizePtr = (char *)&tablesize;
    dbw.writeData(tablesizePtr, 1 * sizeof(size_t),
                  PrefilteringIndexReader::SEQCOUNT, 0);
    dbw.alignToPageSize();

    out->info("Write SCOREMATRIXNAME ({})", PrefilteringIndexReader::SCOREMATRIXNAME);
    char *subData =
        BaseMatrix::serialize(subMat->matrixName, subMat->matrixData);
    dbw.writeData(
        subData, BaseMatrix::memorySize(subMat->matrixName, subMat->matrixData),
        PrefilteringIndexReader::SCOREMATRIXNAME, 0);
    dbw.alignToPageSize();
    free(subData);

    if (par.spacedKmerPattern.empty() != false) {
      out->info("Write SPACEDPATTERN ({})", PrefilteringIndexReader::SPACEDPATTERN);
      dbw.writeData(par.spacedKmerPattern.c_str(),
                    par.spacedKmerPattern.length(),
                    PrefilteringIndexReader::SPACEDPATTERN, 0);
      dbw.alignToPageSize();
    }

    seqDbr.close();

    DBReader<unsigned int> dbr1(
        out, par.db1.c_str(), par.db1Index.c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    dbr1.open(DBReader<unsigned int>::NOSORT);
    out->info("Write DBR1INDEX ({})", PrefilteringIndexReader::DBR1INDEX);
    char *data = DBReader<unsigned int>::serialize(dbr1);
    size_t offsetIndex = dbw.getOffset(0);
    dbw.writeData(data, DBReader<unsigned int>::indexMemorySize(dbr1),
                  PrefilteringIndexReader::DBR1INDEX, 0);
    dbw.alignToPageSize();

    out->info("Write DBR1DATA ({})", PrefilteringIndexReader::DBR1DATA);
    size_t offsetData = dbw.getOffset(0);
    dbw.writeStart(0);
    for (size_t fileIdx = 0; fileIdx < dbr1.getDataFileCnt(); fileIdx++) {
      dbw.writeAdd(dbr1.getDataForFile(fileIdx),
                   dbr1.getDataSizeForFile(fileIdx), 0);
    }
    dbw.writeEnd(PrefilteringIndexReader::DBR1DATA, 0);
    dbw.alignToPageSize();
    free(data);

    if (sameDB == true) {
      dbw.writeIndexEntry(PrefilteringIndexReader::DBR2INDEX, offsetIndex,
                          DBReader<unsigned int>::indexMemorySize(dbr1) + 1, 0);
      dbw.writeIndexEntry(PrefilteringIndexReader::DBR2DATA, offsetData,
                          dbr1.getTotalDataSize() + 1, 0);
      dbr1.close();
    } else {
      dbr1.close();
      DBReader<unsigned int> dbr2(
          out, par.db2.c_str(), par.db2Index.c_str(), par.threads,
          DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
      dbr2.open(DBReader<unsigned int>::NOSORT);
      out->info("Write DBR2INDEX ({})", PrefilteringIndexReader::DBR2INDEX);
      data = DBReader<unsigned int>::serialize(dbr2);
      dbw.writeData(data, DBReader<unsigned int>::indexMemorySize(dbr2),
                    PrefilteringIndexReader::DBR2INDEX, 0);
      dbw.alignToPageSize();
      out->info("Write DBR2DATA ({})", PrefilteringIndexReader::DBR2DATA);
      dbw.writeStart(0);
      for (size_t fileIdx = 0; fileIdx < dbr2.getDataFileCnt(); fileIdx++) {
        dbw.writeAdd(dbr2.getDataForFile(fileIdx),
                     dbr2.getDataSizeForFile(fileIdx), 0);
      }
      dbw.writeEnd(PrefilteringIndexReader::DBR2DATA, 0);
      dbw.alignToPageSize();
      free(data);
      dbr2.close();
    }

    {
      out->info("Write HDR1INDEX ({})", PrefilteringIndexReader::HDR1INDEX);
      DBReader<unsigned int> hdbr1(
          out, par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads,
          DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
      hdbr1.open(DBReader<unsigned int>::NOSORT);

      data = DBReader<unsigned int>::serialize(hdbr1);
      size_t offsetIndex = dbw.getOffset(0);
      dbw.writeData(data, DBReader<unsigned int>::indexMemorySize(hdbr1),
                    PrefilteringIndexReader::HDR1INDEX, 0);
      dbw.alignToPageSize();
      out->info("Write HDR1DATA ({})", PrefilteringIndexReader::HDR1DATA);
      size_t offsetData = dbw.getOffset(0);
      dbw.writeStart(0);
      for (size_t fileIdx = 0; fileIdx < hdbr1.getDataFileCnt(); fileIdx++) {
        dbw.writeAdd(hdbr1.getDataForFile(fileIdx),
                     hdbr1.getDataSizeForFile(fileIdx), 0);
      }
      dbw.writeEnd(PrefilteringIndexReader::HDR1DATA, 0);
      dbw.alignToPageSize();
      free(data);
      if (sameDB == true) {
        dbw.writeIndexEntry(PrefilteringIndexReader::HDR2INDEX, offsetIndex,
                            DBReader<unsigned int>::indexMemorySize(hdbr1) + 1,
                            0);
        dbw.writeIndexEntry(PrefilteringIndexReader::HDR2DATA, offsetData,
                            hdbr1.getTotalDataSize() + 1, 0);
        hdbr1.close();
      } else {
        hdbr1.close();
        DBReader<unsigned int> hdbr2(out, par.hdr2.c_str(), par.hdr2Index.c_str(),
                                     par.threads,
                                     DBReader<unsigned int>::USE_INDEX |
                                         DBReader<unsigned int>::USE_DATA);
        hdbr2.open(DBReader<unsigned int>::NOSORT);
        out->info("Write HDR2INDEX ({})", PrefilteringIndexReader::HDR2INDEX);
        data = DBReader<unsigned int>::serialize(hdbr2);
        dbw.writeData(data, DBReader<unsigned int>::indexMemorySize(hdbr2),
                      PrefilteringIndexReader::HDR2INDEX, 0);
        dbw.alignToPageSize();
        out->info("Write HDR2DATA ({})", PrefilteringIndexReader::HDR2DATA);
        dbw.writeStart(0);
        for (size_t fileIdx = 0; fileIdx < hdbr2.getDataFileCnt(); fileIdx++) {
          dbw.writeAdd(hdbr2.getDataForFile(fileIdx),
                       hdbr2.getDataSizeForFile(fileIdx), 0);
        }
        dbw.writeEnd(PrefilteringIndexReader::HDR2DATA, 0);
        dbw.alignToPageSize();
        hdbr2.close();
        free(data);
      }
    }

    out->info("Write GENERATOR ({})", PrefilteringIndexReader::GENERATOR);
    dbw.writeData(version, strlen(version), PrefilteringIndexReader::GENERATOR,
                  0);
    dbw.alignToPageSize();

    out->info("Time for fill: {}\n", timer.lap());
    // add missing entries to the result (needed for clustering)
    dbw.close();
  }
  // free memory
  delete subMat;
  if (mpiRank != 0) {
    if (hashSeqPair) {
      delete[] hashSeqPair;
    }
    seqDbr.close();
  }

  return EXIT_SUCCESS;
}
