//
// Created by Martin Steinegger on 2019-01-04.
//

#include <mmseqs/linclust/kmersearch.h>
#include <mmseqs/commons/command.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/fastSort.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/linclust/kmerIndex.h>
#include <mmseqs/linclust/linsearchIndexReader.h>
#include <mmseqs/commons/nucleotideMatrix.h>
#include <mmseqs/prefiltering/prefilteringIndexReader.h>
#include <mmseqs/prefiltering/queryMatcher.h>
#include <mmseqs/prefiltering/reducedMatrix.h>
#include <mmseqs/commons/substitutionMatrix.h>
#include <mmseqs/commons/timer.h>
#include <mmseqs/output.h>

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t)-1)
#endif

KmerSearch::ExtractKmerAndSortResult KmerSearch::extractKmerAndSort(
    size_t totalKmers, size_t hashStartRange, size_t hashEndRange,
    DBReader<unsigned int> &seqDbr, Parameters &par, BaseMatrix *subMat) {
  KmerPosition<short> *hashSeqPair = initKmerPositionMemory<short>(totalKmers);
  Timer timer;
  size_t elementsToSort;
  if (par.pickNbest > 1) {
    std::pair<size_t, size_t> ret =
        fillKmerPositionArray<Parameters::DBTYPE_HMM_PROFILE, short>(
            hashSeqPair, totalKmers, seqDbr, par, subMat, false, hashStartRange,
            hashEndRange, NULL);
    elementsToSort = ret.first;
  } else if (Parameters::isEqualDbtype(seqDbr.getDbtype(),
                                       Parameters::DBTYPE_NUCLEOTIDES)) {
    std::pair<size_t, size_t> ret =
        fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, short>(
            hashSeqPair, totalKmers, seqDbr, par, subMat, false, hashStartRange,
            hashEndRange, NULL);
    elementsToSort = ret.first;
    par.kmerSize = ret.second;
    out->info("\nAdjusted k-mer length {}\n", par.kmerSize);
  } else {
    std::pair<size_t, size_t> ret =
        fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, short>(
            hashSeqPair, totalKmers, seqDbr, par, subMat, false, hashStartRange,
            hashEndRange, NULL);
    elementsToSort = ret.first;
  }
  out->info("\nTime for fill: {}\n", timer.lap());
  if (hashEndRange == SIZE_T_MAX) {
    seqDbr.unmapData();
  }

  out->info("Sort kmer ... ");
  timer.reset();
  if (Parameters::isEqualDbtype(seqDbr.getDbtype(),
                                Parameters::DBTYPE_NUCLEOTIDES)) {
    SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort,
                  KmerPosition<short>::compareRepSequenceAndIdAndPosReverse);
  } else {
    SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort,
                  KmerPosition<short>::compareRepSequenceAndIdAndPos);
  }

  out->info("Time for sort: {}\n", timer.lap());

  return ExtractKmerAndSortResult(elementsToSort, hashSeqPair, par.kmerSize);
}

template <int TYPE>
void KmerSearch::writeResult(DBWriter &dbw, KmerPosition<short> *kmers,
                             size_t kmerCount) {
  size_t repSeqId = SIZE_T_MAX;
  unsigned int prevHitId;
  char buffer[100];
  std::string prefResultsOutString;
  prefResultsOutString.reserve(100000000);
  for (size_t i = 0; i < kmerCount; i++) {
    size_t currId = kmers[i].kmer;
    int reverMask = 0;
    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
      reverMask = BIT_CHECK(kmers[i].kmer, 63) == false;
      currId = BIT_CLEAR(currId, 63);
    }
    if (repSeqId != currId) {
      if (repSeqId != SIZE_T_MAX) {
        dbw.writeData(prefResultsOutString.c_str(),
                      prefResultsOutString.length(),
                      static_cast<unsigned int>(repSeqId), 0);
      }
      repSeqId = currId;
      prefResultsOutString.clear();
    }
    //        std::cout << kmers[i].id << "\t" << kmers[i].pos << std::endl;
    // find maximal diagonal and top score
    short prevDiagonal;
    int cnt = 0;
    int bestDiagonalCnt = 0;
    int bestRevertMask = reverMask;
    short bestDiagonal = kmers[i].pos;
    int topScore = 0;
    unsigned int tmpCurrId = currId;

    unsigned int hitId;
    do {
      prevHitId = kmers[i].id;
      prevDiagonal = kmers[i].pos;

      cnt = (kmers[i].pos == prevDiagonal) ? cnt + 1 : 1;
      if (cnt > bestDiagonalCnt) {
        bestDiagonalCnt = cnt;
        bestDiagonal = kmers[i].pos;
        bestRevertMask = reverMask;
      }
      topScore++;
      i++;
      hitId = kmers[i].id;
      tmpCurrId = kmers[i].kmer;
      if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
        reverMask = BIT_CHECK(kmers[i].kmer, 63) == false;
        tmpCurrId = BIT_CLEAR(tmpCurrId, 63);
      }
    } while (hitId == prevHitId && currId == tmpCurrId && i < kmerCount);
    i--;

    hit_t h;
    h.seqId = prevHitId;
    h.prefScore = (bestRevertMask) ? -topScore : topScore;
    h.diagonal = bestDiagonal;
    int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
    prefResultsOutString.append(buffer, len);
  }
  // last element
  if (prefResultsOutString.size() > 0) {
    if (repSeqId != SIZE_T_MAX) {
      dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(),
                    static_cast<unsigned int>(repSeqId), 0);
    }
  }
}

template void KmerSearch::writeResult<0>(DBWriter &dbw,
                                         KmerPosition<short> *kmers,
                                         size_t kmerCount);
template void KmerSearch::writeResult<1>(DBWriter &dbw,
                                         KmerPosition<short> *kmers,
                                         size_t kmerCount);

int kmersearch(mmseqs_output *out, Parameters &par) {
  //    Parameters &par = Parameters::getInstance();
  //    setLinearFilterDefault(&par);
  //    par.parseParameters(argc, argv, command, true, 0,
  //    MMseqsParameter::COMMAND_CLUSTLINEAR);
  int targetSeqType;
  int adjustedKmerSize = 0;
  if (Parameters::isEqualDbtype(FileUtil::parseDbType(par.db2.c_str()),
                                Parameters::DBTYPE_INDEX_DB) == false) {
    out->failure("Create index before calling kmersearch with mmseqs createlinindex");
  }

  DBReader<unsigned int> tidxdbr(
      par.db2.c_str(), par.db2Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  tidxdbr.open(DBReader<unsigned int>::NOSORT);
  PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&tidxdbr);
  if (par.PARAM_K.wasSet) {
    if (par.kmerSize != 0 && data.kmerSize != par.kmerSize) {
      out->failure("Index was created with -k {} but the prefilter was called with -k {}. Please execute createindex -k {}", data.kmerSize, par.kmerSize, par.kmerSize);
    }
  }
  if (par.PARAM_ALPH_SIZE.wasSet) {
    if (data.alphabetSize !=
        (Parameters::isEqualDbtype(data.seqType, Parameters::DBTYPE_AMINO_ACIDS)
             ? par.alphabetSize.aminoacids
             : par.alphabetSize.nucleotides)) {
      out->failure("Index was created with --alph-size {} but the prefilter was called with --alph-size {}. Please execute createindex --alph-size {}", data.alphabetSize, MultiParam<int>::format(par.alphabetSize), MultiParam<int>::format(par.alphabetSize));
    }
  }
  if (par.PARAM_SPACED_KMER_MODE.wasSet) {
    if (data.spacedKmer != par.spacedKmer) {
      out->failure("Index was created with --spaced-kmer-mode {} but the prefilter was called with --spaced-kmer-mode {}. Please execute createindex --spaced-kmer-mode {}", data.spacedKmer, par.spacedKmer, par.spacedKmer);
    }
  }
  par.kmerSize = data.kmerSize;
  par.alphabetSize = data.alphabetSize;
  targetSeqType = data.seqType;
  par.spacedKmer = (data.spacedKmer == 1) ? true : false;
  par.maxSeqLen = data.maxSeqLength;
  // Reuse the compBiasCorr field to store the adjustedKmerSize, It is not
  // needed in the linsearch
  adjustedKmerSize = data.compBiasCorr;

  DBReader<unsigned int> queryDbr(
      par.db1.c_str(), par.db1Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  queryDbr.open(DBReader<unsigned int>::NOSORT);
  int querySeqType = queryDbr.getDbtype();
  if (Parameters::isEqualDbtype(querySeqType, targetSeqType) == false) {
    out->failure("Dbtype of query and target database do not match");
  }

  setKmerLengthAndAlphabet(par, queryDbr.getAminoAcidDBSize(), querySeqType);

  // par.printParameters(command.cmd, argc, argv, *command.params);

  // queryDbr.readMmapedDataInMemory();

  // memoryLimit in bytes
  size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);

  float kmersPerSequenceScale =
      (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES))
          ? par.kmersPerSequenceScale.nucleotides
          : par.kmersPerSequenceScale.aminoacids;
  size_t totalKmers = computeKmerCount(
      queryDbr, par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
  size_t totalSizeNeeded = computeMemoryNeededLinearfilter<short>(totalKmers);

  BaseMatrix *subMat;
  if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
    subMat =
        new NucleotideMatrix(par.seedScoringMatrixFile.nucleotides, 1.0, 0.0);
  } else {
    if (par.alphabetSize.aminoacids == 21) {
      subMat = new SubstitutionMatrix(par.seedScoringMatrixFile.aminoacids, 8.0,
                                      -0.2);
    } else {
      SubstitutionMatrix sMat(par.seedScoringMatrixFile.aminoacids, 8.0, -0.2);
      subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts,
                                 sMat.aa2num, sMat.num2aa, sMat.alphabetSize,
                                 par.alphabetSize.aminoacids, 8.0);
    }
  }

  // compute splits
  size_t splits = static_cast<size_t>(
      std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
  size_t totalKmersPerSplit =
      std::max(static_cast<size_t>(1024 + 1),
               static_cast<size_t>(std::min(totalSizeNeeded, memoryLimit) /
                                   sizeof(KmerPosition<short>)) +
                   1);

  std::vector<std::pair<size_t, size_t>> hashRanges =
      setupKmerSplits<short>(par, subMat, queryDbr, totalKmersPerSplit, splits);

  int outDbType = (Parameters::isEqualDbtype(queryDbr.getDbtype(),
                                             Parameters::DBTYPE_NUCLEOTIDES))
                      ? Parameters::DBTYPE_PREFILTER_REV_RES
                      : Parameters::DBTYPE_PREFILTER_RES;
  out->info("Process file into {} parts\n", hashRanges.size());

  std::vector<std::string> splitFiles;
  for (size_t split = 0; split < hashRanges.size(); split++) {
    tidxdbr.remapData();
    char *entriesData = tidxdbr.getDataUncompressed(
        tidxdbr.getId(PrefilteringIndexReader::ENTRIES));
    char *entriesOffsetsData = tidxdbr.getDataUncompressed(
        tidxdbr.getId(PrefilteringIndexReader::ENTRIESOFFSETS));
    int64_t entriesNum = *((int64_t *)tidxdbr.getDataUncompressed(
        tidxdbr.getId(PrefilteringIndexReader::ENTRIESNUM)));
    int64_t entriesGridSize = *((int64_t *)tidxdbr.getDataUncompressed(
        tidxdbr.getId(PrefilteringIndexReader::ENTRIESGRIDSIZE)));
    int alphabetSize =
        (Parameters::isEqualDbtype(queryDbr.getDbtype(),
                                   Parameters::DBTYPE_NUCLEOTIDES))
            ? par.alphabetSize.nucleotides
            : par.alphabetSize.aminoacids;
    KmerIndex kmerIndex(alphabetSize, adjustedKmerSize, entriesData,
                        entriesOffsetsData, entriesNum, entriesGridSize);
    //        kmerIndex.printIndex<Parameters::DBTYPE_NUCLEOTIDES>(subMat);
    std::pair<std::string, std::string> tmpFiles;
    if (splits > 1) {
      tmpFiles = Util::createTmpFileNames(par.db3.c_str(), par.db3Index.c_str(),
                                          split);
    } else {
      tmpFiles = std::make_pair(par.db3, par.db3Index);
    }
    splitFiles.push_back(tmpFiles.first);

    std::string splitFileNameDone = tmpFiles.first + ".done";
    if (FileUtil::fileExists(splitFileNameDone.c_str()) == false) {
      KmerSearch::ExtractKmerAndSortResult sortedKmers =
          KmerSearch::extractKmerAndSort(
              totalKmersPerSplit, hashRanges[split].first,
              hashRanges[split].second, queryDbr, par, subMat);
      std::pair<KmerPosition<short> *, size_t> result;
      if (Parameters::isEqualDbtype(queryDbr.getDbtype(),
                                    Parameters::DBTYPE_NUCLEOTIDES)) {
        result = KmerSearch::searchInIndex<Parameters::DBTYPE_NUCLEOTIDES>(
            sortedKmers.kmers, sortedKmers.kmerCount, kmerIndex,
            par.resultDirection);
      } else {
        result = KmerSearch::searchInIndex<Parameters::DBTYPE_AMINO_ACIDS>(
            sortedKmers.kmers, sortedKmers.kmerCount, kmerIndex,
            par.resultDirection);
      }

      KmerPosition<short> *kmers = result.first;
      size_t kmerCount = result.second;
      if (splits == 1) {
        DBWriter dbw(tmpFiles.first.c_str(), tmpFiles.second.c_str(), 1,
                     par.compressed, outDbType);
        dbw.open();
        if (Parameters::isEqualDbtype(queryDbr.getDbtype(),
                                      Parameters::DBTYPE_NUCLEOTIDES)) {
          KmerSearch::writeResult<Parameters::DBTYPE_NUCLEOTIDES>(dbw, kmers,
                                                                  kmerCount);
        } else {
          KmerSearch::writeResult<Parameters::DBTYPE_AMINO_ACIDS>(dbw, kmers,
                                                                  kmerCount);
        }
        dbw.close();
      } else {
        if (Parameters::isEqualDbtype(queryDbr.getDbtype(),
                                      Parameters::DBTYPE_NUCLEOTIDES)) {
          writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short>(
              tmpFiles.first, kmers, kmerCount);
        } else {
          writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short>(
              tmpFiles.first, kmers, kmerCount);
        }
      }
      delete[] kmers;
    }
  }
  delete subMat;
  tidxdbr.close();
  queryDbr.close();
  if (splitFiles.size() > 1) {
    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed,
                    outDbType);
    writer.open();  // 1 GB buffer
    std::vector<char> empty;
    if (Parameters::isEqualDbtype(querySeqType,
                                  Parameters::DBTYPE_NUCLEOTIDES)) {
      mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev>(
          writer, splitFiles, empty);
    } else {
      mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry>(
          writer, splitFiles, empty);
    }
    for (size_t i = 0; i < splitFiles.size(); i++) {
      FileUtil::remove(splitFiles[i].c_str());
      std::string splitFilesDone = splitFiles[i] + ".done";
      FileUtil::remove(splitFilesDone.c_str());
    }
    writer.close();
  }
  return EXIT_SUCCESS;
}
template <int TYPE>
std::pair<KmerPosition<short> *, size_t> KmerSearch::searchInIndex(
    KmerPosition<short> *kmers, size_t kmersSize, KmerIndex &kmerIndex,
    int resultDirection) {
  Timer timer;
  bool queryTargetSwitched =
      (resultDirection == Parameters::PARAM_RESULT_DIRECTION_TARGET);
  kmerIndex.reset();
  KmerIndex::KmerEntry currTargetKmer;
  bool isDone = false;
  if (kmerIndex.hasNextEntry()) {
    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
      currTargetKmer = kmerIndex.getNextEntry<Parameters::DBTYPE_NUCLEOTIDES>();
    } else {
      currTargetKmer = kmerIndex.getNextEntry<Parameters::DBTYPE_AMINO_ACIDS>();
    }
  } else {
    isDone = true;
  }

  size_t kmerPos = 0;
  size_t writePos = 0;
  // this is IO bound, optimisation does not make much sense here.
  size_t queryKmer;
  size_t targetKmer;

  while (isDone == false) {
    KmerPosition<short> *currQueryKmer = &kmers[kmerPos];
    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
      queryKmer = BIT_SET(currQueryKmer->kmer, 63);
      targetKmer = BIT_SET(currTargetKmer.kmer, 63);
    } else {
      queryKmer = currQueryKmer->kmer;
      targetKmer = currTargetKmer.kmer;
    }

    if (queryKmer < targetKmer) {
      while (queryKmer < targetKmer) {
        if (kmerPos + 1 < kmersSize) {
          kmerPos++;
        } else {
          isDone = true;
          break;
        }
        KmerPosition<short> *currQueryKmer = &kmers[kmerPos];
        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
          queryKmer = BIT_SET(currQueryKmer->kmer, 63);
        } else {
          queryKmer = currQueryKmer->kmer;
        }
      }
    } else if (targetKmer < queryKmer) {
      while (targetKmer < queryKmer) {
        if (kmerIndex.hasNextEntry()) {
          if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
            currTargetKmer =
                kmerIndex.getNextEntry<Parameters::DBTYPE_NUCLEOTIDES>();
          } else {
            currTargetKmer =
                kmerIndex.getNextEntry<Parameters::DBTYPE_AMINO_ACIDS>();
          }
          if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
            targetKmer = BIT_SET(currTargetKmer.kmer, 63);
          } else {
            targetKmer = currTargetKmer.kmer;
          }
          // TODO remap logic to speed things up
        } else {
          isDone = true;
          break;
        }
      }
    } else {
      if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
        //  00 No problem here both are forward
        //  01 We can revert the query of target, lets invert the query.
        //  10 Same here, we can revert query to match the not inverted target
        //  11 Both are reverted so no problem!
        //  So we need just 1 bit of information to encode all four states
        bool targetIsReverse =
            (queryTargetSwitched)
                ? (BIT_CHECK(currQueryKmer->kmer, 63) == false)
                : (BIT_CHECK(currTargetKmer.kmer, 63) == false);
        bool repIsReverse = (queryTargetSwitched)
                                ? (BIT_CHECK(currTargetKmer.kmer, 63) == false)
                                : (BIT_CHECK(currQueryKmer->kmer, 63) == false);
        bool queryNeedsToBeRev = false;
        // we now need 2 byte of information (00),(01),(10),(11)
        // we need to flip the coordinates of the query
        short queryPos = currTargetKmer.pos;
        short targetPos = currQueryKmer->pos;
        ;
        // revert kmer in query hits normal kmer in target
        // we need revert the query
        if (repIsReverse == true && targetIsReverse == false) {
          queryNeedsToBeRev = true;
          // both k-mers were extracted on the reverse strand
          // this is equal to both are extract on the forward strand
          // we just need to offset the position to the forward strand
        } else if (repIsReverse == true && targetIsReverse == true) {
          queryPos = (currTargetKmer.seqLen - 1) - currTargetKmer.pos;
          targetPos = (currQueryKmer->seqLen - 1) - currQueryKmer->pos;
          queryNeedsToBeRev = false;
          // query is not revers but target k-mer is reverse
          // instead of reverting the target, we revert the query and offset the
          // the query/target position
        } else if (repIsReverse == false && targetIsReverse == true) {
          queryPos = (currTargetKmer.seqLen - 1) - currTargetKmer.pos;
          targetPos = (currQueryKmer->seqLen - 1) - currQueryKmer->pos;
          queryNeedsToBeRev = true;
          // both are forward, everything is good here
        }
        (kmers + writePos)->pos =
            (queryTargetSwitched) ? queryPos - targetPos : targetPos - queryPos;
        size_t id =
            (queryTargetSwitched) ? currTargetKmer.id : currQueryKmer->id;
        id = (queryNeedsToBeRev) ? BIT_CLEAR(static_cast<size_t>(id), 63)
                                 : BIT_SET(static_cast<size_t>(id), 63);
        (kmers + writePos)->kmer = id;
        (kmers + writePos)->id =
            (queryTargetSwitched) ? currQueryKmer->id : currTargetKmer.id;
      } else {
        // i - j
        (kmers + writePos)->kmer =
            (queryTargetSwitched) ? currTargetKmer.id : currQueryKmer->id;
        (kmers + writePos)->id =
            (queryTargetSwitched) ? currQueryKmer->id : currTargetKmer.id;
        //                std::cout << currTargetKmer.pos - currQueryKmer->pos
        //                << "\t" << currTargetKmer.pos << "\t" <<
        //                currQueryKmer->pos << std::endl;
        (kmers + writePos)->pos = (queryTargetSwitched)
                                      ? currTargetKmer.pos - currQueryKmer->pos
                                      : currQueryKmer->pos - currTargetKmer.pos;
      }
      (kmers + writePos)->seqLen = currQueryKmer->seqLen;

      writePos++;
      if (kmerPos + 1 < kmersSize) {
        kmerPos++;
      }
    }
  }
  out->info("Time to find k-mers: {}\n", timer.lap());
  timer.reset();
  if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
    SORT_PARALLEL(kmers, kmers + writePos,
                  KmerPosition<short>::compareRepSequenceAndIdAndDiagReverse);
  } else {
    SORT_PARALLEL(kmers, kmers + writePos,
                  KmerPosition<short>::compareRepSequenceAndIdAndDiag);
  }

  out->info("Time to sort: {}\n", timer.lap());
  return std::make_pair(kmers, writePos);
}

template std::pair<KmerPosition<short> *, size_t> KmerSearch::searchInIndex<0>(
    KmerPosition<short> *kmers, size_t kmersSize, KmerIndex &kmerIndex,
    int resultDirection);
template std::pair<KmerPosition<short> *, size_t> KmerSearch::searchInIndex<1>(
    KmerPosition<short> *kmers, size_t kmersSize, KmerIndex &kmerIndex,
    int resultDirection);

#undef SIZE_T_MAX
