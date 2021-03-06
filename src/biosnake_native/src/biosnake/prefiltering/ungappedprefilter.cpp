//
// Created by Martin Steinegger on 17.09.18.
//

#include <biosnake/commons/dBReader.h>
#include <biosnake/commons/dBWriter.h>
#include <biosnake/output.h>
#include <biosnake/alignment/distanceCalculator.h>
#include <biosnake/commons/fastSort.h>
#include <biosnake/alignment/matcher.h>
#include <biosnake/commons/nucleotideMatrix.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/prefiltering/queryMatcher.h>
#include <biosnake/commons/substitutionMatrixProfileStates.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#ifdef OPENMP
#include <omp.h>
#endif

int doRescorealldiagonal(biosnake_output* out, Parameters &par, DBReader<unsigned int> &qdbr,
                         DBWriter &resultWriter, size_t dbStart,
                         size_t dbSize) {
  int querySeqType = qdbr.getDbtype();
  DBReader<unsigned int> *tdbr = NULL;
  bool sameDB = false;
  if (par.db1.compare(par.db2) == 0) {
    sameDB = true;
    tdbr = &qdbr;
  } else {
    tdbr = new DBReader<unsigned int>(
        out,
        par.db2.c_str(), par.db2Index.c_str(), par.threads,
        DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    tdbr->open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
      tdbr->readMmapedDataInMemory();
    }
  }
  const int targetSeqType = tdbr->getDbtype();
  if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE) &&
      Parameters::isEqualDbtype(targetSeqType,
                                Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
    querySeqType = Parameters::DBTYPE_PROFILE_STATE_PROFILE;
  }

  BaseMatrix *subMat;
  EvalueComputation *evaluer;
  int8_t *tinySubMat;
  if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
    subMat = new NucleotideMatrix(out, par.scoringMatrixFile.nucleotides, 1.0, 0.0);
    evaluer = new EvalueComputation(out, tdbr->getAminoAcidDBSize(), subMat);
    tinySubMat = new int8_t[subMat->alphabetSize * subMat->alphabetSize];
    for (int i = 0; i < subMat->alphabetSize; i++) {
      for (int j = 0; j < subMat->alphabetSize; j++) {
        tinySubMat[i * subMat->alphabetSize + j] = subMat->subMatrix[i][j];
      }
    }
  } else if (Parameters::isEqualDbtype(targetSeqType,
                                       Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
    SubstitutionMatrix sMat(out, par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    evaluer = new EvalueComputation(out, tdbr->getAminoAcidDBSize(), &sMat);
    subMat = new SubstitutionMatrixProfileStates(
        out,
        sMat.matrixName, sMat.probMatrix, sMat.pBack,
        sMat.subMatrixPseudoCounts, 2.0, 0.0, 219);
    tinySubMat = new int8_t[sMat.alphabetSize * sMat.alphabetSize];
    for (int i = 0; i < sMat.alphabetSize; i++) {
      for (int j = 0; j < sMat.alphabetSize; j++) {
        tinySubMat[i * sMat.alphabetSize + j] = sMat.subMatrix[i][j];
      }
    }
  } else {
    // keep score bias at 0.0 (improved ROC)
    subMat = new SubstitutionMatrix(out, par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    evaluer = new EvalueComputation(out, tdbr->getAminoAcidDBSize(), subMat);
    tinySubMat = new int8_t[subMat->alphabetSize * subMat->alphabetSize];
    for (int i = 0; i < subMat->alphabetSize; i++) {
      for (int j = 0; j < subMat->alphabetSize; j++) {
        tinySubMat[i * subMat->alphabetSize + j] = subMat->subMatrix[i][j];
      }
    }
  }

  Log::Progress progress(dbSize);

#pragma omp parallel
  {
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = (unsigned int)omp_get_thread_num();
#endif
    char buffer[1024 + 32768];
    std::vector<hit_t> shortResults;
    shortResults.reserve(std::max(static_cast<size_t>(1), tdbr->getSize() / 5));
    Sequence qSeq(out, par.maxSeqLen, querySeqType, subMat, 0, false,
                  par.compBiasCorrection);
    Sequence tSeq(out, par.maxSeqLen, targetSeqType, subMat, 0, false,
                  par.compBiasCorrection);
    SmithWaterman aligner(out, par.maxSeqLen, subMat->alphabetSize,
                          par.compBiasCorrection);

    std::string resultBuffer;
    resultBuffer.reserve(262144);
#pragma omp for schedule(dynamic, 1)
    for (size_t id = dbStart; id < (dbStart + dbSize); id++) {
      progress.updateProgress();
      char *querySeqData = qdbr.getData(id, thread_idx);
      size_t queryKey = qdbr.getDbKey(id);
      unsigned int querySeqLen = qdbr.getSeqLen(id);

      qSeq.mapSequence(id, queryKey, querySeqData, querySeqLen);
      //            qSeq.printProfileStatePSSM();
      if (Parameters::isEqualDbtype(qSeq.getSeqType(),
                                    Parameters::DBTYPE_HMM_PROFILE) ||
          Parameters::isEqualDbtype(qSeq.getSeqType(),
                                    Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
        aligner.ssw_init(&qSeq, qSeq.getAlignmentProfile(), subMat, 0);
      } else {
        aligner.ssw_init(&qSeq, tinySubMat, subMat, 0);
      }

      for (size_t tId = 0; tId < tdbr->getSize(); tId++) {
        unsigned int targetKey = tdbr->getDbKey(tId);
        const bool isIdentity =
            (queryKey == targetKey && (par.includeIdentity || sameDB)) ? true
                                                                       : false;
        char *targetSeq = tdbr->getData(tId, thread_idx);
        unsigned int targetSeqLen = tdbr->getSeqLen(tId);
        tSeq.mapSequence(tId, targetKey, targetSeq, targetSeqLen);
        //                tSeq.print();
        float queryLength = qSeq.L;
        float targetLength = tSeq.L;
        if (Util::canBeCovered(par.covThr, par.covMode, queryLength,
                               targetLength) == false) {
          continue;
        }

        int score = aligner.ungapped_alignment(tSeq.numSequence, tSeq.L);
        bool hasDiagScore = (score > par.minDiagScoreThr);
        double evalue = evaluer->computeEvalue(score, qSeq.L);
        bool hasEvalue = (evalue <= par.evalThr);
        // --filter-hits
        if (isIdentity || (hasDiagScore && hasEvalue)) {
          hit_t hit;
          hit.seqId = targetKey;
          hit.prefScore = score;
          hit.diagonal = 0;
          shortResults.emplace_back(hit);
        }
      }

      SORT_SERIAL(shortResults.begin(), shortResults.end(),
                  hit_t::compareHitsByScoreAndId);
      size_t maxSeqs = std::min(par.maxResListLen, shortResults.size());
      for (size_t i = 0; i < maxSeqs; ++i) {
        size_t len =
            QueryMatcher::prefilterHitToBuffer(buffer, shortResults[i]);
        resultBuffer.append(buffer, len);
      }

      resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(),
                             queryKey, thread_idx);
      resultBuffer.clear();
      shortResults.clear();
    }
  }

  qdbr.close();
  if (sameDB == false) {
    tdbr->close();
    delete tdbr;
  }

  delete[] tinySubMat;
  delete subMat;
  delete evaluer;
  return 0;
}

int ungappedprefilter(biosnake_output *out, Parameters &par) {
  //    BiosnakeMPI::init(argc, argv);
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  DBReader<unsigned int> qdbr(
      out,
      par.db1.c_str(), par.db1Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
  qdbr.open(DBReader<unsigned int>::NOSORT);
  if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
    qdbr.readMmapedDataInMemory();
  }

#ifdef HAVE_MPI
  size_t dbFrom = 0;
  size_t dbSize = 0;

  qdbr.decomposeDomainByAminoAcid(BiosnakeMPI::rank, BiosnakeMPI::numProc, &dbFrom,
                                  &dbSize);
  std::pair<std::string, std::string> tmpOutput =
      Util::createTmpFileNames(par.db3, par.db3Index, BiosnakeMPI::rank);
  DBWriter resultWriter(out, tmpOutput.first.c_str(), tmpOutput.second.c_str(),
                        par.threads, par.compressed,
                        Parameters::DBTYPE_PREFILTER_RES);
  resultWriter.open();
  int status = doRescorealldiagonal(par, qdbr, resultWriter, dbFrom, dbSize);
  resultWriter.close();

  MPI_Barrier(MPI_COMM_WORLD);
  if (BiosnakeMPI::rank == 0) {
    std::vector<std::pair<std::string, std::string>> splitFiles;
    for (int proc = 0; proc < BiosnakeMPI::numProc; ++proc) {
      std::pair<std::string, std::string> tmpFile =
          Util::createTmpFileNames(par.db3, par.db3Index, proc);
      splitFiles.push_back(std::make_pair(tmpFile.first, tmpFile.second));
    }
    DBWriter::mergeResults(out, par.db3, par.db3Index, splitFiles);
  }
#else
  DBWriter resultWriter(out, par.db3.c_str(), par.db3Index.c_str(), par.threads,
                        par.compressed, Parameters::DBTYPE_PREFILTER_RES);
  resultWriter.open();
  int status = doRescorealldiagonal(out, par, qdbr, resultWriter, 0, qdbr.getSize());
  resultWriter.close();
#endif
  return status;
}
