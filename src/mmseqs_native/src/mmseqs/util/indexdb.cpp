#include <mmseqs/commons/dBReader.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/prefiltering/prefiltering.h>
#include <mmseqs/prefiltering/prefilteringIndexReader.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

#ifdef OPENMP
#include <omp.h>
#endif

void setIndexDbDefaults(Parameters *p) { p->sensitivity = 5.7; }

std::string findIncompatibleParameter(DBReader<unsigned int> &index,
                                      const Parameters &par, const int dbtype) {
  PrefilteringIndexData meta = PrefilteringIndexReader::getMetadata(&index);
  if (meta.compBiasCorr != par.compBiasCorrection) return "compBiasCorrection";
  if (meta.maxSeqLength != static_cast<int>(par.maxSeqLen)) return "maxSeqLen";
  if (meta.seqType != dbtype) return "seqType";
  if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_NUCLEOTIDES) ==
          false &&
      par.searchType != Parameters::SEARCH_TYPE_NUCLEOTIDES &&
      meta.alphabetSize != par.alphabetSize.aminoacids)
    return "alphabetSize";
  if (meta.kmerSize != par.kmerSize) return "kmerSize";
  if (meta.mask != par.maskMode) return "maskMode";
  if (meta.kmerThr != par.kmerScore) return "kmerScore";
  if (meta.spacedKmer != par.spacedKmer) return "spacedKmer";
  if (BaseMatrix::unserializeName(par.seedScoringMatrixFile.aminoacids) !=
          PrefilteringIndexReader::getSubstitutionMatrixName(&index) &&
      BaseMatrix::unserializeName(par.seedScoringMatrixFile.nucleotides) !=
          PrefilteringIndexReader::getSubstitutionMatrixName(&index))
    return "seedScoringMatrixFile";
  if (par.spacedKmerPattern !=
      PrefilteringIndexReader::getSpacedPattern(&index))
    return "spacedKmerPattern";
  return "";
}

int indexdb(mmseqs_output *out, Parameters &par) {
  //    Parameters &par = Parameters::getInstance();
  //    setIndexDbDefaults(&par);
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  const bool sameDB = (par.db1 == par.db2);

  std::string path = FileUtil::getRealPathFromSymLink(out, par.db2dbtype);
  DBReader<unsigned int> dbr(
      out, par.db1.c_str(), par.db1Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  dbr.open(DBReader<unsigned int>::NOSORT);

  DBReader<unsigned int> *dbr2 = NULL;
  if (sameDB == false) {
    dbr2 = new DBReader<unsigned int>(
        out, par.db2.c_str(), par.db2Index.c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    dbr2->open(DBReader<unsigned int>::NOSORT);
  }

  const bool db1IsNucl = Parameters::isEqualDbtype(
      dbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
  const bool db2IsNucl =
      dbr2 != NULL && Parameters::isEqualDbtype(dbr2->getDbtype(),
                                                Parameters::DBTYPE_NUCLEOTIDES);
  BaseMatrix *seedSubMat = Prefiltering::getSubstitutionMatrix(
      out, par.seedScoringMatrixFile, par.alphabetSize, 8.0f, false,
      (db1IsNucl && db2IsNucl));

  // memoryLimit in bytes
  size_t memoryLimit = Util::computeMemory(out, par.splitMemoryLimit);

  int splitMode = Parameters::TARGET_DB_SPLIT;
  par.maxResListLen = std::min(dbr.getSize(), par.maxResListLen);
  Prefiltering::setupSplit(
      out, dbr, seedSubMat->alphabetSize - 1, dbr.getDbtype(), par.threads, false,
      memoryLimit, 1, par.maxResListLen, par.kmerSize, par.split, splitMode);

  bool kScoreSet = false;
  for (size_t i = 0; i < par.indexdb.size(); i++) {
    if (par.indexdb[i]->uniqid == par.PARAM_K_SCORE.uniqid &&
        par.indexdb[i]->wasSet) {
      kScoreSet = true;
    }
  }

  const bool isProfileSearch = (Parameters::isEqualDbtype(
      dbr.getDbtype(), Parameters::DBTYPE_HMM_PROFILE));
  if (isProfileSearch && kScoreSet == false) {
    par.kmerScore = 0;
  }

  // TODO: investigate if it makes sense to mask the profile consensus sequence
  if (isProfileSearch) {
    par.maskMode = 0;
  }

  // query seq type is actually unknown here, but if we pass DBTYPE_HMM_PROFILE
  // then its +20 k-score
  par.kmerScore = Prefiltering::getKmerThreshold(
      out, par.sensitivity, isProfileSearch, par.kmerScore, par.kmerSize);

  std::string indexDB = PrefilteringIndexReader::indexName(par.db2);

  int status = EXIT_SUCCESS;
  bool recreate = true;
  std::string indexDbType = indexDB + ".dbtype";
  if (par.checkCompatible > 0 && FileUtil::fileExists(out, indexDbType.c_str())) {
    out->info("Check index {}\n", indexDB);
    DBReader<unsigned int> index(
        out, indexDB.c_str(), (indexDB + ".index").c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    index.open(DBReader<unsigned int>::NOSORT);

    if (Parameters::isEqualDbtype(dbr.getDbtype(),
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES &&
        par.PARAM_ALPH_SIZE.wasSet) {
      out->warn("Alphabet size is not taken into account for compatibility check in nucleotide search");
    }

    std::string check;
    const bool compatible =
        PrefilteringIndexReader::checkIfIndexFile(&index) &&
        (check = findIncompatibleParameter(index, par, dbr.getDbtype())) == "";
    index.close();
    if (compatible) {
      out->info("Index is up to date and compatible. Force recreation with --check-compatibility 0 parameter");
      recreate = false;
    } else {
      if (par.checkCompatible == 2) {
        out->error("Index is incompatible. Incompatible parameter: {}", check);
        recreate = false;
        status = EXIT_FAILURE;
      } else {
        out->warn("Index is incompatible and will be recreated. Incompatible parameter: {}", check);
        recreate = true;
      }
    }
  }

  if (recreate) {
    DBReader<unsigned int> hdbr1(
        out, par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    hdbr1.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *hdbr2 = NULL;
    if (sameDB == false) {
      hdbr2 = new DBReader<unsigned int>(
          out, par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads,
          DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
      hdbr2->open(DBReader<unsigned int>::NOSORT);
    }

    DBReader<unsigned int>::removeDb(out, indexDB);
    PrefilteringIndexReader::createIndexFile(
        out, indexDB, &dbr, dbr2, &hdbr1, hdbr2, seedSubMat, par.maxSeqLen,
        par.spacedKmer, par.spacedKmerPattern, par.compBiasCorrection,
        seedSubMat->alphabetSize, par.kmerSize, par.maskMode,
        par.maskLowerCaseMode, par.kmerScore, par.split);

    if (hdbr2 != NULL) {
      hdbr2->close();
      delete hdbr2;
    }

    hdbr1.close();
  }

  if (dbr2 != NULL) {
    dbr2->close();
    delete dbr2;
  }

  delete seedSubMat;
  dbr.close();

  return status;
}
