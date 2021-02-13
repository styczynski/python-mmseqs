#include <biosnake/commons/application.h>
#include <biosnake/commons/commandCaller.h>
#include <biosnake/commons/dBReader.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>

#include "createindex.sh.h"
#include <biosnake/output.h>

#include <cassert>
#include <climits>
#include <string>

int createindex(biosnake_output *out, Parameters &par,
                const std::string &indexerModule, const std::string &flag) {
  bool sensitivity = false;
  // only set kmerScore  to INT_MAX if -s was used
  for (size_t i = 0; i < par.createindex.size(); i++) {
    if (par.createindex[i]->uniqid == par.PARAM_S.uniqid &&
        par.createindex[i]->wasSet) {
      par.kmerScore = INT_MAX;
      sensitivity = true;
      break;
    }
  }

  int dbType = FileUtil::parseDbType(out, par.db1.c_str());
  if (Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_HMM_PROFILE) &&
      sensitivity == false) {
    out->error("Please adjust the sensitivity of your target profile index with -s. Be aware that this searches can take huge amount of memory.");
    return EXIT_FAILURE;
  }

  std::string tmpDir = par.db2;
  // TODO: Implement correct hash?
  std::string hash =
      "some_hash__1";  // SSTR(par.hashParameter(par.databases_types,
                       // par.filenames, par.createindex));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();
  par.filenames.push_back(tmpDir);

  out->output_string("INDEXER", indexerModule);
  out->output_string("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : "");
  par.translate = 1;
  out->output_string("ORF_PAR", par.createParameterString(out, par.extractorfs));
  out->output_string("EXTRACT_FRAMES_PAR",
                     par.createParameterString(out, par.extractframes));
  out->output_string("SPLIT_SEQ_PAR",
                     par.createParameterString(out, par.splitsequence));
  if (indexerModule == "kmerindexdb") {
    out->output_string("INDEX_PAR", par.createParameterString(out, par.kmerindexdb));
  } else {
    out->output_string("INDEX_PAR", par.createParameterString(out, par.indexdb));
  }
  if (flag.size() > 0) {
    out->output_string(flag, "1");
  }

  std::string tmp_db_path = "";
  if (flag == "TRANSLATED") {
    tmp_db_path = tmpDir + "/orfs_aa";
    if (!FileUtil::fileExists(out, (tmp_db_path + ".dbtype").c_str())) {
      Parameters extractorfs_par;
      extractorfs_par.setDBFields(1, par.db1);
      extractorfs_par.setDBFields(2, tmp_db_path + ".dbtype");
      extractorfs_par.orfMinLength = 30;
      extractorfs_par.orfMaxLength = 32734;
      extractorfs_par.orfMaxGaps = 2147483647;
      extractorfs_par.contigStartMode = 2;
      extractorfs_par.contigEndMode = 2;
      extractorfs_par.orfStartMode = 1;
      extractorfs_par.forwardFrames = "1,2,3";
      extractorfs_par.reverseFrames = "1,2,3";
      extractorfs_par.translationTable = 1;
      extractorfs_par.translate = 0;
      extractorfs_par.useAllTableStarts = true;
      extractorfs_par.identifierOffset = 0;
      extractorfs_par.createLookup = 0;
      extractorfs_par.threads = 1;
      extractorfs_par.compressed = 0;
      subcall_biosnake(out, "extractorfs", extractorfs_par);
    }
  } else if (flag == "LIN_NUCL") {
    tmp_db_path = tmpDir + "/nucl_split_seq";
    if (!FileUtil::fileExists(out, (tmp_db_path + ".dbtype").c_str())) {
      Parameters extractorfs_par;
      extractorfs_par.maxSeqLen = 65535;
      extractorfs_par.sequenceOverlap = 0;
      extractorfs_par.sequenceSplitMode = 1;
      extractorfs_par.headerSplitMode = 0;
      extractorfs_par.createLookup = 0;
      extractorfs_par.threads = 1;
      extractorfs_par.compressed = 0;
      subcall_biosnake(out, "extractorfs", extractorfs_par);
    }
  }

  Parameters indexer_par;
  indexer_par.setDBFields(1, tmp_db_path);
  indexer_par.setDBFields(2, par.db1);
  indexer_par.setSeedSubstitutionMatrices("VTML80.out", "nucleotide.out");
  indexer_par.kmerSize = 0;
  indexer_par.alphabetSize = MultiParam<int>(21, 5);
  indexer_par.compBiasCorrection = 1;
  indexer_par.maxSeqLen = 65535;
  indexer_par.maxResListLen = 300;
  indexer_par.maskMode = 1;
  indexer_par.maskLowerCaseMode = 0;
  indexer_par.spacedKmer = 1;
  indexer_par.sensitivity = 7.5;
  indexer_par.kmerScore = 0;
  indexer_par.checkCompatible = 0;
  indexer_par.searchType = 3;
  indexer_par.split = 0;
  indexer_par.splitMemoryLimit = 0;
  indexer_par.threads = 1;

  if (flag != "TRANSLATED" && flag != "LIN_NUCL") {
    indexer_par.setDBFields(1, par.db1);
  }

  subcall_biosnake(out, indexerModule, indexer_par);

  if (par.removeTmpFiles) {
    Parameters rmdb_par;
    rmdb_par.setDBFields(1, tmp_db_path);
    subcall_biosnake(out, "rmdb", rmdb_par);
  }

  return 0;
}

int createlinindex(biosnake_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.orfStartMode = 1;
  //    par.orfMinLength = 30;
  //    par.orfMaxLength = 32734;
  //    par.kmerScore = 0; // extract all k-mers
  //    par.maskMode = 0;
  //    par.spacedKmer = false;
  //    // VTML has a slightly lower sensitivity in the regression test
  //    par.seedScoringMatrixFile = MultiParam<char*>("blosum62.out",
  //    "nucleotide.out");
  //
  //    par.PARAM_COV_MODE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_C.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MIN_SEQ_ID.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    for (size_t i = 0; i < par.extractorfs.size(); i++) {
  //        par.extractorfs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.translatenucs.size(); i++) {
  //        par.translatenucs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    par.PARAM_COMPRESSED.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_THREADS.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_V.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  int dbType = FileUtil::parseDbType(out, par.db1.c_str());
  bool isNucl =
      Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);
  if (isNucl && par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES &&
      par.PARAM_MAX_SEQ_LEN.wasSet == false) {
    if (par.PARAM_MAX_SEQ_LEN.wasSet == false) {
      par.maxSeqLen = 10000;
    }
  }
  //// par.printParameters(command.cmd, argc, argv, *command.params);

  if (isNucl && par.searchType == Parameters::SEARCH_TYPE_AUTO) {
    out->warn("Database {} is a nucleotide database. Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)", par.db1);
    return EXIT_FAILURE;
  }
  return createindex(
      out, par, "kmerindexdb",
      (isNucl == false)
          ? ""
          : (par.searchType == Parameters::SEARCH_TYPE_TRANSLATED ||
             par.searchType == Parameters::SEARCH_TYPE_TRANS_NUCL_ALN)
                ? "TRANSLATED"
                : "LIN_NUCL");
}

int createindex(biosnake_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.orfStartMode = 1;
  //    par.orfMinLength = 30;
  //    par.orfMaxLength = 32734;
  //    par.kmerScore = 0; // extract all k-mers
  //    par.sensitivity = 7.5;
  //    par.maskMode = 1;
  //
  //    par.PARAM_COV_MODE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_C.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MIN_SEQ_ID.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MAX_SEQS.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_SPLIT.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    for (size_t i = 0; i < par.splitsequence.size(); i++) {
  //        par.splitsequence[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.extractorfs.size(); i++) {
  //        par.extractorfs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.splitsequence.size(); i++) {
  //        par.splitsequence[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.translatenucs.size(); i++) {
  //        par.translatenucs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    par.PARAM_COMPRESSED.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_THREADS.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_V.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  int dbType = FileUtil::parseDbType(out, par.db1.c_str());
  bool isNucl =
      Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);

  if (par.PARAM_STRAND.wasSet == false) {
    par.strand = 1;
  }
  if (isNucl && par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES) {
    if (par.PARAM_K.wasSet == false) {
      par.kmerSize = 15;
    }
    if (par.PARAM_MAX_SEQ_LEN.wasSet == false) {
      par.maxSeqLen = 10000;
    }

    //  0: reverse, 1: forward, 2: both
    switch (par.strand) {
      case 0:
        par.forwardFrames = "";
        par.reverseFrames = "1";
        break;
      case 1:
        par.forwardFrames = "1";
        par.reverseFrames = "";
        break;
      case 2:
        par.forwardFrames = "1";
        par.reverseFrames = "1";
        break;
    }
  }
  //// par.printParameters(command.cmd, argc, argv, *command.params);
  if (isNucl && par.searchType == Parameters::SEARCH_TYPE_AUTO) {
    out->warn("Database {} is a nucleotide database. Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)", par.db1);
    return EXIT_FAILURE;
  }
  return createindex(
      out, par, "indexdb",
      (isNucl == false)
          ? ""
          : (par.searchType == Parameters::SEARCH_TYPE_TRANSLATED ||
             par.searchType == Parameters::SEARCH_TYPE_TRANS_NUCL_ALN)
                ? "TRANSLATED"
                : "NUCL");
}
