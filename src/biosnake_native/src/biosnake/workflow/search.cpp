#include <biosnake/commons/application.h>
#include <biosnake/commons/commandCaller.h>
#include <biosnake/commons/dBReader.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/prefiltering/prefilteringIndexReader.h>
#include <biosnake/commons/util.h>
#include "blastn.sh.h"
#include "blastp.sh.h"
#include "blastpgp.sh.h"
#include <biosnake/output.h>
#include "searchslicedtargetprofile.sh.h"
#include "searchtargetprofile.sh.h"
#include "translated_search.sh.h"

#include <cassert>
#include <climits>
#include <iomanip>
#include <iostream>
#include <set>

/*
std::string blastp_input = query;
std::string blastp_target = target;
std::string blastp_tmp = searchTmpDir;
std::string blastp_out = tmpDir + "/aln";
*/
void call_blastp(biosnake_output *out, Parameters &par, int no_steps,
                 std::vector<float> senses, std::string align_module,
                 std::string blastp_input, std::string blastp_target,
                 std::string blastp_out, std::string blastp_tmp,
                 std::string script_path, CommandCaller cmd) {
  if (false) {
    out->info("[blastn.sh] Executing native script: {}", script_path);
    std::vector<std::string> args = {blastp_input, blastp_target, blastp_out,
                                     blastp_tmp};
    cmd.addVar("BIOSNAKE", "biosnake");
    cmd.callProgram(script_path.c_str(), args);
    out->info("[blastn.sh] Exit");
  } else {
    out->info("step_search L");

    // "$SEARCH" "${QUERY}" "${TARGET}" "$4/aln" "$4/search"
    int blastp_steps = no_steps;
    int current_step = 0;

    std::set<std::string> merged_senses;

    std::string aln_res_merge = blastp_tmp + "/aln_0";
    while (current_step < blastp_steps) {
      out->info("step_search K_1");
      const float sens = senses[current_step];
      std::string sens_str = std::to_string(((float)((int)(sens * 10))) / 10);
      // Call prefilter module
      std::string step_str = std::to_string(current_step);
      if (!FileUtil::fileExists(
              out, (blastp_tmp + "/pref_" + step_str + ".dbtype").c_str())) {
        out->info("step_search K_2: prefilter {} {} {}", blastp_input, blastp_target, (blastp_tmp + "/pref_" + step_str));
        Parameters prefilter_par(par);
        std::vector<std::string> prefilter_filenames = {
            blastp_input, blastp_target, (blastp_tmp + "/pref_" + step_str)};
        prefilter_par.filenames = prefilter_filenames;
        prefilter_par.setDBFields(1, blastp_input);
        prefilter_par.setDBFields(2, blastp_target);
        prefilter_par.setDBFields(3, blastp_tmp + "/pref_" + step_str);
        prefilter_par.setSubstitutionMatrices("blosum62.out", "nucleotide.out");
        prefilter_par.setSeedSubstitutionMatrices("VTML80.out",
                                                  "nucleotide.out");
        prefilter_par.kmerSize = 0;
        prefilter_par.kmerScore = 2147483647;
        prefilter_par.alphabetSize = MultiParam<int>(21, 5);
        prefilter_par.split = 0;
        prefilter_par.splitMode = 2;
        prefilter_par.splitMemoryLimit = 0;
        prefilter_par.covThr = 0;
        prefilter_par.covMode = 0;
        prefilter_par.compBiasCorrection = 1;
        prefilter_par.diagonalScoring = 1;
        prefilter_par.exactKmerMatching = 0;
        prefilter_par.maskMode = 1;
        prefilter_par.maskLowerCaseMode = 0;
        prefilter_par.minDiagScoreThr = 15;
        prefilter_par.includeIdentity = 0;
        prefilter_par.spacedKmer = 1;
        prefilter_par.preloadMode = 0;
        prefilter_par.pca = 1;
        prefilter_par.pcb = 1.5;
        prefilter_par.compressed = 0;
        prefilter_par.sensitivity = sens;
        subcall_biosnake(out, "prefilter", prefilter_par);
        out->info("step_search K_3");
      }

      // call alignment module
      std::string align_path = blastp_out;
      if (blastp_steps > 1) {
        align_path = blastp_tmp + "/aln_" + step_str;
      }
      out->info("step_search K_4");
      if (!FileUtil::fileExists(out, (align_path + ".dbtype").c_str())) {
        out->info("step_search K_5");
        Parameters align_par(par);
        std::vector<std::string> align_filenames = {
            blastp_input,
            blastp_target,
            (blastp_tmp + "/pref_" + step_str),
            align_path,
        };
        align_par.setDBFields(1, blastp_input);
        align_par.setDBFields(2, blastp_target);
        align_par.setDBFields(3, blastp_tmp + "/pref_" + step_str);
        align_par.setDBFields(4, align_path);
        align_par.filenames = align_filenames;
        align_par.compressed = 0;
        align_par.setSubstitutionMatrices("blosum62.out", "nucleotide.out");
        align_par.addBacktrace = 0;
        align_par.alignmentMode = 3;
        align_par.wrappedScoring = 0;
        align_par.evalThr = 0.001;
        align_par.seqIdThr = 0;
        align_par.alnLenThr = 0;
        align_par.seqIdMode = 0;
        align_par.altAlignment = 0;
        align_par.covThr = 0;
        align_par.covMode = 0;
        align_par.compBiasCorrection = 1;
        align_par.maxRejected = 2147483647;
        align_par.maxAccept = 2147483647;
        align_par.includeIdentity = 0;
        align_par.preloadMode = 0;
        align_par.pca = 1;
        align_par.pcb = 1.5;
        align_par.scoreBias = 0;
        align_par.realign = 0;
        align_par.realignScoreBias = -0.2;
        align_par.realignMaxSeqs = 2147483647;
        align_par.gapOpen = MultiParam<int>(11, 5);
        align_par.gapExtend = MultiParam<int>(1, 2);
        align_par.zdrop = 40;
        align_par.addBacktrace = true;

        subcall_biosnake(out, align_module, align_par);
        out->info("step_search K_6");
      }

      // Only merge results after first step
      out->info("step_search K_7");
      if (current_step > 0) {
        if (merged_senses.find(sens_str) == merged_senses.end()) {
          if (current_step < blastp_steps - 1) {
            out->info("step_search K_8");
            Parameters mergedbs_par(par);
            std::vector<std::string> mergedbs_filenames = {
                blastp_input,
                (blastp_tmp + "/aln_merge_new"),
                (blastp_tmp + "/pref_" + step_str),
                align_path,
            };
            mergedbs_par.filenames = mergedbs_filenames;
            mergedbs_par.setDBFields(1, blastp_input);
            mergedbs_par.setDBFields(2, blastp_tmp + "/aln_merge_new");
            mergedbs_par.setDBFields(3, blastp_tmp + "/pref_" + step_str);
            mergedbs_par.setDBFields(4, align_path);
            mergedbs_par.compressed = 0;
            subcall_biosnake(out, "mergedbs", mergedbs_par);

            Parameters rmdb_par(par);
            rmdb_par.setDBFields(1, blastp_tmp + "/aln_merge");
            subcall_biosnake(out, "rmdb", rmdb_par);

            Parameters mvdb_par(par);
            mvdb_par.setDBFields(1, blastp_tmp + "/aln_merge_new");
            mvdb_par.setDBFields(2, blastp_tmp + "/aln_merge");
            subcall_biosnake(out, "mvdb", mvdb_par);
            out->info("step_search K_9");
          } else {
            out->info("step_search K_10");
            Parameters mergedbs_par(par);
            std::vector<std::string> mergedbs_filenames = {
                blastp_input, align_path, (blastp_tmp + "/pref_" + step_str),
                align_path};
            mergedbs_par.filenames = mergedbs_filenames;
            mergedbs_par.setDBFields(1, blastp_input);
            mergedbs_par.setDBFields(2, align_path);
            mergedbs_par.setDBFields(3, blastp_tmp + "/pref_" + step_str);
            mergedbs_par.setDBFields(4, align_path);
            mergedbs_par.compressed = 0;
            subcall_biosnake(out, "mergedbs", mergedbs_par);
            out->info("step_search K_11");
          }
          merged_senses.insert(sens_str);
        }
        aln_res_merge = blastp_tmp + "/aln_merge";
      }
      out->info("step_search K_12");

      std::string next_input = blastp_tmp + "/input_" + step_str;
      // do not create subdb at last step
      if (current_step < blastp_steps - 1) {
        out->info("step_search K_13");
        if (!FileUtil::fileExists(
                out, (blastp_tmp + "/order_" + step_str + ".dbtype").c_str())) {
          // awk '$3 < 2 { print $1 }' "$TMP_PATH/aln_$STEP.index" >
          // "$TMP_PATH/order_$STEP"
          // TODO: Implement this blastp code
          out->failure("Reached code branch that is yet uninmplemented (TODO: Implement this blastp code)");
        }
        out->info("step_search K_14");

        if (!FileUtil::fileExistsAndIsNotEmpty(
                out, (blastp_tmp + "/order_" + step_str).c_str())) {
          // "$BIOSNAKE" mvdb "$ALN_RES_MERGE" "$3" ${VERBOSITY}
          Parameters mvdb_par(par);
          mvdb_par.setDBFields(1, aln_res_merge);
          mvdb_par.setDBFields(2, align_path);
          subcall_biosnake(out, "mvdb", mvdb_par);
        }
        out->info("step_search K_15");

        if (!FileUtil::fileExists(out, (next_input + ".dbtype").c_str())) {
          Parameters createsubdb_par(par);
          std::vector<std::string> createsubdb_filenames = {
              (blastp_tmp + "/order_" + step_str),
              blastp_input,
              next_input,
          };
          createsubdb_par.filenames = createsubdb_filenames;
          createsubdb_par.setDBFields(1, blastp_tmp + "/order_" + step_str);
          createsubdb_par.setDBFields(2, blastp_input);
          createsubdb_par.setDBFields(3, next_input);
          createsubdb_par.subDbMode = 1;
          subcall_biosnake(out, "createsubdb", createsubdb_par);
        }
        out->info("step_search K_16");
      }

      blastp_input = next_input;
      current_step++;
      out->info("step_search K_17");
    }
    out->info("step_search K_18");

    if (par.removeTmpFiles) {
      for (int i = 0; i < blastp_steps; ++i) {
        std::string step_str = std::to_string(i);
        Parameters rmdb_par(par);
        rmdb_par.setDBFields(1, blastp_tmp + "/pref_" + step_str);
        subcall_biosnake(out, "rmdb", rmdb_par);

        rmdb_par.setDBFields(1, blastp_tmp + "/aln_" + step_str);
        subcall_biosnake(out, "rmdb", rmdb_par);

        rmdb_par.setDBFields(1, blastp_tmp + "/input_" + step_str);
        subcall_biosnake(out, "rmdb", rmdb_par);
      }

      Parameters rmdb_par(par);
      rmdb_par.setDBFields(1, blastp_tmp + "/aln_merge");
      subcall_biosnake(out, "rmdb", rmdb_par);
    }
    out->info("step_search K_19");
  }
}

void setSearchDefaults(Parameters *p) {
  p->spacedKmer = true;
  p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
  p->sensitivity = 5.7;
  p->evalThr = 0.001;
  p->orfStartMode = 1;
  p->orfMinLength = 30;
  p->orfMaxLength = 32734;
  p->evalProfile = 0.1;
}

int computeSearchMode(biosnake_output* out, int queryDbType, int targetDbType, int targetSrcDbType,
                      int searchType) {
  // reject unvalid search
  if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
      Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)) {
    out->failure("Profile-Profile searches are not supported");
  }
  // index was used
  if (targetSrcDbType == -1) {
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES)) {
      if (searchType == Parameters::SEARCH_TYPE_AUTO) {
        // WARNING because its not really an error, just a req. parameter
        out->failure("It is unclear from the input if a translated or nucleotide search should be performed. Please provide the parameter --search-type 2 (translated), 3 (nucleotide) or 4 (translated nucleotide backtrace)");
      }
      // nucl/nucl
      // nucl/nucl translated
      if (searchType == Parameters::SEARCH_TYPE_TRANSLATED ||
          searchType == Parameters::SEARCH_TYPE_TRANS_NUCL_ALN) {
        return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED |
               Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
      } else if (searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES) {
        return Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE |
               Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE;
      } else {
        out->failure("--search-type 1 (amino acid) can not used in combination with a nucleotide database. The only possible options --search-types 2 (translated), 3 (nucleotide) or 4 (translated nucleotide backtrace)");
      }
    }
    // protein/protein
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID |
             Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
    }

    // profile/protein
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE |
             Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
    }

    // protein/profile
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_HMM_PROFILE)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID |
             Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
    }

    // profile/nucleotide
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE |
             Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
    }

    // nucleotide/profile
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_HMM_PROFILE)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED |
             Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
    }

    // nucleotide/protein
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED |
             Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
    }

    // protein/nucleotide
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID |
             Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
    }
  } else {
    // protein/protein
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetSrcDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID |
             Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
    }

    // profile/protein
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetSrcDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE |
             Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
    }

    // protein/profile
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetSrcDbType,
                                  Parameters::DBTYPE_HMM_PROFILE)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID |
             Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
    }

    // profile/nucleotide
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetSrcDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE |
             Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
    }

    // nucleotide/profile
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetSrcDbType,
                                  Parameters::DBTYPE_HMM_PROFILE)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED |
             Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
    }

    // nucleotide/protein
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        (Parameters::isEqualDbtype(targetSrcDbType,
                                   Parameters::DBTYPE_AMINO_ACIDS))) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED |
             Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
    }

    // protein/nucleotide
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        (Parameters::isEqualDbtype(targetSrcDbType,
                                   Parameters::DBTYPE_NUCLEOTIDES))) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID |
             Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
    }

    // nucl/nucl
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        Parameters::isEqualDbtype(targetSrcDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE |
             Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE;
    }

    // nucl/nucl translated
    if (Parameters::isEqualDbtype(queryDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES) &&
        Parameters::isEqualDbtype(targetDbType,
                                  Parameters::DBTYPE_AMINO_ACIDS) &&
        Parameters::isEqualDbtype(targetSrcDbType,
                                  Parameters::DBTYPE_NUCLEOTIDES)) {
      return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED |
             Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
    }
  }
  Parameters::getDbTypeName(queryDbType),
  out->failure("Invalid input database and --search-type combination. queryDbType: {}, targetSrcDbType:: {}, searchMode: {}", Parameters::getDbTypeName(targetDbType), Parameters::getDbTypeName(targetSrcDbType), searchType);
}

void setNuclSearchDefaults(Parameters *p) {
  // leave ungapped alignment untouched
  if (p->alignmentMode != Parameters::ALIGNMENT_MODE_UNGAPPED) {
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
  }
  // p->orfLongest = true;
  p->exactKmerMatching = true;
  //    if ( p->PARAM_DIAGONAL_SCORING.wasSet == false) {
  //        p->diagonalScoring = 0;
  //    }
  if (p->PARAM_STRAND.wasSet == false) {
    p->strand = 2;
  }
  if (p->PARAM_K.wasSet == false) {
    p->kmerSize = 15;
  }
  if (p->PARAM_MAX_SEQ_LEN.wasSet == false) {
    p->maxSeqLen = 10000;
  }
}

int search(biosnake_output *out, Parameters &par) {
  //    Parameters &par = Parameters::getInstance();
  //    setSearchDefaults(&par);
  //    par.PARAM_COV_MODE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_C.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MIN_SEQ_ID.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    for (size_t i = 0; i < par.extractorfs.size(); i++) {
  //        par.extractorfs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.translatenucs.size(); i++) {
  //        par.translatenucs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.splitsequence.size(); i++) {
  //        par.splitsequence[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    par.PARAM_COMPRESSED.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_THREADS.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_V.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //
  //    par.parseParameters(argc, argv, command, false, 0,
  //    BiosnakeParameter::COMMAND_ALIGN | BiosnakeParameter::COMMAND_PREFILTER);

  out->info("Search: {} {}", par.db1, par.db2);

  std::string indexStr = PrefilteringIndexReader::searchForIndex(out, par.db2);

  int targetDbType = FileUtil::parseDbType(out, par.db2.c_str());
  std::string targetDB = (indexStr == "") ? par.db2.c_str() : indexStr.c_str();
  int targetSrcDbType = -1;
  if (indexStr != "" ||
      Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_INDEX_DB)) {
    indexStr = par.db2;
    DBReader<unsigned int> dbr(
        out, targetDB.c_str(), (targetDB + ".index").c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::NOSORT);
    PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&dbr);
    targetSrcDbType = data.srcSeqType;
    targetDbType = data.seqType;
  }
  const int queryDbType = FileUtil::parseDbType(out, par.db1.c_str());
  if (queryDbType == -1 || targetDbType == -1) {
    out->failure("Please recreate your database or add a .dbtype file to your sequence/profile database");
  }

  int searchMode = computeSearchMode(out, queryDbType, targetDbType, targetSrcDbType,
                                     par.searchType);
  if ((searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE) &&
      (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE)) {
    setNuclSearchDefaults(&par);
  } else {
    par.PARAM_STRAND.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  }
  // FIXME: use larger default k-mer size in target-profile case if memory is
  // available overwrite default kmerSize for target-profile searches and parse
  // parameters again
  if (par.sliceSearch == false &&
      (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) &&
      par.PARAM_K.wasSet == false) {
    par.kmerSize = 5;
  }

  const bool isUngappedMode =
      par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
  if (isUngappedMode &&
      (searchMode & (Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE |
                     Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE))) {
    // par.printUsageMessage(command, BiosnakeParameter::COMMAND_ALIGN |
    // BiosnakeParameter::COMMAND_PREFILTER);
    out->failure("Cannot use ungapped alignment mode with profile databases");
  }

  if (isUngappedMode && par.lcaSearch) {
    // par.printUsageMessage(command, BiosnakeParameter::COMMAND_ALIGN |
    // BiosnakeParameter::COMMAND_PREFILTER);
    out->failure("Cannot use ungapped alignment mode with lca search");
  }

  // validate and set parameters for iterative search
  if (par.numIterations > 1) {
    if (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) {
      // par.printUsageMessage(command, BiosnakeParameter::COMMAND_ALIGN |
      // BiosnakeParameter::COMMAND_PREFILTER);
      out->failure("Iterative target-profile searches are not supported");
    }

    par.addBacktrace = true;
    if (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE) {
      for (size_t i = 0; i < par.searchworkflow.size(); i++) {
        if (par.searchworkflow[i]->uniqid == par.PARAM_REALIGN.uniqid &&
            par.searchworkflow[i]->wasSet) {
          // par.printUsageMessage(command,
          //                      BiosnakeParameter::COMMAND_ALIGN |
          //                      BiosnakeParameter::COMMAND_PREFILTER);
          out->failure("Cannot realign query profiles");
        }
      }

      par.realign = false;
    }
  }
  // par.printParameters(command.cmd, argc, argv, par.searchworkflow);

  std::string tmpDir = par.db4;
  std::string hash = SSTR(par.hashParameter(out, par.databases_types, par.filenames,
                                            par.searchworkflow));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  std::string parentTmpDir = tmpDir;
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();
  par.filenames.push_back(tmpDir);

  const int originalRescoreMode = par.rescoreMode;
  bool needTargetSplit = false;
  bool extract_frames = false;
  bool need_query_split = false;
  int no_steps = 0;
  std::string align_module = "";
  std::vector<float> senses;
  std::string search_prog = "";

  CommandCaller cmd(out);
  cmd.addVar("VERBOSITY", par.createParameterString(out, par.onlyverbosity));
  cmd.addVar("THREADS_COMP_PAR",
             par.createParameterString(out, par.threadsandcompression));
  cmd.addVar("VERB_COMP_PAR",
             par.createParameterString(out, par.verbandcompression));
  if (isUngappedMode) {
    cmd.addVar("ALIGN_MODULE", "rescorediagonal");
    align_module = "rescorediagonal";
  } else if (par.lcaSearch) {
    cmd.addVar("ALIGN_MODULE", "lcaalign");
    align_module = "lcaalign";
  } else {
    cmd.addVar("ALIGN_MODULE", "align");
    align_module = "align";
  }
  cmd.addVar("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : "");
  std::string program;
  cmd.addVar("RUNNER", par.runner);
  //    cmd.addVar("ALIGNMENT_DB_EXT", Parameters::isEqualDbtype(targetDbType,
  //    Parameters::DBTYPE_PROFILE_STATE_SEQ) ? ".255" : "");
  par.filenames[1] = targetDB;
  if (par.sliceSearch == true) {
    // By default (0), diskSpaceLimit (in bytes) will be set in the workflow to
    // use as much as possible
    cmd.addVar("AVAIL_DISK", SSTR(static_cast<size_t>(par.diskSpaceLimit)));

    // correct Eval threshold for inverted search
    const size_t queryDbSize = FileUtil::countLines(out, par.db1Index.c_str());
    const size_t targetDbSize = FileUtil::countLines(out, par.db2Index.c_str());
    par.evalThr *= ((float)queryDbSize) / targetDbSize;

    int originalCovMode = par.covMode;
    par.covMode = Util::swapCoverageMode(out, par.covMode);
    size_t maxResListLen = par.maxResListLen;
    par.maxResListLen = INT_MAX;
    cmd.addVar("PREFILTER_PAR", par.createParameterString(out, par.prefilter));
    par.maxResListLen = maxResListLen;
    double originalEvalThr = par.evalThr;
    par.evalThr = std::numeric_limits<double>::max();
    cmd.addVar("SWAP_PAR", par.createParameterString(out, par.swapresult));
    par.evalThr = originalEvalThr;
    if (isUngappedMode) {
      par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
      cmd.addVar("ALIGNMENT_PAR",
                 par.createParameterString(out, par.rescorediagonal));
      par.rescoreMode = originalRescoreMode;
    } else {
      cmd.addVar("ALIGNMENT_PAR", par.createParameterString(out, par.align));
    }
    cmd.addVar("SORTRESULT_PAR", par.createParameterString(out, par.sortresult));
    par.covMode = originalCovMode;

    program = tmpDir + "/searchslicedtargetprofile.sh";
    FileUtil::writeFile(out, program, searchslicedtargetprofile_sh,
                        searchslicedtargetprofile_sh_len);
  } else if (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) {
    cmd.addVar("PREFILTER_PAR", par.createParameterString(out, par.prefilter));
    // we need to align all hits in case of target Profile hits
    size_t maxResListLen = par.maxResListLen;
    par.maxResListLen = INT_MAX;
    int originalCovMode = par.covMode;
    par.covMode = Util::swapCoverageMode(out, par.covMode);
    if (isUngappedMode) {
      par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
      cmd.addVar("ALIGNMENT_PAR",
                 par.createParameterString(out, par.rescorediagonal));
      par.rescoreMode = originalRescoreMode;
    } else {
      cmd.addVar("ALIGNMENT_PAR", par.createParameterString(out, par.align));
    }
    par.covMode = originalCovMode;
    par.maxResListLen = maxResListLen;
    cmd.addVar("SWAP_PAR", par.createParameterString(out, par.swapresult));
    FileUtil::writeFile(out, tmpDir + "/searchtargetprofile.sh",
                        searchtargetprofile_sh, searchtargetprofile_sh_len);
    program = std::string(tmpDir + "/searchtargetprofile.sh");
  } else if (par.numIterations > 1) {
    cmd.addVar("NUM_IT", SSTR(par.numIterations));
    cmd.addVar("SUBSTRACT_PAR", par.createParameterString(out, par.subtractdbs));
    cmd.addVar("VERBOSITY_PAR", par.createParameterString(out, par.onlyverbosity));

    double originalEval = par.evalThr;
    par.evalThr =
        (par.evalThr < par.evalProfile) ? par.evalThr : par.evalProfile;
    for (int i = 0; i < par.numIterations; i++) {
      if (i == 0 &&
          (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) == false) {
        par.realign = true;
      }

      if (i > 0) {
        //                par.queryProfile = true;
        par.realign = false;
      }

      if (i == (par.numIterations - 1)) {
        par.evalThr = originalEval;
      }

      cmd.addVar(std::string("PREFILTER_PAR_" + SSTR(i)),
                 par.createParameterString(out, par.prefilter));
      if (isUngappedMode) {
        par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
        cmd.addVar(std::string("ALIGNMENT_PAR_" + SSTR(i)),
                   par.createParameterString(out, par.rescorediagonal));
        par.rescoreMode = originalRescoreMode;
      } else {
        cmd.addVar(std::string("ALIGNMENT_PAR_" + SSTR(i)),
                   par.createParameterString(out, par.align));
      }
      par.pca = 0.0;
      cmd.addVar(std::string("PROFILE_PAR_" + SSTR(i)),
                 par.createParameterString(out, par.result2profile));
      par.pca = 1.0;
    }

    FileUtil::writeFile(out, tmpDir + "/blastpgp.sh", blastpgp_sh, blastpgp_sh_len);
    program = std::string(tmpDir + "/blastpgp.sh");
  } else {
    if (par.sensSteps > 1) {
      if (par.startSens > par.sensitivity) {
        out->failure("--start-sens should not be greater -s");
      }
      cmd.addVar("SENSE_0", SSTR(par.startSens));
      senses.push_back(par.startSens);
      float sensStepSize = (par.sensitivity - par.startSens) /
                           (static_cast<float>(par.sensSteps) - 1);
      for (int step = 1; step < par.sensSteps; step++) {
        std::string stepKey = "SENSE_" + SSTR(step);
        float stepSense = par.startSens + sensStepSize * step;
        senses.push_back(stepSense);
        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << stepSense;
        std::string value = stream.str();
        cmd.addVar(stepKey.c_str(), value);
      }
      no_steps = par.sensSteps;
      cmd.addVar("STEPS", SSTR((int)par.sensSteps));
    } else {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(1) << par.sensitivity;
      std::string sens = stream.str();
      cmd.addVar("SENSE_0", sens);
      cmd.addVar("STEPS", SSTR(1));
      senses.push_back(par.sensitivity);
      no_steps = 1;
    }

    std::vector<BiosnakeParameter *> prefilterWithoutS;
    for (size_t i = 0; i < par.prefilter.size(); i++) {
      if (par.prefilter[i]->uniqid != par.PARAM_S.uniqid) {
        prefilterWithoutS.push_back(par.prefilter[i]);
      }
    }
    cmd.addVar("PREFILTER_PAR", par.createParameterString(out, prefilterWithoutS));
    if (isUngappedMode) {
      par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
      cmd.addVar("ALIGNMENT_PAR",
                 par.createParameterString(out, par.rescorediagonal));
      par.rescoreMode = originalRescoreMode;
    } else {
      cmd.addVar("ALIGNMENT_PAR", par.createParameterString(out, par.align));
    }
    FileUtil::writeFile(out, tmpDir + "/blastp.sh", blastp_sh, blastp_sh_len);
    program = std::string(tmpDir + "/blastp.sh");
  }

  if (searchMode & (Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED |
                    Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)) {
    cmd.addVar("NO_TARGET_INDEX", (indexStr == "") ? "TRUE" : "");
    FileUtil::writeFile(out, tmpDir + "/translated_search.sh", translated_search_sh,
                        translated_search_sh_len);
    cmd.addVar("QUERY_NUCL",
               (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED)
                   ? "TRUE"
                   : "");
    cmd.addVar("TARGET_NUCL",
               (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)
                   ? "TRUE"
                   : "");
    cmd.addVar("THREAD_COMP_PAR",
               par.createParameterString(out, par.threadsandcompression).c_str());
    par.subDbMode = 1;
    cmd.addVar("CREATESUBDB_PAR", par.createParameterString(out, par.createsubdb));
    par.translate = 1;
    cmd.addVar("ORF_PAR", par.createParameterString(out, par.extractorfs));
    cmd.addVar("OFFSETALIGNMENT_PAR",
               par.createParameterString(out, par.offsetalignment));
    cmd.addVar("SEARCH", program);
    search_prog = program;
    program = std::string(tmpDir + "/translated_search.sh");
  } else if (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE &&
             searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE) {
    FileUtil::writeFile(out, tmpDir + "/blastn.sh", blastn_sh, blastn_sh_len);
    //  0: reverse, 1: forward, 2: both
    switch (par.strand) {
      case 0:
        par.forwardFrames = "";
        par.reverseFrames = "1";
        cmd.addVar("EXTRACTFRAMES", "TRUE");
        extract_frames = true;
        break;
      case 1:
        par.forwardFrames = "1";
        par.reverseFrames = "";
        break;
      case 2:
        par.forwardFrames = "1";
        par.reverseFrames = "1";
        cmd.addVar("EXTRACTFRAMES", "TRUE");
        extract_frames = true;
        break;
    }
    cmd.addVar("SPLITSEQUENCE_PAR",
               par.createParameterString(out, par.splitsequence));
    if (indexStr == "") {
      cmd.addVar("NEEDTARGETSPLIT", "TRUE");
      needTargetSplit = true;
    }
    cmd.addVar("NEEDQUERYSPLIT", "TRUE");
    need_query_split = true;
    cmd.addVar("EXTRACT_FRAMES_PAR",
               par.createParameterString(out, par.extractframes));
    cmd.addVar("OFFSETALIGNMENT_PAR",
               par.createParameterString(out, par.offsetalignment));
    cmd.addVar("SEARCH", program);
    search_prog = program;
    program = std::string(tmpDir + "/blastn.sh");
  }

  std::string query = par.filenames[0];
  std::string target = par.filenames[1];
  std::string result = par.filenames[2];
  std::string searchTmpDir = par.baseTmpPath + parentTmpDir + "/search";
  // mkdir tmpDir
  FileUtil::makeDir(out, searchTmpDir.c_str());

  if (program == tmpDir + "/blastp.sh") {
    call_blastp(out, par, no_steps, senses, align_module, query, target, result,
                searchTmpDir, tmpDir + "/blastp.sh", cmd);
    return 0;
  }

  //////////

  out->info("step_search A");
  if (needTargetSplit) {
    if (!FileUtil::fileExists(out, (tmpDir + "/target_seqs_split.dbtype").c_str())) {
      out->info("step_search B");
      Parameters splitsequence_par(par);

      // OFFSETALIGNMENT_PAR: [--chain-alignments 0 --merge-query 1
      // --search-type 0 --threads 1 --compressed 0 --db-load-mode 0 -v 3 ]
      // EXTRACT_FRAMES_PAR: [--forward-frames 1,2,3 --reverse-frames 1,2,3
      // --create-lookup 0 --threads 1 --compressed 0 -v 3 ]
      // SPLITSEQUENCE_PAR: [--max-seq-len 65535 --sequence-overlap 0
      // --sequence-split-mode 1 --headers-split-mode 0 --create-lookup 0
      // --threads 1 --compressed 0 -v 3 ]
      splitsequence_par.setDBFields(1, target);
      splitsequence_par.setDBFields(2, tmpDir + "/target_seqs_split");
      splitsequence_par.maxSeqLen = 65535;
      splitsequence_par.sequenceOverlap = 0;
      splitsequence_par.sequenceSplitMode = 1;
      splitsequence_par.headerSplitMode = 0;
      splitsequence_par.createLookup = 0;
      splitsequence_par.compressed = 0;
      subcall_biosnake(out, "splitsequence", splitsequence_par);
    }
    target = tmpDir + "/target_seqs_split";
    out->info("step_search C");
  }

  out->info("step_search D");
  if (extract_frames) {
    if (!FileUtil::fileExists(out, (tmpDir + "/query_seqs.dbtype").c_str())) {
      out->info("step_search E");
      Parameters extractframes_par(par);
      std::vector<std::string> extractframes_filenames = {
          query, tmpDir + "/query_seqs"};
      extractframes_par.filenames = extractframes_filenames;
      // EXTRACT_FRAMES_PAR: [--forward-frames 1,2,3 --reverse-frames 1,2,3
      // --create-lookup 0 --threads 1 --compressed 0 -v 3 ]
      extractframes_par.setDBFields(1, query);
      extractframes_par.setDBFields(2, tmpDir + "/query_seqs");
      extractframes_par.forwardFrames = par.forwardFrames;
      extractframes_par.reverseFrames = par.reverseFrames;
      extractframes_par.createLookup = 0;
      extractframes_par.compressed = 0;

      subcall_biosnake(out, "extractframes", extractframes_par);
      out->info("Call extract frames [{}] [{}]", query, (tmpDir + "/query_seqs"));
    }
    query = tmpDir + "/query_seqs";
    out->info("step_search F");
  }

  out->info("step_search G");
  if (need_query_split) {
    if (!FileUtil::fileExists(out, (tmpDir + "/query_seqs_split.dbtype").c_str())) {
      out->info("step_search H [{}] [{}]", query, (tmpDir + "/query_seqs_split"));
      Parameters splitsequence_par(par);
      splitsequence_par.setDBFields(1, query);
      splitsequence_par.setDBFields(2, tmpDir + "/query_seqs_split");
      splitsequence_par.maxSeqLen = 65535;
      splitsequence_par.sequenceOverlap = 0;
      splitsequence_par.sequenceSplitMode = 1;
      splitsequence_par.headerSplitMode = 0;
      splitsequence_par.createLookup = 0;
      splitsequence_par.compressed = 0;
      subcall_biosnake(out, "splitsequence", splitsequence_par);
    }
    query = tmpDir + "/query_seqs_split";
    out->info("step_search I");
  }

  out->info("step_search J");
  if (!FileUtil::fileExists(out, (tmpDir + "/aln.dbtype").c_str())) {
    call_blastp(out, par, no_steps, senses, align_module, query, target,
                tmpDir + "/aln", searchTmpDir, "", cmd);
  }

  out->info("step_search M");
  if (!FileUtil::fileExists(out, (result + ".dbtype").c_str())) {
    Parameters offsetalignment_par(par);
    std::vector<std::string> alignment_filenames = {
        par.filenames[0], query,           par.filenames[1],
        target,           tmpDir + "/aln", result,
    };
    offsetalignment_par.filenames = alignment_filenames;
    offsetalignment_par.setDBFields(1, par.filenames[0]);
    offsetalignment_par.setDBFields(2, query);
    offsetalignment_par.setDBFields(3, par.filenames[1]);
    offsetalignment_par.setDBFields(4, target);
    offsetalignment_par.setDBFields(5, tmpDir + "/aln");
    offsetalignment_par.setDBFields(6, result);
    offsetalignment_par.searchType = par.searchType;
    offsetalignment_par.baseTmpPath = par.baseTmpPath;
    subcall_biosnake(out, "offsetalignment", offsetalignment_par);
    out->info("biosnake offsetalignment {} {} {} {} {} {}", par.filenames[0], query, par.filenames[1], target, (tmpDir + "/aln"), result);
    out->info("step_search O");
  }

  out->info("step_search P");
  if (par.removeTmpFiles) {
    out->info("step_search Q");
    Parameters rmdb_par(par);
    rmdb_par.setDBFields(1, tmpDir + "/q_orfs");
    subcall_biosnake(out, "rmdb", rmdb_par);
    rmdb_par.setDBFields(1, tmpDir + "/q_orfs_aa");
    subcall_biosnake(out, "rmdb", rmdb_par);
    rmdb_par.setDBFields(1, tmpDir + "/t_orfs");
    subcall_biosnake(out, "rmdb", rmdb_par);
    rmdb_par.setDBFields(1, tmpDir + "/t_orfs_aa");
    subcall_biosnake(out, "rmdb", rmdb_par);
  }
  out->info("step_search R");

  // cmd.execProgram(program.c_str(), par.filenames);
  // Should never get here
  // assert(false);

  return 0;
}
