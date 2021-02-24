#include <cassert>
#include <biosnake/commons/application.h>
#include <biosnake/commons/commandCaller.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/linclust/linsearchIndexReader.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/prefiltering/prefilteringIndexReader.h>
#include <biosnake/commons/util.h>
#include "easysearch.sh.h"
#include <biosnake/output.h>

#include <iostream>

void setEasySearchDefaults(Parameters *p, bool linsearch) {
  if (linsearch) {
    p->shuffleDatabase = false;
  }
  p->sensitivity = 5.7;
  p->removeTmpFiles = true;
  p->writeLookup = false;
  p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}

void setEasySearchMustPassAlong(Parameters *p, bool linsearch) {
  if (linsearch) {
    p->PARAM_SHUFFLE.wasSet = true;
  }
  p->PARAM_S.wasSet = true;
  p->PARAM_REMOVE_TMP_FILES.wasSet = true;
  p->PARAM_ALIGNMENT_MODE.wasSet = true;
}

int doeasysearch(biosnake_output *out, Parameters &par, bool linsearch) {
  std::cout << "CMDDEBUG biosnake doeasysearch " << par.createParameterString(out, par.easysearchworkflow);
  //    Parameters &par = Parameters::getInstance();
  //    par.PARAM_ADD_BACKTRACE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MAX_REJECTED.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_ZDROP.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_DB_OUTPUT.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_OVERLAP.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_DB_OUTPUT.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_RESCORE_MODE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    for (size_t i = 0; i < par.createdb.size(); i++){
  //        par.createdb[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.extractorfs.size(); i++){
  //        par.extractorfs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.translatenucs.size(); i++){
  //        par.translatenucs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.splitsequence.size(); i++) {
  //        par.splitsequence[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.result2profile.size(); i++){
  //        par.result2profile[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    par.PARAM_COMPRESSED.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_THREADS.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_V.removeCategory(BiosnakeParameter::COMMAND_EXPERT);

  //    setEasySearchDefaults(&par, linsearch);
  //    par.parseParameters(argc, argv, command, true,
  //    Parameters::PARSE_VARIADIC, 0);
  //setEasySearchMustPassAlong(&par, linsearch);

  bool needBacktrace = false;
  bool needTaxonomy = false;
  bool needTaxonomyMapping = false;
  bool needLookup = false;
  std::string results;

  {
    bool needSequenceDB = false;
    bool needFullHeaders = false;
    bool needSource = false;
    Parameters::getOutputFormat(out, par.formatAlignmentMode, par.outfmt,
                                needSequenceDB, needBacktrace, needFullHeaders,
                                needLookup, needSource, needTaxonomyMapping,
                                needTaxonomy);
  }

  if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM ||
      par.greedyBestHits) {
    needBacktrace = true;
  }

  if (needBacktrace) {
    out->info("Alignment backtraces will be computed, since they were requested by output format");
    par.addBacktrace = true;
    par.PARAM_ADD_BACKTRACE.wasSet = true;
  }

  if (needLookup) {
    par.writeLookup = true;
  }

  std::string tmpDir = par.filenames.back();
  // TODO: Fix

  std::string hash = "abc";  // SSTR(par.hashParameter(par.databases_types,
                             // par.filenames, *command.params));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  std::string originalTmpDir = tmpDir;
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();

  CommandCaller cmd(out);
  cmd.addVariableStr("BIOSNAKE", "mmseqs");
  cmd.addVariableStr("TMP_PATH", tmpDir);
  results = par.filenames.back();
  cmd.addVariableStr("RESULTS", results);
  par.filenames.pop_back();
  std::string target = par.filenames.back();
  cmd.addVariableStr("TARGET", target);
  par.filenames.pop_back();

  if (needTaxonomy || needTaxonomyMapping) {
    std::vector<std::string> missingFiles =
        Parameters::findMissingTaxDbFiles(out, target);
    if (missingFiles.empty() == false) {
      Parameters::printTaxDbError(out, target, missingFiles);
      return 1;
    }
  }

  std::string search_module = "";
  std::string index_ext = "";
  if (linsearch) {
    search_module = "linsearch";
    const bool isIndex =
        LinsearchIndexReader::searchForIndex(out, target).empty() == false;
    cmd.addVariableStr("INDEXEXT", isIndex ? ".linidx" : "");
    index_ext = isIndex ? ".linidx" : "";
    cmd.addVariableStr("SEARCH_MODULE", "linsearch");
    cmd.addVariableStr("LINSEARCH", "TRUE");
    cmd.addVariableStr("CREATELININDEX_PAR",
                       par.createParameterString(out, par.createlinindex));
    auto search_par = par.createParameterString(out, par.linsearchworkflow, true);
    std::cout << "LIN_SEARCH_PAR = [" << search_par << "]\n";
    cmd.addVariableStr("SEARCH_PAR", search_par);
  } else {
    search_module = "search";
    const bool isIndex =
        PrefilteringIndexReader::searchForIndex(out, target).empty() == false;
    cmd.addVariableStr("INDEXEXT", isIndex ? ".idx" : "");
    index_ext = isIndex ? ".idx" : "";
    cmd.addVariableStr("SEARCH_MODULE", "search");
    cmd.addVariableStr("LINSEARCH", "");
    cmd.addVariableStr("CREATELININDEX_PAR", "");
    auto search_par = par.createParameterString(out, par.searchworkflow, true);
    std::cout << "LIN_SEARCH_PAR = [" << search_par << "]\n";
    cmd.addVariableStr("SEARCH_PAR", search_par);
  }
  cmd.addVariableStr("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : "");
  cmd.addVariableStr("GREEDY_BEST_HITS", par.greedyBestHits ? "TRUE" : "");
  cmd.addVariableStr("LEAVE_INPUT", par.dbOut ? "TRUE" : "");

  cmd.addVariableStr("RUNNER", par.runner);
  cmd.addVariableStr("VERBOSITY", par.createParameterString(out, par.onlyverbosity));

  cmd.addVariableStr("CREATEDB_QUERY_PAR",
                     par.createParameterString(out, par.createdb));
  par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
  cmd.addVariableStr("CREATEDB_PAR", par.createParameterString(out, par.createdb));
  cmd.addVariableStr("CONVERT_PAR",
                     par.createParameterString(out, par.convertalignments));
  cmd.addVariableStr("SUMMARIZE_PAR",
                     par.createParameterString(out, par.summarizeresult));

  std::string program = tmpDir + "/easysearch.sh";
  FileUtil::writeFile(out, program, easysearch_sh, easysearch_sh_len);
  std::vector<std::string> args;
  args.push_back(program);
  for (auto path : par.filenames) {
    args.push_back(path);
  }
//  cmd.execProgram("/bin/bash", args);
//  return 0;

  std::string query_file_path = par.filenames.back();

  if (!FileUtil::fileExists(out, (tmpDir + "/query.dbtype").c_str())) {
    Parameters createdb_par(par);
    std::vector<std::string> createdb_filenames = {query_file_path,
                                                   tmpDir + "/query"};
    createdb_par.filenames = createdb_filenames;
    createdb_par.dbType = 0;
    createdb_par.shuffleDatabase = 0;
    createdb_par.createdbMode = 0;
    createdb_par.identifierOffset = 0;
    createdb_par.writeLookup = 0;
    createdb_par.compressed = 0;
    subcall_biosnake(out, "createdb", createdb_par);
  }

  if (!FileUtil::fileExists(out, (target + ".dbtype").c_str())) {
    if (!FileUtil::fileExists(out, (tmpDir + "/target").c_str())) {
      Parameters createdb_par(par);
      std::vector<std::string> createdb_filenames = {target,
                                                     tmpDir + "/target"};
      createdb_par.filenames = createdb_filenames;
      createdb_par.dbType = 0;
      createdb_par.shuffleDatabase = 0;
      createdb_par.createdbMode = 0;
      createdb_par.identifierOffset = 0;
      createdb_par.writeLookup = 0;
      createdb_par.compressed = 0;

      subcall_biosnake(out, "createdb", createdb_par);
    }
    target = tmpDir + "/target";
  }

  if (linsearch) {
    if (!FileUtil::fileExists(out, (target + ".linidx").c_str())) {
      Parameters createlinindex_par(par);
      std::vector<std::string> createlinindex_filenames = {tmpDir +
                                                           "/index_tmp"};
      createlinindex_par.filenames = createlinindex_filenames;
      createlinindex_par.setDBFields(1, target);
      createlinindex_par.setDBFields(2, tmpDir + "/index_tmp");

      subcall_biosnake(out, "createlinindex", createlinindex_par);
    }
  }

  std::string intermediate = tmpDir + "/result";
  if (!FileUtil::fileExists(out, (intermediate + ".dbtype").c_str())) {
    // search_module
    Parameters search_par(par);
    std::vector<std::string> search_filenames = {
        tmpDir + "/query",
        target,
        intermediate,
        originalTmpDir,
    };
    search_par.filenames = search_filenames;
    search_par.setDBFields(1, tmpDir + "/query");
    search_par.setDBFields(2, target);
    search_par.setDBFields(3, intermediate);
    search_par.setDBFields(4, originalTmpDir);
    //search_par.maxSeqLen = 10000;

    out->info("Call search (subcall): {}", search_module);
    subcall_biosnake(out, search_module, search_par);
    out->info("Call search terminted (subcall): {}", search_module);
  }

  /*
  if [ -n "${GREEDY_BEST_HITS}" ]; then
      if notExists "${TMP_PATH}/result_best.dbtype"; then
          # shellcheck disable=SC2086
          $RUNNER "$BIOSNAKE" summarizeresult "${TMP_PATH}/result"
  "${TMP_PATH}/result_best" ${SUMMARIZE_PAR} \
              || fail "Search died"
      fi
      INTERMEDIATE="${TMP_PATH}/result_best"
  fi
  */
  if (false) {
    out->info("Call summarizeresult");
    Parameters summarizeresult_par(par);
    std::vector<std::string> summarizeresult_filenames = {
        tmpDir + "/result", tmpDir + "/result_best"};
    summarizeresult_par.filenames = summarizeresult_filenames;
    summarizeresult_par.setDBFields(1, tmpDir + "/result");
    summarizeresult_par.setDBFields(2, tmpDir + "/result_best");
    subcall_biosnake(out, "summarizeresult", summarizeresult_par);
    out->info("Call summarizeresult ended");
  }

  // "$BIOSNAKE" convertalis "${TMP_PATH}/query" "${TARGET}${INDEXEXT}"
  // "${INTERMEDIATE}" "${RESULTS}" ${CONVERT_PAR}
  if (true) {
    out->info("Call convertalis");
    //--db-output 0 --db-load-mode 0 --search-type 3 --threads 1 --compressed 0
    //-v 3 ]
    Parameters convertalis_par(par);
    std::vector<std::string> convertalis_filenames = {
        tmpDir + "/query",
        target + index_ext,
        intermediate,
        results,
    };
    out->info("convertalis {} {} {} {}", (tmpDir + "/query"), (target + index_ext), intermediate, results);
    convertalis_par.filenames = convertalis_filenames;
    convertalis_par.setDBFields(1, tmpDir + "/query");
    convertalis_par.setDBFields(2, target + index_ext);
    convertalis_par.setDBFields(3, intermediate);
    convertalis_par.setDBFields(4, results);

    subcall_biosnake(out, "convertalis", convertalis_par);
  }

  return 0;

  //    std::string program = tmpDir + "/easysearch.sh";
  //    FileUtil::writeFile(out, program, easysearch_sh, easysearch_sh_len);
  //    cmd.execProgram(program.c_str(), par.filenames);
  //
  //    // Should never get here
  //    assert(false);
  //    return EXIT_FAILURE;
}

int easysearch(biosnake_output *out, Parameters &par) {
  return doeasysearch(out, par, false);
}

int easylinsearch(biosnake_output *out, Parameters &par) {
  return doeasysearch(out, par, true);
}
