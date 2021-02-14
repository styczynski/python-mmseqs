#include <cassert>
#include <biosnake/commons/commandCaller.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/linclust/linsearchIndexReader.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/prefiltering/prefilteringIndexReader.h>
#include <biosnake/commons/util.h>
#include "easyrbh.sh.h"
#include <biosnake/output.h>

int easyrbh(biosnake_output* out, Parameters& par) {
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
  //    for (size_t i = 0; i < par.result2profile.size(); i++){
  //        par.result2profile[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    par.PARAM_COMPRESSED.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_THREADS.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_V.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //
  //    par.sensitivity = 5.7;
  //    par.removeTmpFiles = true;
  //    par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
  //    par.writeLookup = false;
  //    par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_SOFT;
  //    par.parseParameters(argc, argv, command, true,
  //    Parameters::PARSE_VARIADIC, 0);
  par.PARAM_S.wasSet = true;
  par.PARAM_REMOVE_TMP_FILES.wasSet = true;
  par.PARAM_ALIGNMENT_MODE.wasSet = true;

  bool needBacktrace = false;
  bool needTaxonomy = false;
  bool needTaxonomyMapping = false;
  {
    bool needSequenceDB = false;
    bool needFullHeaders = false;
    bool needLookup = false;
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
    out->info("Alignment backtraces will be computed, since they were requested by output format.");
    par.addBacktrace = true;
    par.PARAM_ADD_BACKTRACE.wasSet = true;
  }

  std::string tmpDir = par.filenames.back();
  // TODO: Fix
  std::string hash = "abc";  // SSTR(par.hashParameter(par.databases_types,
                             // par.filenames, *command.params));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();

  CommandCaller cmd(out);
  cmd.addVariable("TMP_PATH", tmpDir.c_str());
  cmd.addVariable("RESULTS", par.filenames.back().c_str());
  par.filenames.pop_back();
  std::string target = par.filenames.back().c_str();
  cmd.addVariable("TARGET", target.c_str());
  par.filenames.pop_back();

  if (needTaxonomy || needTaxonomyMapping) {
    std::vector<std::string> missingFiles =
        Parameters::findMissingTaxDbFiles(out, target);
    if (missingFiles.empty() == false) {
      Parameters::printTaxDbError(out, target, missingFiles);
      out->failure("Missing taxonomy files for {}", target);
    }
  }

  cmd.addVariable("QUERY", par.filenames.back().c_str());

  cmd.addVariable("SEARCH_PAR",
                  par.createParameterString(out, par.searchworkflow, true).c_str());
  cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
  cmd.addVariable("LEAVE_INPUT", par.dbOut ? "TRUE" : NULL);

  cmd.addVariable("RUNNER", par.runner.c_str());
  cmd.addVariable("VERBOSITY",
                  par.createParameterString(out, par.onlyverbosity).c_str());

  cmd.addVariable("CREATEDB_QUERY_PAR",
                  par.createParameterString(out, par.createdb).c_str());
  par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
  cmd.addVariable("CREATEDB_PAR",
                  par.createParameterString(out, par.createdb).c_str());
  cmd.addVariable("CONVERT_PAR",
                  par.createParameterString(out, par.convertalignments).c_str());

  std::string program = tmpDir + "/easyrbh.sh";
  FileUtil::writeFile(out, program, easyrbh_sh, easyrbh_sh_len);
  cmd.execProgram(program.c_str(), par.filenames);

  // Should never get here
  assert(false);
  return EXIT_FAILURE;
}
