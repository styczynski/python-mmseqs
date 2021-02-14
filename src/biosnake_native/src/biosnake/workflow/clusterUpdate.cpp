#include <biosnake/commons/commandCaller.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>

#include <cassert>

#include <biosnake/output.h>
#include "update_clustering.sh.h"

extern void setClusterAutomagicParameters(biosnake_output* out, Parameters &par);

void setClusterUpdateDefaults(Parameters *p) {
  p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}
void setClusterUpdateMustPassAlong(Parameters *p) {
  p->PARAM_ALIGNMENT_MODE.wasSet = true;
}
int clusterupdate(biosnake_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    setClusterUpdateDefaults(&par);
  //    par.PARAM_ADD_BACKTRACE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_ALT_ALIGNMENT.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_RESCORE_MODE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MAX_REJECTED.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MAX_ACCEPT.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_KMER_PER_SEQ_SCALE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_KMER_PER_SEQ.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_START_SENS.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_SENS_STEPS.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_CLUSTER_REASSIGN.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //
  //    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_NUM_ITERATIONS.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    for (size_t i = 0; i < par.extractorfs.size(); i++) {
  //        par.extractorfs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.translatenucs.size(); i++) {
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
  //
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  setClusterUpdateMustPassAlong(&par);
  setClusterAutomagicParameters(out, par);

  CommandCaller cmd(out);
  cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
  cmd.addVariable("RECOVER_DELETED", par.recoverDeleted ? "TRUE" : NULL);

  cmd.addVariable("RUNNER", par.runner.c_str());
  cmd.addVariable("DIFF_PAR", par.createParameterString(out, par.diff).c_str());
  cmd.addVariable("VERBOSITY",
                  par.createParameterString(out, par.onlyverbosity).c_str());
  cmd.addVariable("THREADS_PAR",
                  par.createParameterString(out, par.onlythreads).c_str());
  cmd.addVariable("RESULT2REPSEQ_PAR",
                  par.createParameterString(out, par.result2repseq).c_str());

  cmd.addVariable("CLUST_PAR",
                  par.createParameterString(out, par.clusterworkflow, true).c_str());

  int maxAccept = par.maxAccept;
  par.maxAccept = 1;
  par.PARAM_MAX_ACCEPT.wasSet = true;
  cmd.addVariable(
      "SEARCH_PAR",
      par.createParameterString(out, par.clusterUpdateSearch, true).c_str());
  par.maxAccept = maxAccept;

  std::string tmpDir = par.db6;
  std::string hash = SSTR(
      par.hashParameter(out, par.databases_types, par.filenames, par.clusterUpdate));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();
  par.filenames.push_back(tmpDir);

  std::string program = tmpDir + "/update_clustering.sh";
  FileUtil::writeFile(out, program, update_clustering_sh, update_clustering_sh_len);
  cmd.execProgram(program.c_str(), par.filenames);

  // Should never get here
  assert(false);
  return EXIT_FAILURE;
}
