#include <mmseqs/commons/commandCaller.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>

#include <cassert>

#include <mmseqs/output.h>
#include "update_clustering.sh.h"

extern void setClusterAutomagicParameters(Parameters &par);

void setClusterUpdateDefaults(Parameters *p) {
  p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
}
void setClusterUpdateMustPassAlong(Parameters *p) {
  p->PARAM_ALIGNMENT_MODE.wasSet = true;
}
int clusterupdate(mmseqs_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    setClusterUpdateDefaults(&par);
  //    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_ALT_ALIGNMENT.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_KMER_PER_SEQ_SCALE.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_KMER_PER_SEQ.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_START_SENS.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_SENS_STEPS.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_CLUSTER_REASSIGN.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //
  //    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_NUM_ITERATIONS.addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    for (size_t i = 0; i < par.extractorfs.size(); i++) {
  //        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.translatenucs.size(); i++) {
  //        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.splitsequence.size(); i++) {
  //        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.result2profile.size(); i++){
  //        par.result2profile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    }
  //    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);
  //
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  setClusterUpdateMustPassAlong(&par);
  setClusterAutomagicParameters(par);

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
