#include <biosnake/commons/commandCaller.h>
#include <biosnake/commons/dBWriter.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#include "linclust.sh.h"

#include <cassert>

void setLinclustWorkflowDefaults(Parameters *p) {
  p->spacedKmer = false;
  p->covThr = 0.8;
  p->maskMode = 0;
  p->evalThr = 0.001;
  p->seqIdThr = 0.9;
  p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
}

int linclust(biosnake_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    setLinclustWorkflowDefaults(&par);
  //    par.PARAM_ADD_BACKTRACE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_ALT_ALIGNMENT.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_ZDROP.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_RESCORE_MODE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MAX_REJECTED.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_MAX_ACCEPT.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.overrideParameterDescription(par.PARAM_S, "Sensitivity will be
  //    automatically determined but can be adjusted", NULL,
  //    par.PARAM_S.category | BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  std::string tmpDir = par.db3;
  std::string hash = SSTR(par.hashParameter(out, par.databases_types, par.filenames,
                                            par.linclustworkflow));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();
  par.filenames.push_back(tmpDir);

  CommandCaller cmd(out);
  cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
  cmd.addVariable("RUNNER", par.runner.c_str());

  // save some values to restore them later
  MultiParam<int> alphabetSize = par.alphabetSize;
  size_t kmerSize = par.kmerSize;
  // # 1. Finding exact $k$-mer matches.
  bool kmerSizeWasSet = false;
  bool alphabetSizeWasSet = false;
  bool clusterModeSet = false;
  for (size_t i = 0; i < par.linclustworkflow.size(); i++) {
    if (par.linclustworkflow[i]->uniqid == par.PARAM_K.uniqid &&
        par.linclustworkflow[i]->wasSet) {
      kmerSizeWasSet = true;
    }
    if (par.linclustworkflow[i]->uniqid == par.PARAM_ALPH_SIZE.uniqid &&
        par.linclustworkflow[i]->wasSet) {
      alphabetSizeWasSet = true;
    }
    if (par.linclustworkflow[i]->uniqid == par.PARAM_CLUSTER_MODE.uniqid &&
        par.linclustworkflow[i]->wasSet) {
      clusterModeSet = true;
    }
  }

  const bool nonSymetric = (par.covMode == Parameters::COV_MODE_TARGET ||
                            par.covMode == Parameters::COV_MODE_QUERY);
  if (clusterModeSet == false) {
    if (nonSymetric) {
      par.clusteringMode = Parameters::GREEDY_MEM;
    } else {
      par.clusteringMode = Parameters::SET_COVER;
    }
    std::string cluMode = (par.clusteringMode == Parameters::GREEDY_MEM)
                              ? "GREEDY MEM"
                              : "SET COVER";
    out->info("Set cluster mode {}.", cluMode);
  }

  if (kmerSizeWasSet == false) {
    par.kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
  }
  if (alphabetSizeWasSet == false) {
    par.alphabetSize =
        MultiParam<int>(Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE, 5);
  }

  const int dbType = FileUtil::parseDbType(out, par.db1.c_str());
  const bool isUngappedMode =
      par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
  if (isUngappedMode &&
      Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_HMM_PROFILE)) {
    // par.printUsageMessage(command,
    // BiosnakeParameter::COMMAND_ALIGN|BiosnakeParameter::COMMAND_PREFILTER);
    out->failure("Cannot use ungapped alignment mode with profile databases");
  }

  cmd.addVariable("ALIGN_MODULE", isUngappedMode ? "rescorediagonal" : "align");
  // filter by diagonal in case of AA (do not filter for nucl, profiles, ...)
  cmd.addVariable(
      "FILTER",
      Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_AMINO_ACIDS) ? "1"
                                                                        : NULL);
  cmd.addVariable("KMERMATCHER_PAR",
                  par.createParameterString(out, par.kmermatcher).c_str());
  cmd.addVariable("VERBOSITY",
                  par.createParameterString(out, par.onlyverbosity).c_str());
  cmd.addVariable("VERBOSITYANDCOMPRESS",
                  par.createParameterString(out, par.threadsandcompression).c_str());

  par.alphabetSize = alphabetSize;
  par.kmerSize = kmerSize;

  // # 2. Hamming distance pre-clustering
  par.rescoreMode = Parameters::RESCORE_MODE_HAMMING;
  par.filterHits = false;
  float prevSeqId = par.seqIdThr;
  // hamming distance does not work well with seq. id < 0.5 since it does not
  // have an e-value criteria
  par.seqIdThr = std::max(0.5f, par.seqIdThr);
  // also coverage should not be under 0.5
  float prevCov = par.covThr;
  par.covThr = std::max(0.5f, par.covThr);
  cmd.addVariable("HAMMING_PAR",
                  par.createParameterString(out, par.rescorediagonal).c_str());
  // set it back to old value
  par.covThr = prevCov;
  par.seqIdThr = prevSeqId;
  par.rescoreMode = Parameters::RESCORE_MODE_SUBSTITUTION;

  // # 3. Ungapped alignment filtering
  par.filterHits = true;
  cmd.addVariable("UNGAPPED_ALN_PAR",
                  par.createParameterString(out, par.rescorediagonal).c_str());

  // # 4. Local gapped sequence alignment.
  if (isUngappedMode) {
    const int originalRescoreMode = par.rescoreMode;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("ALIGNMENT_PAR",
                    par.createParameterString(out, par.rescorediagonal).c_str());
    par.rescoreMode = originalRescoreMode;
  } else {
    cmd.addVariable("ALIGNMENT_PAR",
                    par.createParameterString(out, par.align).c_str());
  }
  // # 5. Clustering using greedy set cover.
  cmd.addVariable("CLUSTER_PAR", par.createParameterString(out, par.clust).c_str());
  cmd.addVariable("MERGECLU_PAR",
                  par.createParameterString(out, par.threadsandcompression).c_str());

  std::string program = tmpDir + "/linclust.sh";
  FileUtil::writeFile(out, program, linclust_sh, linclust_sh_len);
  cmd.execProgram(program.c_str(), par.filenames);

  // Unreachable
  assert(false);
  return 0;
}
