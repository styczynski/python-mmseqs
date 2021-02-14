#include <biosnake/commons/commandCaller.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#include "enrich.sh.h"

#include <cassert>

void setEnrichWorkflowDefaults(Parameters *p) {
  p->numIterations = 3;
  p->expansionMode = 1;
}

int enrich(biosnake_output *out, Parameters &par) {
  //    Parameters &par = Parameters::getInstance();
  //    setEnrichWorkflowDefaults(&par);
  ////    par.parseParameters(argc, argv, command, true, 0, 0);

  std::string tmpDir = par.db6;
  std::string hash = SSTR(par.hashParameter(out, par.databases_types, par.filenames,
                                            par.enrichworkflow));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();
  par.filenames.push_back(tmpDir);

  CommandCaller cmd(out);
  cmd.addVariable("RUNNER", par.runner.c_str());
  cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
  cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);

  par.addBacktrace = true;

  int originalNumIterations = par.numIterations;
  par.numIterations = 1;
  par.sliceSearch = true;
  cmd.addVariable("PROF_SEARCH_PAR",
                  par.createParameterString(out, par.searchworkflow).c_str());
  par.sliceSearch = false;
  par.numIterations = originalNumIterations;

  cmd.addVariable("PROF_PROF_PAR",
                  par.createParameterString(out, par.result2profile).c_str());
  cmd.addVariable("SUBSTRACT_PAR",
                  par.createParameterString(out, par.subtractdbs).c_str());
  cmd.addVariable("VERBOSITY_PAR",
                  par.createParameterString(out, par.onlyverbosity).c_str());

  // change once rescorediagonal supports profiles
  const bool isUngappedMode = false;
  cmd.addVariable("ALIGN_MODULE", "align");

  double originalEval = par.evalThr;
  par.evalThr = par.evalProfile;
  par.realign = false;
  for (int i = 0; i < par.numIterations; i++) {
    if (i == (par.numIterations - 1)) {
      par.evalThr = originalEval;
    }

    cmd.addVariable(std::string("PREFILTER_PAR_" + SSTR(i)).c_str(),
                    par.createParameterString(out, par.prefilter).c_str());
    if (isUngappedMode) {
      int originalRescoreMode = par.rescoreMode;
      par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
      cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(),
                      par.createParameterString(out, par.rescorediagonal).c_str());
      par.rescoreMode = originalRescoreMode;
    } else {
      cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(),
                      par.createParameterString(out, par.align).c_str());
    }

    cmd.addVariable(std::string("EXPAND_PAR_" + SSTR(i)).c_str(),
                    par.createParameterString(out, par.expandaln).c_str());

    par.pca = 0.0;
    cmd.addVariable(std::string("PROFILE_PAR_" + SSTR(i)).c_str(),
                    par.createParameterString(out, par.result2profile).c_str());
    par.pca = 1.0;
  }

  std::string program = tmpDir + "/enrich.sh";
  FileUtil::writeFile(out, program, enrich_sh, enrich_sh_len);
  cmd.execProgram(program.c_str(), par.filenames);

  // Should never get here
  assert(false);
  return 0;
}
