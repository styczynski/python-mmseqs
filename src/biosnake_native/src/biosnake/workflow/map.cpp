#include <biosnake/commons/commandCaller.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#include "map.sh.h"

#include <cassert>

void setMapWorkflowDefaults(Parameters *p) {
  p->compBiasCorrection = 0;
  p->maskMode = 0;
  p->covThr = 0.95;
  p->covMode = 2;
  p->seqIdThr = 0.9;
  p->sensitivity = 2;
  p->rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
  p->sortResults = true;
  // p->orfLongest = true;
  p->orfStartMode = 1;
  p->orfMinLength = 10;
  p->orfMaxLength = 32734;
}

int map(biosnake_output *out, Parameters &par) {
  //    Parameters &par = Parameters::getInstance();
  //    setMapWorkflowDefaults(&par);
  //
  //    par.PARAM_OVERLAP.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_DB_OUTPUT.addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    for (size_t i = 0; i < par.extractorfs.size(); i++){
  //        par.extractorfs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    for (size_t i = 0; i < par.translatenucs.size(); i++){
  //        par.translatenucs[i]->addCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    }
  //    par.PARAM_COMPRESSED.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_V.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //    par.PARAM_THREADS.removeCategory(BiosnakeParameter::COMMAND_EXPERT);
  //
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  std::string tmpDir = par.db4;
  std::string hash = SSTR(
      par.hashParameter(out, par.databases_types, par.filenames, par.mapworkflow));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();
  par.filenames.push_back(tmpDir);

  CommandCaller cmd(out);
  cmd.addVariable("RUNNER", par.runner.c_str());

  par.mapworkflow.push_back(&(par.PARAM_ALIGNMENT_MODE));
  par.alignmentMode = 4;
  cmd.addVariable("SEARCH_PAR",
                  par.createParameterString(out, par.mapworkflow).c_str());

  std::string program = tmpDir + "/map.sh";
  FileUtil::writeFile(out, program, map_sh, map_sh_len);
  cmd.execProgram(program.c_str(), par.filenames);

  // Should never get here
  assert(false);
  return 0;
}
