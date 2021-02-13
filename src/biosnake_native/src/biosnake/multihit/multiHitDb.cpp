#include <biosnake/commons/commandCaller.h>
#include <biosnake/commons/dBReader.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#include "multihitdb.sh.h"

void setMultiHitDbWorkflowDefaults(Parameters *p) { p->orfMinLength = 30; }

int multihitdb(biosnake_output *out, Parameters &par) {
  //    Parameters &par = Parameters::getInstance();
  //    setMultiHitDbWorkflowDefaults(&par);
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  std::string tmpDir = par.filenames.back();
  par.filenames.pop_back();

  if (FileUtil::directoryExists(out, tmpDir.c_str()) == false) {
    out->info("Tmp {} folder does not exist or is not a directory.", tmpDir
                      );
    if (FileUtil::makeDir(out, tmpDir.c_str()) == false) {
      out->failure("Can not create tmp folder {}.", tmpDir);
    } else {
      out->info("Created dir {}", tmpDir);
    }
  }
  std::string hash = SSTR(
      par.hashParameter(out, par.databases_types, par.filenames, par.multihitdb));
  if (par.reuseLatest == true) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = tmpDir + "/" + hash;
  if (FileUtil::directoryExists(out, tmpDir.c_str()) == false) {
    if (FileUtil::makeDir(out, tmpDir.c_str()) == false) {
      out->failure("Can not create sub tmp folder {}", tmpDir);
    }
  }
  FileUtil::symlinkAlias(out, tmpDir, "latest");

  std::string outDb = par.filenames.back();
  par.filenames.pop_back();

  CommandCaller cmd(out);
  cmd.addVariable("OUTDB", outDb.c_str());
  cmd.addVariable("TMP_PATH", tmpDir.c_str());

  if (par.removeTmpFiles) {
    cmd.addVariable("REMOVE_TMP", "TRUE");
  }

  cmd.addVariable("CREATEDB_PAR",
                  par.createParameterString(out, par.createdb).c_str());
  cmd.addVariable("EXTRACTORFS_PAR",
                  par.createParameterString(out, par.extractorfs).c_str());
  cmd.addVariable("TRANSLATENUCS_PAR",
                  par.createParameterString(out, par.translatenucs).c_str());
  cmd.addVariable("SWAPDB_PAR", par.createParameterString(out, par.swapdb).c_str());
  par.stat = "linecount";
  cmd.addVariable("RESULT2STATS_PAR",
                  par.createParameterString(out, par.result2stats).c_str());
  cmd.addVariable("THREADS_PAR",
                  par.createParameterString(out, par.onlythreads).c_str());

  FileUtil::writeFile(out, tmpDir + "/multihitdb.sh", multihitdb_sh,
                      multihitdb_sh_len);
  std::string program(tmpDir + "/multihitdb.sh");
  cmd.execProgram(program.c_str(), par.filenames);

  return EXIT_SUCCESS;
}
