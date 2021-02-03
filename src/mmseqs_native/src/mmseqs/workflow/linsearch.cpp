#include <mmseqs/commons/commandCaller.h>
#include <mmseqs/commons/dBReader.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/linclust/linsearchIndexReader.h>
#include <mmseqs/prefiltering/prefilteringIndexReader.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

#include "linsearch.sh.h"

namespace Linsearch {
#include "translated_search.sh.h"
}

#include <cassert>
#include <climits>

void setLinsearchDefaults(Parameters *p) {
  p->spacedKmer = false;
  p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
  p->sensitivity = 5.7;
  p->evalThr = 0.001;
  p->maskMode = 0;
  p->orfStartMode = 1;
  p->orfMinLength = 30;
  p->orfMaxLength = 32734;
  p->evalProfile = 0.1;

  // VTML has a slightly lower sensitivity in the regression test
  p->seedScoringMatrixFile =
      MultiParam<char *>("blosum62.out", "nucleotide.out");
}

int linsearch(mmseqs_output *out, Parameters &par) {
  const int queryDbType = FileUtil::parseDbType(out, par.db1.c_str());
  std::string indexStr = LinsearchIndexReader::searchForIndex(par.db2);
  if (indexStr.empty()) {
    out->failure("Database {} needs to be index: createlinindex {}", par.db2, par.db2);
  }
  int targetDbType = 0;
  if (indexStr != "") {
    DBReader<unsigned int> dbr(
        indexStr.c_str(), (indexStr + ".index").c_str(), par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::NOSORT);
    PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&dbr);
    targetDbType = data.seqType;
    dbr.close();
  }

  if (queryDbType == -1 || targetDbType == -1) {
    out->failure("Please recreate your database or add a .dbtype file to your sequence/profile database.");
  }

  if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
      Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)) {
    out->failure("Profile-Profile searches are not supported.");
  }

  const bool isNuclSearch =
      (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
       Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES));
  //    if(isNuclSearch == true){
  //        setNuclSearchDefaults(&par);
  //    }else{
  //        par.overrideParameterDescription((Command &) command,
  //        par.PARAM_STRAND.uniqid, NULL, NULL,
  //                                         par.PARAM_STRAND.category |
  //                                         MMseqsParameter::COMMAND_EXPERT);
  //    }

  par.filenames[1] = indexStr;
  const bool isTranslatedNuclSearch =
      isNuclSearch == false &&
      (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) ||
       Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES));

  const bool isUngappedMode =
      par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
  if (isUngappedMode &&
      (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) ||
       Parameters::isEqualDbtype(targetDbType,
                                 Parameters::DBTYPE_HMM_PROFILE))) {
    // par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN |
    // MMseqsParameter::COMMAND_PREFILTER);
    out->failure("Cannot use ungapped alignment mode with profile databases.");
  }

  // par.printParameters(command.cmd, argc, argv, par.searchworkflow);

  std::string tmpDir = par.db4;
  std::string hash = SSTR(par.hashParameter(par.databases_types, par.filenames,
                                            par.linsearchworkflow));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();
  par.filenames.push_back(tmpDir);

  CommandCaller cmd;
  cmd.addVariable("FILTER", "1");
  int oldCovMode = par.covMode;
  par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
  if (par.PARAM_COV_MODE.wasSet == false) {
    par.covMode = Parameters::COV_MODE_TARGET;
  }
  float oldCov = par.covThr;
  par.covThr = std::max(par.covThr, 0.9f);
  cmd.addVariable("RESCORE_FILTER_PAR",
                  par.createParameterString(par.rescorediagonal).c_str());
  par.covMode = oldCovMode;
  par.covThr = oldCov;

  cmd.addVariable("ALIGN_MODULE", isUngappedMode ? "rescorediagonal" : "align");
  cmd.addVariable("KMERSEARCH_PAR",
                  par.createParameterString(par.kmersearch).c_str());
  double oldEval = par.evalThr;
  par.evalThr = 100000;
  cmd.addVariable("ALIGNMENT_PAR",
                  par.createParameterString(par.align).c_str());
  par.evalThr = oldEval;
  cmd.addVariable("SWAPRESULT_PAR",
                  par.createParameterString(par.swapresult).c_str());
  cmd.addVariable("NUCL", isNuclSearch ? "1" : NULL);

  std::string program = tmpDir + "/linsearch.sh";
  FileUtil::writeFile(out, program, linsearch_sh, linsearch_sh_len);

  if (isTranslatedNuclSearch == true) {
    cmd.addVariable("NO_TARGET_INDEX", (indexStr == "") ? "TRUE" : NULL);
    cmd.addVariable(
        "QUERY_NUCL",
        Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES)
            ? "TRUE"
            : NULL);
    cmd.addVariable(
        "TARGET_NUCL",
        Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES)
            ? "TRUE"
            : NULL);
    par.translate = 1;
    cmd.addVariable("ORF_PAR",
                    par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("OFFSETALIGNMENT_PAR",
                    par.createParameterString(par.offsetalignment).c_str());
    cmd.addVariable("TRANSLATE_PAR",
                    par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("SEARCH", program.c_str());
    program = std::string(tmpDir + "/translated_search.sh");
    FileUtil::writeFile(out, program, Linsearch::translated_search_sh,
                        Linsearch::translated_search_sh_len);
  }
  cmd.execProgram(program.c_str(), par.filenames);

  // Should never get here
  assert(false);
  return 0;
}
