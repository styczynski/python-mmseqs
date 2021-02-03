#include <mmseqs/commons/commandCaller.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/prefiltering/prefilteringIndexReader.h>
#include <mmseqs/commons/util.h>

#include <mmseqs/output.h>
#include "taxonomy.sh.h"
#include "taxpercontig.sh.h"

extern int computeSearchMode(int queryDbType, int targetDbType,
                             int targetSrcDbType, int searchType);

void setTaxonomyDefaults(Parameters *p) {
  p->spacedKmer = true;
  p->sensitivity = 2;
  p->evalThr = 1;
  p->maxAccept = 30;
  p->maxRejected = 5;
  p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
  p->orfStartMode = 1;
  p->orfMinLength = 30;
  p->orfMaxLength = 32734;
  p->orfFilter = true;
}
void setTaxonomyMustPassAlong(Parameters *p) {
  p->PARAM_SPACED_KMER_MODE.wasSet = true;
  p->PARAM_S.wasSet = true;
  p->PARAM_E.wasSet = true;
  p->PARAM_MAX_ACCEPT.wasSet = true;
  p->PARAM_MAX_REJECTED.wasSet = true;
  p->PARAM_ALIGNMENT_MODE.wasSet = true;
  p->PARAM_ORF_START_MODE.wasSet = true;
  p->PARAM_ORF_MIN_LENGTH.wasSet = true;
  p->PARAM_ORF_MAX_LENGTH.wasSet = true;
}

int taxonomy(mmseqs_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //
  //    for (size_t i = 0; i < par.searchworkflow.size(); i++) {
  //        par.searchworkflow[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
  //    }
  //    par.PARAM_S.removeCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_E.removeCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
  //    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);
  //
  //    setTaxonomyDefaults(&par);
  //    par.parseParameters(argc, argv, command, true, 0, 0);
  setTaxonomyMustPassAlong(&par);

  if (par.taxonomySearchMode == Parameters::TAXONOMY_2BLCA) {
    out->warn("2bLCA was replaced by Accelerated 2bLCA");
    par.taxonomySearchMode = Parameters::TAXONOMY_ACCEL_2BLCA;
  }

  std::string indexStr = PrefilteringIndexReader::searchForIndex(par.db2);
  int targetDbType = FileUtil::parseDbType(out, par.db2.c_str());
  std::string targetDB = (indexStr == "") ? par.db2.c_str() : indexStr.c_str();
  int targetSrcDbType = -1;
  if (indexStr != "" ||
      Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_INDEX_DB)) {
    indexStr = par.db2;
    DBReader<unsigned int> dbr(
        targetDB.c_str(), (targetDB + ".index").c_str(), par.threads,
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

  int searchMode = computeSearchMode(queryDbType, targetDbType, targetSrcDbType,
                                     par.searchType);
  if ((searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE) &&
      (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE)) {
    if (par.taxonomySearchMode == Parameters::TAXONOMY_ACCEL_2BLCA) {
      out->warn("Accel. 2bLCA cannot be used with nucl-nucl taxonomy, using top-hit instead");
      par.taxonomySearchMode = Parameters::TAXONOMY_TOP_HIT;
    }
  }

  std::string tmpDir = par.db4;
  std::string hash =
      SSTR(par.hashParameter(par.databases_types, par.filenames, par.taxonomy));
  if (par.reuseLatest) {
    hash = FileUtil::getHashFromSymLink(out, tmpDir + "/latest");
  }
  tmpDir = FileUtil::createTemporaryDirectory(out, par.baseTmpPath, tmpDir, hash);
  par.filenames.pop_back();
  par.filenames.push_back(tmpDir);

  CommandCaller cmd;
  std::string program;
  cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
  cmd.addVariable("RUNNER", par.runner.c_str());
  cmd.addVariable("THREADS_COMP_PAR",
                  par.createParameterString(par.threadsandcompression).c_str());
  cmd.addVariable("VERBOSITY",
                  par.createParameterString(par.onlyverbosity).c_str());

  if (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED &&
      !(searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)) {
    cmd.addVariable("TARGETDB_IDX", targetDB.c_str());
    par.translate = 1;
    cmd.addVariable("EXTRACT_ORFS_PAR",
                    par.createParameterString(par.extractorfs).c_str());
    int showTaxLineageOrig = par.showTaxLineage;
    // never show lineage for the orfs
    par.showTaxLineage = 0;
    par.PARAM_TAXON_ADD_LINEAGE.wasSet = true;
    int taxonomyOutputMode = par.taxonomyOutputMode;
    par.taxonomyOutputMode = Parameters::TAXONOMY_OUTPUT_BOTH;
    par.PARAM_TAX_OUTPUT_MODE.wasSet = true;
    cmd.addVariable("TAXONOMY_PAR",
                    par.createParameterString(par.taxonomy, true).c_str());
    par.showTaxLineage = showTaxLineageOrig;
    par.taxonomyOutputMode = taxonomyOutputMode;
    cmd.addVariable("AGGREGATETAX_PAR",
                    par.createParameterString(par.aggregatetax).c_str());
    cmd.addVariable("SWAPDB_PAR",
                    par.createParameterString(par.swapdb).c_str());

    cmd.addVariable(
        "ORF_FILTER",
        par.orfFilter && par.orfFilterSens <= par.sensitivity ? "TRUE" : NULL);
    par.minDiagScoreThr = 3;
    par.sensitivity = par.orfFilterSens;
    par.diagonalScoring = false;
    par.maxResListLen = 1;
    cmd.addVariable("ORF_FILTER_PREFILTER",
                    par.createParameterString(par.prefilter).c_str());

    par.evalThr = par.orfFilterEval;
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    cmd.addVariable("ORF_FILTER_RESCOREDIAGONAL",
                    par.createParameterString(par.rescorediagonal).c_str());

    par.subDbMode = Parameters::SUBDB_MODE_SOFT;
    cmd.addVariable("CREATESUBDB_PAR",
                    par.createParameterString(par.createsubdb).c_str());

    program = tmpDir + "/taxpercontig.sh";
    FileUtil::writeFile(out, program, taxpercontig_sh, taxpercontig_sh_len);
  } else {
    if (par.taxonomySearchMode == Parameters::TAXONOMY_TOP_HIT) {
      cmd.addVariable("TOPHIT_MODE", "1");
    } else if (par.taxonomySearchMode == Parameters::TAXONOMY_ACCEL_2BLCA) {
      par.lcaSearch = true;
      par.PARAM_LCA_SEARCH.wasSet = true;
      cmd.addVariable("TOPHIT_MODE", NULL);
    }
    cmd.addVariable(
        "SEARCH_PAR",
        par.createParameterString(par.searchworkflow, true).c_str());

    program = tmpDir + "/taxonomy.sh";
    FileUtil::writeFile(out, program.c_str(), taxonomy_sh, taxonomy_sh_len);
  }
  if (par.taxonomyOutputMode == Parameters::TAXONOMY_OUTPUT_LCA) {
    cmd.addVariable("TAX_OUTPUT", "0");
    cmd.addVariable("LCA_PAR", par.createParameterString(par.lca).c_str());
  } else if (par.taxonomyOutputMode == Parameters::TAXONOMY_OUTPUT_BOTH) {
    cmd.addVariable("TAX_OUTPUT", "2");
    cmd.addVariable("LCA_PAR", par.createParameterString(par.lca).c_str());
  } else {
    cmd.addVariable("TAX_OUTPUT", "1");
  }
  cmd.execProgram(program.c_str(), par.filenames);

  return EXIT_SUCCESS;
}
