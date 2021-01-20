#include <cassert>
#include "LinsearchIndexReader.h"
#include "PrefilteringIndexReader.h"
#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Application.h"
#include "Debug.h"
#include "Parameters.h"
#include "easysearch.sh.h"
#include "output.h"

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

int doeasysearch(mmseqs_output* out, Parameters &par, bool linsearch) {
//    Parameters &par = Parameters::getInstance();
//    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    for (size_t i = 0; i < par.createdb.size(); i++){
//        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
//    }
//    for (size_t i = 0; i < par.extractorfs.size(); i++){
//        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
//    }
//    for (size_t i = 0; i < par.translatenucs.size(); i++){
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

//    setEasySearchDefaults(&par, linsearch);
//    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    //setEasySearchMustPassAlong(&par, linsearch);
    
    bool needBacktrace = false;
    bool needTaxonomy = false;
    bool needTaxonomyMapping = false;
    bool needLookup = false;

    {
        bool needSequenceDB = false;
        bool needFullHeaders = false;
        bool needSource = false;
        Parameters::getOutputFormat(par.formatAlignmentMode, par.outfmt, needSequenceDB, needBacktrace, needFullHeaders,
                needLookup, needSource, needTaxonomyMapping, needTaxonomy);
    }
    

    if (par.formatAlignmentMode == Parameters::FORMAT_ALIGNMENT_SAM || par.greedyBestHits) {
        needBacktrace = true;
    }
    
    if (needBacktrace) {
        Debug(Debug::INFO) << "Alignment backtraces will be computed, since they were requested by output format.\n";
        par.addBacktrace = true;
        par.PARAM_ADD_BACKTRACE.wasSet = true;
    }
    
    if(needLookup){
        par.writeLookup = true;
    }

    std::string tmpDir = par.filenames.back();
    // TODO: Fix
    
    std::string hash = "abc"; //SSTR(par.hashParameter(par.databases_types, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(par.baseTmpPath, tmpDir, hash);
    par.filenames.pop_back();

    out->output_string("TMP_PATH", tmpDir);
    out->output_string("RESULTS", par.filenames.back());
    par.filenames.pop_back();
    std::string target = par.filenames.back();
    out->output_string("TARGET", target);
    par.filenames.pop_back();

    if (needTaxonomy || needTaxonomyMapping) {
        std::vector<std::string> missingFiles = Parameters::findMissingTaxDbFiles(target);
        if (missingFiles.empty() == false) {
            Parameters::printTaxDbError(target, missingFiles);
            return 1;
        }
    }

    std::string search_module = "";
    if (linsearch) {
        search_module = "linsearch";
        const bool isIndex = LinsearchIndexReader::searchForIndex(target).empty() == false;
        out->output_string("INDEXEXT", isIndex ? ".linidx" : "");
        out->output_string("SEARCH_MODULE", "linsearch");
        out->output_string("LINSEARCH", "TRUE");
        out->output_string("CREATELININDEX_PAR", par.createParameterString(par.createlinindex));
        out->output_string("SEARCH_PAR", par.createParameterString(par.linsearchworkflow, true));
    } else {
        search_module = "search";
        const bool isIndex = PrefilteringIndexReader::searchForIndex(target).empty() == false;
        out->output_string("INDEXEXT", isIndex ? ".idx" : "");
        out->output_string("SEARCH_MODULE", "search");
        out->output_string("LINSEARCH", "");
        out->output_string("CREATELININDEX_PAR", "");
        out->output_string("SEARCH_PAR", par.createParameterString(par.searchworkflow, true));

    }
    out->output_string("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : "");
    out->output_string("GREEDY_BEST_HITS", par.greedyBestHits ? "TRUE" : "");
    out->output_string("LEAVE_INPUT", par.dbOut ? "TRUE" : "");

    out->output_string("RUNNER", par.runner);
    out->output_string("VERBOSITY", par.createParameterString(par.onlyverbosity));

    out->output_string("CREATEDB_QUERY_PAR", par.createParameterString(par.createdb));
    par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
    out->output_string("CREATEDB_PAR", par.createParameterString(par.createdb));
    out->output_string("CONVERT_PAR", par.createParameterString(par.convertalignments));
    out->output_string("SUMMARIZE_PAR", par.createParameterString(par.summarizeresult));

    std::string query_file_path = par.filenames.back();

    if (!FileUtil::fileExists((tmpDir + "/query.dbtype").c_str())) {
        Parameters createdb_par;
        std::vector<std::string> createdb_filenames = {query_file_path, tmpDir + "/query"};
        createdb_par.filenames = createdb_filenames;
        createdb_par.dbType = 0;
        createdb_par.shuffleDatabase = 0;
        createdb_par.createdbMode = 0;
        createdb_par.identifierOffset = 0;
        createdb_par.writeLookup = 0;
        createdb_par.compressed = 0;
        subcall_mmseqs(out, "createdb", createdb_par);
    }

    if (!FileUtil::fileExists((target + ".dbtype").c_str())) {
        if (!FileUtil::fileExists((tmpDir + "/target").c_str())) {
            Parameters createdb_par;
            std::vector<std::string> createdb_filenames = {target, tmpDir + "/target"};
            createdb_par.filenames = createdb_filenames;
            createdb_par.dbType = 0;
            createdb_par.shuffleDatabase = 0;
            createdb_par.createdbMode = 0;
            createdb_par.identifierOffset = 0;
            createdb_par.writeLookup = 0;
            createdb_par.compressed = 0;
            subcall_mmseqs(out, "createdb", createdb_par);
        }
        target = tmpDir + "/target";
    }

    if (linsearch) {
        if (!FileUtil::fileExists((target + ".linidx").c_str())) {
            Parameters createlinindex_par;
            std::vector<std::string> createlinindex_filenames = {tmpDir + "/index_tmp"};
            createlinindex_par.filenames = createlinindex_filenames;
            createlinindex_par.setDBFields(1, target);
            createlinindex_par.setDBFields(2, tmpDir + "/index_tmp");
            createlinindex_par.orfStartMode=1;
            createlinindex_par.orfMinLength=30;
            createlinindex_par.orfMaxLength=32734;
            createlinindex_par.kmerScore=0;
            createlinindex_par.maskMode=1;
            createlinindex_par.sensitivity=7.5;
            createlinindex_par.removeTmpFiles=true;
            subcall_mmseqs(out, "createlinindex", createlinindex_par);
        }
    }

    std::string intermediate = tmpDir + "/result";
    if (!FileUtil::fileExists((intermediate + ".dbtype").c_str())) {
        //search_module
        Parameters search_par;
        std::vector<std::string> search_filenames = {
            tmpDir + "/query",
            target,
            intermediate,
            tmpDir + "/search_tmp",
        };
        search_par.filenames = search_filenames;
        search_par.setDBFields(1, tmpDir + "/query");
        search_par.setDBFields(2, target);
        search_par.setDBFields(3, intermediate);
        search_par.setDBFields(4, tmpDir + "/search_tmp");
        subcall_mmseqs(out, search_module, search_par);
    }

    return 0;

//    std::string program = tmpDir + "/easysearch.sh";
//    FileUtil::writeFile(program, easysearch_sh, easysearch_sh_len);
//    cmd.execProgram(program.c_str(), par.filenames);
//
//    // Should never get here
//    assert(false);
//    return EXIT_FAILURE;
}

int easysearch(mmseqs_output* out, Parameters &par) {
    return doeasysearch(out, par, false);
}

int easylinsearch(mmseqs_output* out, Parameters &par) {
    return doeasysearch(out, par, true);
}
