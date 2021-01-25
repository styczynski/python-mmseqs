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

#include <iostream>

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
    std::string results;

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
    std::string originalTmpDir = tmpDir;
    tmpDir = FileUtil::createTemporaryDirectory(par.baseTmpPath, tmpDir, hash);
    par.filenames.pop_back();

    out->output_string("TMP_PATH", tmpDir);
    results = par.filenames.back();
    out->output_string("RESULTS", results);
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
    std::string index_ext = "";
    if (linsearch) {
        search_module = "linsearch";
        const bool isIndex = LinsearchIndexReader::searchForIndex(target).empty() == false;
        out->output_string("INDEXEXT", isIndex ? ".linidx" : "");
        index_ext = isIndex ? ".linidx" : "";
        out->output_string("SEARCH_MODULE", "linsearch");
        out->output_string("LINSEARCH", "TRUE");
        out->output_string("CREATELININDEX_PAR", par.createParameterString(par.createlinindex));
        out->output_string("SEARCH_PAR", par.createParameterString(par.linsearchworkflow, true));
    } else {
        search_module = "search";
        const bool isIndex = PrefilteringIndexReader::searchForIndex(target).empty() == false;
        out->output_string("INDEXEXT", isIndex ? ".idx" : "");
        index_ext = isIndex ? ".idx" : "";
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
        Parameters createdb_par(par);
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
            Parameters createdb_par(par);
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
            Parameters createlinindex_par(par);
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
            //createlinindex_par.removeTmpFiles=true;
            subcall_mmseqs(out, "createlinindex", createlinindex_par);
        }
    }

    std::string intermediate = tmpDir + "/result";
    if (!FileUtil::fileExists((intermediate + ".dbtype").c_str())) {
        //search_module
        Parameters search_par(par);
        std::vector<std::string> search_filenames = {
            tmpDir + "/query",
            target,
            intermediate,
            originalTmpDir,
        };
        search_par.filenames = search_filenames;
        search_par.setDBFields(1, tmpDir + "/query");
        search_par.setDBFields(2, target);
        search_par.setDBFields(3, intermediate);
        search_par.setDBFields(4, originalTmpDir);
        search_par.spacedKmer = true;
        search_par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
        search_par.sensitivity = 5.7;
        search_par.evalThr = 0.001;
        search_par.orfStartMode = 1;
        search_par.orfMinLength = 30;
        search_par.orfMaxLength = 32734;
        search_par.evalProfile = 0.1;
        search_par.baseTmpPath = par.baseTmpPath;
        search_par.searchType = par.searchType;
        search_par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
        search_par.exactKmerMatching = true;
        search_par.strand = 2;
        search_par.kmerSize = 15;
        search_par.maxSeqLen = 10000;

        std::cout << "Call search\n" << std::flush;
        subcall_mmseqs(out, search_module, search_par);
        std::cout << "Calling search doneA\n" << std::flush;
    }
    
    /*
    if [ -n "${GREEDY_BEST_HITS}" ]; then
        if notExists "${TMP_PATH}/result_best.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" summarizeresult "${TMP_PATH}/result" "${TMP_PATH}/result_best" ${SUMMARIZE_PAR} \
                || fail "Search died"
        fi
        INTERMEDIATE="${TMP_PATH}/result_best"
    fi
    */
    out->print();
    if (false) {
        std::cout << "Call summarizeresult\n" << std::flush;
        Parameters summarizeresult_par(par);
        std::vector<std::string> summarizeresult_filenames = {
            tmpDir + "/result",
            tmpDir + "/result_best"
        };
        summarizeresult_par.filenames = summarizeresult_filenames;
        summarizeresult_par.setDBFields(1, tmpDir + "/result");
        summarizeresult_par.setDBFields(2, tmpDir + "/result_best");
        subcall_mmseqs(out, "summarizeresult", summarizeresult_par);
        std::cout << "Call summarizeresult ended\n" << std::flush;
    }
    
    // "$MMSEQS" convertalis "${TMP_PATH}/query" "${TARGET}${INDEXEXT}" "${INTERMEDIATE}" "${RESULTS}" ${CONVERT_PAR}
    if (true) {
        std::cout << "Call convertalis\n" << std::flush;
        //--db-output 0 --db-load-mode 0 --search-type 3 --threads 1 --compressed 0 -v 3 ]
        Parameters convertalis_par(par);
        std::vector<std::string> convertalis_filenames = {
            tmpDir + "/query",
            target + index_ext,
            intermediate,
            results,
        };
        out->print();
        std::cout << "convertalis " << (tmpDir + "/query") << " " << (target+index_ext) << " " << intermediate << " " << results << "\n" << std::flush;
        convertalis_par.filenames = convertalis_filenames;
        convertalis_par.setDBFields(1, tmpDir + "/query");
        convertalis_par.setDBFields(2, target+index_ext);
        convertalis_par.setDBFields(3, intermediate);
        convertalis_par.setDBFields(4, results);
        convertalis_par.setSubstitutionMatrices("blosum62.out", "nucleotide.out");
        convertalis_par.formatAlignmentMode = 0;
        convertalis_par.outfmt = par.outfmt;
        convertalis_par.translationTable = 1;
        convertalis_par.gapOpen = MultiParam<int>(11, 5);
        convertalis_par.gapExtend = MultiParam<int>(1, 2);
        convertalis_par.dbOut = 0;
        convertalis_par.preloadMode = 0;
        convertalis_par.threads = 1;
        convertalis_par.compressed = 0;

        subcall_mmseqs(out, "convertalis", convertalis_par);
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
