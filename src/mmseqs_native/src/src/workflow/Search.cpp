#include "DBReader.h"
#include "CommandCaller.h"
#include "Util.h"
#include "FileUtil.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"
#include "searchtargetprofile.sh.h"
#include "searchslicedtargetprofile.sh.h"
#include "blastpgp.sh.h"
#include "translated_search.sh.h"
#include "blastp.sh.h"
#include "blastn.sh.h"
#include "Parameters.h"
#include "output.h"
#include "Application.h"

#include <iostream>
#include <iomanip>
#include <climits>
#include <cassert>


void setSearchDefaults(Parameters *p) {
    p->spacedKmer = true;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
    p->sensitivity = 5.7;
    p->evalThr = 0.001;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->orfMaxLength = 32734;
    p->evalProfile = 0.1;
}


int computeSearchMode(int queryDbType, int targetDbType, int targetSrcDbType, int searchType) {
    // reject unvalid search
    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
        Parameters::isEqualDbtype(targetDbType,Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "Profile-Profile searches are not supported.\n";
        EXIT(EXIT_FAILURE);
    }
    // index was used
    if(targetSrcDbType == -1) {
        if(Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
           Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES))
        {
            if(searchType == Parameters::SEARCH_TYPE_AUTO){
                // WARNING because its not really an error, just a req. parameter
                Debug(Debug::WARNING) << "It is unclear from the input if a translated or nucleotide search should be performed\n"
                                         "Please provide the parameter --search-type 2 (translated), 3 (nucleotide) or 4 (translated nucleotide backtrace)\n";
                EXIT(EXIT_FAILURE);
            }
            // nucl/nucl
            // nucl/nucl translated
            if(searchType == Parameters::SEARCH_TYPE_TRANSLATED||searchType == Parameters::SEARCH_TYPE_TRANS_NUCL_ALN){
                return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED| Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
            }else if (searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES ){
                return Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE| Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE;
            } else {
                Debug(Debug::ERROR) << "--search-type 1 (amino acid) can not used in combination with a nucleotide database\n "
                                       "The only possible options --search-types 2 (translated), 3 (nucleotide) or 4 (translated nucleotide backtrace)\n";
                EXIT(EXIT_FAILURE);
            }
        }
        // protein/protein
        if (Parameters::isEqualDbtype(queryDbType,  Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // profile/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // protein/profile
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
        }

        // profile/nucleotide
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // nucleotide/profile
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
        }

        // nucleotide/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // protein/nucleotide
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
        }
    } else{

        // protein/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_AMINO_ACIDS)){
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // profile/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_AMINO_ACIDS)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // protein/profile
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_HMM_PROFILE)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
        }

        // profile/nucleotide
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE | Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
        }

        // nucleotide/profile
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_HMM_PROFILE) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_HMM_PROFILE)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE;
        }

        // nucleotide/protein
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            (Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_AMINO_ACIDS))) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_AMINOACID;
        }

        // protein/nucleotide
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            (Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_NUCLEOTIDES))) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_AMINOACID | Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
        }

        // nucl/nucl
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE | Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE;
        }

        // nucl/nucl translated
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES) &&
            Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_AMINO_ACIDS) &&
            Parameters::isEqualDbtype(targetSrcDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
            return Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED | Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED;
        }

    }
    Debug(Debug::ERROR) << "Invalid input database and --search-type combination\n"
                        << "queryDbType: " << Parameters::getDbTypeName(queryDbType) << "\n"
                        << "targetDbType: " <<  Parameters::getDbTypeName(targetDbType) << "\n"
                        << "targetSrcDbType: " <<  Parameters::getDbTypeName(targetSrcDbType) << "\n"
                        << "searchMode: " << searchType << "\n";
    EXIT(EXIT_FAILURE);
}



void setNuclSearchDefaults(Parameters *p) {
    // leave ungapped alignment untouched
    if(p->alignmentMode != Parameters::ALIGNMENT_MODE_UNGAPPED){
        p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    //p->orfLongest = true;
    p->exactKmerMatching = true;
//    if ( p->PARAM_DIAGONAL_SCORING.wasSet == false) {
//        p->diagonalScoring = 0;
//    }
    if ( p->PARAM_STRAND.wasSet == false) {
        p->strand = 2;
    }
    if ( p->PARAM_K.wasSet == false) {
        p->kmerSize = 15;
    }
    if (  p->PARAM_MAX_SEQ_LEN.wasSet == false) {
        p->maxSeqLen = 10000;
    }
}


int search(mmseqs_output* out, Parameters &par) {
//    Parameters &par = Parameters::getInstance();
//    setSearchDefaults(&par);
//    par.PARAM_COV_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_C.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_MIN_SEQ_ID.addCategory(MMseqsParameter::COMMAND_EXPERT);
//    for (size_t i = 0; i < par.extractorfs.size(); i++) {
//        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
//    }
//    for (size_t i = 0; i < par.translatenucs.size(); i++) {
//        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
//    }
//    for (size_t i = 0; i < par.splitsequence.size(); i++) {
//        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
//    }
//    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
//    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);
//
//    par.parseParameters(argc, argv, command, false, 0, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);

    std::string indexStr = PrefilteringIndexReader::searchForIndex(par.db2);

    int targetDbType = FileUtil::parseDbType(par.db2.c_str());
    std::string targetDB =  (indexStr == "") ? par.db2.c_str() : indexStr.c_str();
    int targetSrcDbType = -1;
    if(indexStr != "" || Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_INDEX_DB)){
        indexStr = par.db2;
        DBReader<unsigned int> dbr(targetDB.c_str(), (targetDB+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr.open(DBReader<unsigned int>::NOSORT);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&dbr);
        targetSrcDbType = data.srcSeqType;
        targetDbType = data.seqType;
    }
    const int queryDbType = FileUtil::parseDbType(par.db1.c_str());
    if (queryDbType == -1 || targetDbType == -1) {
        Debug(Debug::ERROR)
                << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        EXIT(EXIT_FAILURE);
    }

    int searchMode = computeSearchMode(queryDbType, targetDbType, targetSrcDbType, par.searchType);
    if ((searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE) && (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE)) {
        setNuclSearchDefaults(&par);
    } else{
        par.PARAM_STRAND.addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    // FIXME: use larger default k-mer size in target-profile case if memory is available
    // overwrite default kmerSize for target-profile searches and parse parameters again
    if (par.sliceSearch == false && (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) && par.PARAM_K.wasSet == false) {
        par.kmerSize = 5;
    }

    const bool isUngappedMode = par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
    if (isUngappedMode && (searchMode & (Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE |Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE ))) {
        // par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with profile databases\n";
        EXIT(EXIT_FAILURE);
    }

    if (isUngappedMode && par.lcaSearch) {
        // par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with lca search\n";
        EXIT(EXIT_FAILURE);
    }

    // validate and set parameters for iterative search
    if (par.numIterations > 1) {
        if (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) {
            // par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
            Debug(Debug::ERROR) << "Iterative target-profile searches are not supported.\n";
            EXIT(EXIT_FAILURE);
        }

        par.addBacktrace = true;
        if (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_PROFILE) {
            for (size_t i = 0; i < par.searchworkflow.size(); i++) {
                if (par.searchworkflow[i]->uniqid == par.PARAM_REALIGN.uniqid && par.searchworkflow[i]->wasSet) {
                    // par.printUsageMessage(command,
                    //                      MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_PREFILTER);
                    Debug(Debug::ERROR) << "Cannot realign query profiles.\n";
                    EXIT(EXIT_FAILURE);
                }
            }

            par.realign = false;
        }
    }
    //par.printParameters(command.cmd, argc, argv, par.searchworkflow);

    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(par.databases_types, par.filenames, par.searchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    std::string parentTmpDir = tmpDir;
    tmpDir = FileUtil::createTemporaryDirectory(par.baseTmpPath, tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    const int originalRescoreMode = par.rescoreMode;
    bool needTargetSplit = false;
    bool extract_frames = false;
    bool need_query_split = false;
    std::string search_prog = "";
    
    CommandCaller cmd;
    out->output_string("VERBOSITY", par.createParameterString(par.onlyverbosity));
    out->output_string("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression));
    out->output_string("VERB_COMP_PAR", par.createParameterString(par.verbandcompression));
    if (isUngappedMode) {
        out->output_string("ALIGN_MODULE", "rescorediagonal");
    } else if (par.lcaSearch) {
        out->output_string("ALIGN_MODULE", "lcaalign");
    } else {
        out->output_string("ALIGN_MODULE", "align");
    }
    out->output_string("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : "");
    std::string program;
    out->output_string("RUNNER", par.runner);
//    out->output_string("ALIGNMENT_DB_EXT", Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_PROFILE_STATE_SEQ) ? ".255" : "");
    par.filenames[1] = targetDB;
    if (par.sliceSearch == true) {
        // By default (0), diskSpaceLimit (in bytes) will be set in the workflow to use as much as possible
        out->output_string("AVAIL_DISK", SSTR(static_cast<size_t>(par.diskSpaceLimit)));

        // correct Eval threshold for inverted search
        const size_t queryDbSize = FileUtil::countLines(par.db1Index.c_str());
        const size_t targetDbSize = FileUtil::countLines(par.db2Index.c_str());
        par.evalThr *= ((float) queryDbSize)/targetDbSize;

        int originalCovMode = par.covMode;
        par.covMode = Util::swapCoverageMode(par.covMode);
        size_t maxResListLen = par.maxResListLen;
        par.maxResListLen = INT_MAX;
        out->output_string("PREFILTER_PAR", par.createParameterString(par.prefilter));
        par.maxResListLen = maxResListLen;
        double originalEvalThr = par.evalThr;
        par.evalThr = std::numeric_limits<double>::max();
        out->output_string("SWAP_PAR", par.createParameterString(par.swapresult));
        par.evalThr = originalEvalThr;
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            out->output_string("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal));
            par.rescoreMode = originalRescoreMode;
        } else {
            out->output_string("ALIGNMENT_PAR", par.createParameterString(par.align));
        }
        out->output_string("SORTRESULT_PAR", par.createParameterString(par.sortresult));
        par.covMode = originalCovMode;

        program = tmpDir + "/searchslicedtargetprofile.sh";
        FileUtil::writeFile(program, searchslicedtargetprofile_sh, searchslicedtargetprofile_sh_len);
    } else if (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) {
        out->output_string("PREFILTER_PAR", par.createParameterString(par.prefilter));
        // we need to align all hits in case of target Profile hits
        size_t maxResListLen = par.maxResListLen;
        par.maxResListLen = INT_MAX;
        int originalCovMode = par.covMode;
        par.covMode = Util::swapCoverageMode(par.covMode);
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            out->output_string("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal));
            par.rescoreMode = originalRescoreMode;
        } else {
            out->output_string("ALIGNMENT_PAR", par.createParameterString(par.align));
        }
        par.covMode = originalCovMode;
        par.maxResListLen = maxResListLen;
        out->output_string("SWAP_PAR", par.createParameterString(par.swapresult));
        FileUtil::writeFile(tmpDir + "/searchtargetprofile.sh", searchtargetprofile_sh, searchtargetprofile_sh_len);
        program = std::string(tmpDir + "/searchtargetprofile.sh");
    } else if (par.numIterations > 1) {
        out->output_string("NUM_IT", SSTR(par.numIterations));
        out->output_string("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs));
        out->output_string("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity));

        double originalEval = par.evalThr;
        par.evalThr = (par.evalThr < par.evalProfile) ? par.evalThr  : par.evalProfile;
        for (int i = 0; i < par.numIterations; i++) {
            if (i == 0 && (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) == false) {
                par.realign = true;
            }

            if (i > 0) {
//                par.queryProfile = true;
                par.realign = false;
            }

            if (i == (par.numIterations - 1)) {
                par.evalThr = originalEval;
            }

            out->output_string(std::string("PREFILTER_PAR_" + SSTR(i)),
                            par.createParameterString(par.prefilter));
            if (isUngappedMode) {
                par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
                out->output_string(std::string("ALIGNMENT_PAR_" + SSTR(i)),
                                par.createParameterString(par.rescorediagonal));
                par.rescoreMode = originalRescoreMode;
            } else {
                out->output_string(std::string("ALIGNMENT_PAR_" + SSTR(i)),
                                par.createParameterString(par.align));
            }
            par.pca = 0.0;
            out->output_string(std::string("PROFILE_PAR_" + SSTR(i)),
                            par.createParameterString(par.result2profile));
            par.pca = 1.0;
        }

        FileUtil::writeFile(tmpDir + "/blastpgp.sh", blastpgp_sh, blastpgp_sh_len);
        program = std::string(tmpDir + "/blastpgp.sh");
    } else {
        if (par.sensSteps > 1) {
            if (par.startSens > par.sensitivity) {
                Debug(Debug::ERROR) << "--start-sens should not be greater -s.\n";
                EXIT(EXIT_FAILURE);
            }
            out->output_string("SENSE_0", SSTR(par.startSens));
            float sensStepSize = (par.sensitivity - par.startSens) / (static_cast<float>(par.sensSteps) - 1);
            for (int step = 1; step < par.sensSteps; step++) {
                std::string stepKey = "SENSE_" + SSTR(step);
                float stepSense = par.startSens + sensStepSize * step;
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << stepSense;
                std::string value = stream.str();
                out->output_string(stepKey.c_str(), value);
            }
            out->output_string("STEPS", SSTR((int) par.sensSteps));
        } else {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(1) << par.sensitivity;
            std::string sens = stream.str();
            out->output_string("SENSE_0", sens);
            out->output_string("STEPS", SSTR(1));
        }

        std::vector<MMseqsParameter*> prefilterWithoutS;
        for (size_t i = 0; i < par.prefilter.size(); i++) {
            if (par.prefilter[i]->uniqid != par.PARAM_S.uniqid) {
                prefilterWithoutS.push_back(par.prefilter[i]);
            }
        }
        out->output_string("PREFILTER_PAR", par.createParameterString(prefilterWithoutS));
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            out->output_string("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal));
            par.rescoreMode = originalRescoreMode;
        } else {
            out->output_string("ALIGNMENT_PAR", par.createParameterString(par.align));
        }
        FileUtil::writeFile(tmpDir + "/blastp.sh", blastp_sh, blastp_sh_len);
        program = std::string(tmpDir + "/blastp.sh");
    }

    if (searchMode & (Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED|Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)) {
        out->output_string("NO_TARGET_INDEX", (indexStr == "") ? "TRUE" : "");
        FileUtil::writeFile(tmpDir + "/translated_search.sh", translated_search_sh, translated_search_sh_len);
        out->output_string("QUERY_NUCL", (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED) ? "TRUE" : "");
        out->output_string("TARGET_NUCL", (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)  ? "TRUE" : "");
        out->output_string("THREAD_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
        par.subDbMode = 1;
        out->output_string("CREATESUBDB_PAR", par.createParameterString(par.createsubdb));
        par.translate = 1;
        out->output_string("ORF_PAR", par.createParameterString(par.extractorfs));
        out->output_string("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment));
        out->output_string("SEARCH", program);
        search_prog = program;
        program = std::string(tmpDir + "/translated_search.sh");
    }else if(searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE &&
            searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE){
        FileUtil::writeFile(tmpDir + "/blastn.sh", blastn_sh, blastn_sh_len);
        //  0: reverse, 1: forward, 2: both
        switch (par.strand){
            case 0:
                par.forwardFrames= "";
                par.reverseFrames= "1";
                out->output_string("EXTRACTFRAMES","TRUE");
                extract_frames = true;
                break;
            case 1:
                par.forwardFrames= "1";
                par.reverseFrames= "";
                break;
            case 2:
                par.forwardFrames= "1";
                par.reverseFrames= "1";
                out->output_string("EXTRACTFRAMES","TRUE");
                extract_frames = true;
                break;
        }
        out->output_string("SPLITSEQUENCE_PAR", par.createParameterString(par.splitsequence));
        if(indexStr=="") {
            out->output_string("NEEDTARGETSPLIT", "TRUE");
            needTargetSplit = true;
        }
        out->output_string("NEEDQUERYSPLIT","TRUE");
        need_query_split = true;
        out->output_string("EXTRACT_FRAMES_PAR", par.createParameterString(par.extractframes));
        out->output_string("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment));
        out->output_string("SEARCH", program);
        search_prog = program;
        program = std::string(tmpDir + "/blastn.sh");

    }

    //////////

    std::string query = par.filenames[0];
    std::string target = par.filenames[1];
    std::string result = par.filenames[2];

    std::cout << "step_search A\n" << std::flush;
    if (needTargetSplit) {
        if (!FileUtil::fileExists((tmpDir + "/target_seqs_split.dbtype").c_str())) {
            std::cout << "step_search B\n" << std::flush;
            Parameters splitsequence_par;

// OFFSETALIGNMENT_PAR: [--chain-alignments 0 --merge-query 1 --search-type 0 --threads 1 --compressed 0 --db-load-mode 0 -v 3 ]
//EXTRACT_FRAMES_PAR: [--forward-frames 1,2,3 --reverse-frames 1,2,3 --create-lookup 0 --threads 1 --compressed 0 -v 3 ]
// SPLITSEQUENCE_PAR: [--max-seq-len 65535 --sequence-overlap 0 --sequence-split-mode 1 --headers-split-mode 0 --create-lookup 0 --threads 1 --compressed 0 -v 3 ]
            splitsequence_par.setDBFields(1, target);
            splitsequence_par.setDBFields(2, tmpDir + "/target_seqs_split");
            splitsequence_par.maxSeqLen = 65535;
            splitsequence_par.sequenceOverlap = 0;
            splitsequence_par.sequenceSplitMode = 1;
            splitsequence_par.headerSplitMode = 0;
            splitsequence_par.createLookup = 0;
            splitsequence_par.threads = 1;
            splitsequence_par.compressed = 0;
            subcall_mmseqs(out, "splitsequence", splitsequence_par);
        }
        target = tmpDir + "/target_seqs_split";
        std::cout << "step_search C\n" << std::flush;
    }

    std::cout << "step_search D\n" << std::flush;
    if (extract_frames) {
        if (!FileUtil::fileExists((tmpDir + "/query_seqs.dbtype").c_str())) {
            std::cout << "step_search E\n" << std::flush;
            Parameters extractframes_par;

//EXTRACT_FRAMES_PAR: [--forward-frames 1,2,3 --reverse-frames 1,2,3 --create-lookup 0 --threads 1 --compressed 0 -v 3 ]
            extractframes_par.setDBFields(1, query);
            extractframes_par.setDBFields(2, tmpDir + "/query_seqs");
            extractframes_par.forwardFrames = "1,2,3";
            extractframes_par.reverseFrames = "1,2,3";
            extractframes_par.createLookup = 0;
            extractframes_par.threads = 1;
            extractframes_par.compressed = 0;
            subcall_mmseqs(out, "extractframes", extractframes_par);
        }
        query = tmpDir + "/query_seqs";
        std::cout << "step_search F\n" << std::flush;
    }

    std::cout << "step_search G\n" << std::flush;
    if (need_query_split) {
        if (!FileUtil::fileExists((tmpDir + "/query_seqs_split.dbtype").c_str())) {
            std::cout << "step_search H\n" << std::flush;
            Parameters splitsequence_par;
            splitsequence_par.setDBFields(1, query);
            splitsequence_par.setDBFields(2, tmpDir + "/query_seqs_split");
            splitsequence_par.maxSeqLen = 65535;
            splitsequence_par.sequenceOverlap = 0;
            splitsequence_par.sequenceSplitMode = 1;
            splitsequence_par.headerSplitMode = 0;
            splitsequence_par.createLookup = 0;
            splitsequence_par.threads = 1;
            splitsequence_par.compressed = 0;
            subcall_mmseqs(out, "splitsequence", splitsequence_par);
        }
        query = tmpDir + "/query_seqs_split";
        std::cout << "step_search I\n" << std::flush;
    }
    
    // mkdir tmpDir
    std::string searchTmpDir = parentTmpDir + "/search";
    //FileUtil::makeDir(searchTmpDir.c_str());

    std::cout << "step_search J\n" << std::flush;
    if (!FileUtil::fileExists((tmpDir + "/aln.dbtype").c_str())) {
        Parameters nested_search_par;
        std::cout << "step_search K\n" << std::flush;
        std::vector<std::string> search_filenames = {
            query,
            target,
            tmpDir + "/aln",
            searchTmpDir,
        };
        nested_search_par.filenames = search_filenames;
        nested_search_par.setDBFields(1, query);
        nested_search_par.setDBFields(2, target);
        nested_search_par.setDBFields(3, tmpDir + "/aln");
        nested_search_par.setDBFields(4, searchTmpDir);
        nested_search_par.spacedKmer = true;
        nested_search_par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV;
        nested_search_par.sensitivity = 5.7;
        nested_search_par.evalThr = 0.001;
        nested_search_par.orfStartMode = 1;
        nested_search_par.orfMinLength = 30;
        nested_search_par.orfMaxLength = 32734;
        nested_search_par.evalProfile = 0.1;
        nested_search_par.baseTmpPath = par.baseTmpPath;


        Debug(Debug::ERROR) << "Call: " << program << "\n";
        out->print();
        EXIT(EXIT_FAILURE);

        subcall_mmseqs(out, "search", nested_search_par);
        std::cout << "step_search L\n" << std::flush;
    }

    std::cout << "step_search M\n" << std::flush;
    if (!FileUtil::fileExists((result+".dbtype").c_str())) {
        std::cout << "step_search N\n" << std::flush;
        Parameters offsetalignment_par;
        std::vector<std::string> alignment_filenames = {
            par.filenames[0],
            query,
            par.filenames[1],
            target,
            tmpDir + "/aln",
            result,
        };
        offsetalignment_par.filenames = alignment_filenames;
        offsetalignment_par.setDBFields(1, par.filenames[0]);
        offsetalignment_par.setDBFields(2, query);
        offsetalignment_par.setDBFields(3, par.filenames[1]);
        offsetalignment_par.setDBFields(4, target);
        offsetalignment_par.setDBFields(5, tmpDir + "/aln");
        offsetalignment_par.setDBFields(6, result);
        offsetalignment_par.chainAlignment = 0;
        offsetalignment_par.mergeQuery = 1;
        offsetalignment_par.searchType = 0;
        offsetalignment_par.threads = 1;
        offsetalignment_par.compressed = 0;
        offsetalignment_par.preloadMode = 0;
        subcall_mmseqs(out, "offsetalignment", offsetalignment_par);
        std::cout << "step_search O\n" << std::flush;
    }

    std::cout << "step_search P\n" << std::flush;
    if (par.removeTmpFiles) {
        std::cout << "step_search Q\n" << std::flush;
        Parameters rmdb_par;
        rmdb_par.setDBFields(1, tmpDir + "/q_orfs");
        subcall_mmseqs(out, "rmdb", rmdb_par);
        rmdb_par.setDBFields(1, tmpDir + "/q_orfs_aa");
        subcall_mmseqs(out, "rmdb", rmdb_par);
        rmdb_par.setDBFields(1, tmpDir + "/t_orfs");
        subcall_mmseqs(out, "rmdb", rmdb_par);
        rmdb_par.setDBFields(1, tmpDir + "/t_orfs_aa");
        subcall_mmseqs(out, "rmdb", rmdb_par);
    }
    std::cout << "step_search R\n" << std::flush;

    //cmd.execProgram(program.c_str(), par.filenames);
    // Should never get here
    //assert(false);

    return 0;
}
