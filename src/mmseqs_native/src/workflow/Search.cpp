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
#include <set>

/*
std::string blastp_input = query;
std::string blastp_target = target;
std::string blastp_tmp = searchTmpDir;
std::string blastp_out = tmpDir + "/aln";
*/
void call_blastp(mmseqs_output* out, Parameters &par, int no_steps, std::vector<float> senses, std::string align_module, std::string blastp_input, std::string blastp_target, std::string blastp_out, std::string blastp_tmp, std::string script_path, CommandCaller cmd) {
    if (false) {
        std::cout << "[blastn.sh] Executing native script: " << script_path << "\n" << std::flush;
        std::vector<std::string> args = {
            blastp_input, blastp_target, blastp_out, blastp_tmp
        };
        cmd.addVar("MMSEQS", "mmseqs");
        cmd.callProgram(script_path.c_str(), args);
        std::cout << "[blastn.sh] Exit\n" << std::flush;
    } else {

        std::cout << "step_search L\n" << std::flush;

        // "$SEARCH" "${QUERY}" "${TARGET}" "$4/aln" "$4/search"
        int blastp_steps = no_steps;
        int current_step = 0;

        std::set<std::string> merged_senses;

        std::string aln_res_merge = blastp_tmp + "/aln_0";
        while (current_step < blastp_steps) {
            std::cout << "step_search K_1\n" << std::flush;
            const float sens = senses[current_step];
            std::string sens_str = std::to_string(((float)((int)(sens*10)))/10);
            // Call prefilter module
            std::string step_str = std::to_string(current_step);
            if (!FileUtil::fileExists((blastp_tmp + "/pref_" + step_str + ".dbtype").c_str())) {
                std::cout << "step_search K_2: prefilter " << blastp_input << blastp_target << (blastp_tmp + "/pref_" + step_str) << "\n" << std::flush;
                Parameters prefilter_par(par);
                std::vector<std::string> prefilter_filenames = {
                    blastp_input, blastp_target, (blastp_tmp + "/pref_" + step_str)
                };
                prefilter_par.filenames = prefilter_filenames;
                prefilter_par.setDBFields(1, blastp_input);
                prefilter_par.setDBFields(2, blastp_target);
                prefilter_par.setDBFields(3, blastp_tmp + "/pref_" + step_str);
                prefilter_par.setSubstitutionMatrices("blosum62.out", "nucleotide.out");
                prefilter_par.setSeedSubstitutionMatrices("VTML80.out", "nucleotide.out");
                prefilter_par.kmerSize = 0;
                prefilter_par.kmerScore = 2147483647;
                prefilter_par.alphabetSize = MultiParam<int>(21, 5);
                prefilter_par.maxSeqLen = 65535;
                prefilter_par.maxResListLen = 300;
                prefilter_par.split = 0;
                prefilter_par.splitMode = 2;
                prefilter_par.splitMemoryLimit = 0;
                prefilter_par.covThr = 0;
                prefilter_par.covMode = 0;
                prefilter_par.compBiasCorrection = 1;
                prefilter_par.diagonalScoring = 1;
                prefilter_par.exactKmerMatching = 0;
                prefilter_par.maskMode = 1;
                prefilter_par.maskLowerCaseMode = 0;
                prefilter_par.minDiagScoreThr = 15;
                prefilter_par.includeIdentity = 0;
                prefilter_par.spacedKmer = 1;
                prefilter_par.preloadMode = 0;
                prefilter_par.pca = 1;
                prefilter_par.pcb = 1.5;
                prefilter_par.threads = 1;
                prefilter_par.compressed = 0;
                prefilter_par.sensitivity = sens;
                subcall_mmseqs(out, "prefilter", prefilter_par);
                std::cout << "step_search K_3\n" << std::flush;
            }

            // call alignment module
            std::string align_path = blastp_out;
            if (blastp_steps > 1) {
                align_path = blastp_tmp + "/aln_" + step_str;
            }
            std::cout << "step_search K_4\n" << std::flush;
            if (!FileUtil::fileExists((align_path + ".dbtype").c_str())) {
                std::cout << "step_search K_5\n" << std::flush;
                Parameters align_par(par);
                std::vector<std::string> align_filenames = {
                    blastp_input, blastp_target, (blastp_tmp + "/pref_" + step_str), align_path,
                };
                align_par.setDBFields(1, blastp_input);
                align_par.setDBFields(2, blastp_target);
                align_par.setDBFields(3, blastp_tmp + "/pref_" + step_str);
                align_par.setDBFields(4, align_path);
                align_par.filenames = align_filenames;
                align_par.threads = 1;
                align_par.compressed = 0;
                align_par.setSubstitutionMatrices("blosum62.out", "nucleotide.out");
                align_par.addBacktrace = 0;
                align_par.alignmentMode = 3;
                align_par.wrappedScoring = 0;
                align_par.evalThr = 0.001;
                align_par.seqIdThr = 0;
                align_par.alnLenThr = 0;
                align_par.seqIdMode = 0;
                align_par.altAlignment = 0;
                align_par.covThr = 0;
                align_par.covMode = 0;
                align_par.maxSeqLen = 65535;
                align_par.compBiasCorrection = 1;
                align_par.maxRejected = 2147483647;
                align_par.maxAccept = 2147483647;
                align_par.includeIdentity = 0;
                align_par.preloadMode = 0;
                align_par.pca = 1;
                align_par.pcb = 1.5;
                align_par.scoreBias = 0;
                align_par.realign = 0;
                align_par.realignScoreBias = -0.2;
                align_par.realignMaxSeqs = 2147483647;
                align_par.gapOpen = MultiParam<int>(11, 5);
                align_par.gapExtend = MultiParam<int>(1, 2);
                align_par.zdrop = 40;

                subcall_mmseqs(out, align_module, align_par);
                std::cout << "step_search K_6\n" << std::flush;
            }

            // Only merge results after first step
            std::cout << "step_search K_7\n" << std::flush;
            if (current_step > 0) {
                if (merged_senses.find(sens_str) == merged_senses.end()) {
                    if (current_step < blastp_steps-1) {
                        std::cout << "step_search K_8\n" << std::flush;
                        Parameters mergedbs_par(par);
                        std::vector<std::string> mergedbs_filenames = {
                            blastp_input, (blastp_tmp + "/aln_merge_new"), (blastp_tmp + "/pref_" + step_str), align_path,
                        };
                        mergedbs_par.filenames = mergedbs_filenames;
                        mergedbs_par.setDBFields(1, blastp_input);
                        mergedbs_par.setDBFields(2, blastp_tmp + "/aln_merge_new");
                        mergedbs_par.setDBFields(3, blastp_tmp + "/pref_" + step_str);
                        mergedbs_par.setDBFields(4, align_path);
                        mergedbs_par.compressed = 0;
                        subcall_mmseqs(out, "mergedbs", mergedbs_par);

                        Parameters rmdb_par(par);
                        rmdb_par.setDBFields(1, blastp_tmp + "/aln_merge");
                        subcall_mmseqs(out, "rmdb", rmdb_par);

                        Parameters mvdb_par(par);
                        mvdb_par.setDBFields(1, blastp_tmp + "/aln_merge_new");
                        mvdb_par.setDBFields(2, blastp_tmp + "/aln_merge");
                        subcall_mmseqs(out, "mvdb", mvdb_par);
                        std::cout << "step_search K_9\n" << std::flush;
                    } else {
                        std::cout << "step_search K_10\n" << std::flush;
                        Parameters mergedbs_par(par);
                        std::vector<std::string> mergedbs_filenames = {
                            blastp_input, align_path, (blastp_tmp + "/pref_" + step_str), align_path
                        };
                        mergedbs_par.filenames = mergedbs_filenames;
                        mergedbs_par.setDBFields(1, blastp_input);
                        mergedbs_par.setDBFields(2, align_path);
                        mergedbs_par.setDBFields(3, blastp_tmp + "/pref_" + step_str);
                        mergedbs_par.setDBFields(4, align_path);
                        mergedbs_par.compressed = 0;
                        subcall_mmseqs(out, "mergedbs", mergedbs_par);
                        std::cout << "step_search K_11\n" << std::flush;
                    }
                    merged_senses.insert(sens_str);
                }
                aln_res_merge = blastp_tmp + "/aln_merge";
            }
            std::cout << "step_search K_12\n" << std::flush;

            std::string next_input = blastp_tmp + "/input_" + step_str;
            // do not create subdb at last step
            if (current_step < blastp_steps-1) {
                std::cout << "step_search K_13\n" << std::flush;
                if (!FileUtil::fileExists((blastp_tmp + "/order_" + step_str + ".dbtype").c_str())) {
                    // awk '$3 < 2 { print $1 }' "$TMP_PATH/aln_$STEP.index" > "$TMP_PATH/order_$STEP"
                    std::cout << "SPIERDALAJ?\n";
                    std::cout << "ELO\n";
                    EXIT(EXIT_FAILURE);
                }
                std::cout << "step_search K_14\n" << std::flush;

                if (!FileUtil::fileExistsAndIsNotEmpty((blastp_tmp + "/order_" + step_str).c_str())) {
                   // "$MMSEQS" mvdb "$ALN_RES_MERGE" "$3" ${VERBOSITY}
                    Parameters mvdb_par(par);
                    mvdb_par.setDBFields(1, aln_res_merge);
                    mvdb_par.setDBFields(2, align_path);
                    subcall_mmseqs(out, "mvdb", mvdb_par);
                }
                std::cout << "step_search K_15\n" << std::flush;

                if (!FileUtil::fileExists((next_input + ".dbtype").c_str())) {
                    Parameters createsubdb_par(par);
                    std::vector<std::string> createsubdb_filenames = {
                        (blastp_tmp + "/order_" + step_str),
                        blastp_input,
                        next_input,
                    };
                    createsubdb_par.filenames = createsubdb_filenames;
                    createsubdb_par.setDBFields(1, blastp_tmp + "/order_" + step_str);
                    createsubdb_par.setDBFields(2, blastp_input);
                    createsubdb_par.setDBFields(3, next_input);
                    createsubdb_par.subDbMode = 1;
                    subcall_mmseqs(out, "createsubdb", createsubdb_par);
                }
                std::cout << "step_search K_16\n" << std::flush;
            }

            blastp_input = next_input;
            current_step++;
            std::cout << "step_search K_17\n" << std::flush;
        }
        std::cout << "step_search K_18\n" << std::flush;

        if (par.removeTmpFiles) {
            for(int i=0;i<blastp_steps;++i) {
                std::string step_str = std::to_string(i);
                Parameters rmdb_par(par);
                rmdb_par.setDBFields(1, blastp_tmp + "/pref_" + step_str);
                subcall_mmseqs(out, "rmdb", rmdb_par);

                rmdb_par.setDBFields(1, blastp_tmp + "/aln_" + step_str);
                subcall_mmseqs(out, "rmdb", rmdb_par);

                rmdb_par.setDBFields(1, blastp_tmp + "/input_" + step_str);
                subcall_mmseqs(out, "rmdb", rmdb_par);
            }

            Parameters rmdb_par(par);
            rmdb_par.setDBFields(1, blastp_tmp + "/aln_merge");
            subcall_mmseqs(out, "rmdb", rmdb_par);
        }
        std::cout << "step_search K_19\n" << std::flush;
    }
}


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
    int no_steps = 0;
    std::string align_module = "";
    std::vector<float> senses;
    std::string search_prog = "";

    CommandCaller cmd;
    cmd.addVar("VERBOSITY", par.createParameterString(par.onlyverbosity));
    cmd.addVar("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression));
    cmd.addVar("VERB_COMP_PAR", par.createParameterString(par.verbandcompression));
    if (isUngappedMode) {
        cmd.addVar("ALIGN_MODULE", "rescorediagonal");
        align_module = "rescorediagonal";
    } else if (par.lcaSearch) {
        cmd.addVar("ALIGN_MODULE", "lcaalign");
        align_module = "lcaalign";
    } else {
        cmd.addVar("ALIGN_MODULE", "align");
        align_module = "align";
    }
    cmd.addVar("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : "");
    std::string program;
    cmd.addVar("RUNNER", par.runner);
//    cmd.addVar("ALIGNMENT_DB_EXT", Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_PROFILE_STATE_SEQ) ? ".255" : "");
    par.filenames[1] = targetDB;
    if (par.sliceSearch == true) {
        // By default (0), diskSpaceLimit (in bytes) will be set in the workflow to use as much as possible
        cmd.addVar("AVAIL_DISK", SSTR(static_cast<size_t>(par.diskSpaceLimit)));

        // correct Eval threshold for inverted search
        const size_t queryDbSize = FileUtil::countLines(par.db1Index.c_str());
        const size_t targetDbSize = FileUtil::countLines(par.db2Index.c_str());
        par.evalThr *= ((float) queryDbSize)/targetDbSize;

        int originalCovMode = par.covMode;
        par.covMode = Util::swapCoverageMode(par.covMode);
        size_t maxResListLen = par.maxResListLen;
        par.maxResListLen = INT_MAX;
        cmd.addVar("PREFILTER_PAR", par.createParameterString(par.prefilter));
        par.maxResListLen = maxResListLen;
        double originalEvalThr = par.evalThr;
        par.evalThr = std::numeric_limits<double>::max();
        cmd.addVar("SWAP_PAR", par.createParameterString(par.swapresult));
        par.evalThr = originalEvalThr;
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVar("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal));
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVar("ALIGNMENT_PAR", par.createParameterString(par.align));
        }
        cmd.addVar("SORTRESULT_PAR", par.createParameterString(par.sortresult));
        par.covMode = originalCovMode;

        program = tmpDir + "/searchslicedtargetprofile.sh";
        FileUtil::writeFile(program, searchslicedtargetprofile_sh, searchslicedtargetprofile_sh_len);
    } else if (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_PROFILE) {
        cmd.addVar("PREFILTER_PAR", par.createParameterString(par.prefilter));
        // we need to align all hits in case of target Profile hits
        size_t maxResListLen = par.maxResListLen;
        par.maxResListLen = INT_MAX;
        int originalCovMode = par.covMode;
        par.covMode = Util::swapCoverageMode(par.covMode);
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVar("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal));
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVar("ALIGNMENT_PAR", par.createParameterString(par.align));
        }
        par.covMode = originalCovMode;
        par.maxResListLen = maxResListLen;
        cmd.addVar("SWAP_PAR", par.createParameterString(par.swapresult));
        FileUtil::writeFile(tmpDir + "/searchtargetprofile.sh", searchtargetprofile_sh, searchtargetprofile_sh_len);
        program = std::string(tmpDir + "/searchtargetprofile.sh");
    } else if (par.numIterations > 1) {
        cmd.addVar("NUM_IT", SSTR(par.numIterations));
        cmd.addVar("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs));
        cmd.addVar("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity));

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

            cmd.addVar(std::string("PREFILTER_PAR_" + SSTR(i)),
                            par.createParameterString(par.prefilter));
            if (isUngappedMode) {
                par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
                cmd.addVar(std::string("ALIGNMENT_PAR_" + SSTR(i)),
                                par.createParameterString(par.rescorediagonal));
                par.rescoreMode = originalRescoreMode;
            } else {
                cmd.addVar(std::string("ALIGNMENT_PAR_" + SSTR(i)),
                                par.createParameterString(par.align));
            }
            par.pca = 0.0;
            cmd.addVar(std::string("PROFILE_PAR_" + SSTR(i)),
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
            cmd.addVar("SENSE_0", SSTR(par.startSens));
            senses.push_back(par.startSens);
            float sensStepSize = (par.sensitivity - par.startSens) / (static_cast<float>(par.sensSteps) - 1);
            for (int step = 1; step < par.sensSteps; step++) {
                std::string stepKey = "SENSE_" + SSTR(step);
                float stepSense = par.startSens + sensStepSize * step;
                senses.push_back(stepSense);
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << stepSense;
                std::string value = stream.str();
                cmd.addVar(stepKey.c_str(), value);
            }
            no_steps = par.sensSteps;
            cmd.addVar("STEPS", SSTR((int) par.sensSteps));
        } else {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(1) << par.sensitivity;
            std::string sens = stream.str();
            cmd.addVar("SENSE_0", sens);
            cmd.addVar("STEPS", SSTR(1));
            senses.push_back(par.sensitivity);
            no_steps = 1;
        }

        std::vector<MMseqsParameter*> prefilterWithoutS;
        for (size_t i = 0; i < par.prefilter.size(); i++) {
            if (par.prefilter[i]->uniqid != par.PARAM_S.uniqid) {
                prefilterWithoutS.push_back(par.prefilter[i]);
            }
        }
        cmd.addVar("PREFILTER_PAR", par.createParameterString(prefilterWithoutS));
        if (isUngappedMode) {
            par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
            cmd.addVar("ALIGNMENT_PAR", par.createParameterString(par.rescorediagonal));
            par.rescoreMode = originalRescoreMode;
        } else {
            cmd.addVar("ALIGNMENT_PAR", par.createParameterString(par.align));
        }
        FileUtil::writeFile(tmpDir + "/blastp.sh", blastp_sh, blastp_sh_len);
        program = std::string(tmpDir + "/blastp.sh");
    }

    if (searchMode & (Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED|Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)) {
        cmd.addVar("NO_TARGET_INDEX", (indexStr == "") ? "TRUE" : "");
        FileUtil::writeFile(tmpDir + "/translated_search.sh", translated_search_sh, translated_search_sh_len);
        cmd.addVar("QUERY_NUCL", (searchMode & Parameters::SEARCH_MODE_FLAG_QUERY_TRANSLATED) ? "TRUE" : "");
        cmd.addVar("TARGET_NUCL", (searchMode & Parameters::SEARCH_MODE_FLAG_TARGET_TRANSLATED)  ? "TRUE" : "");
        cmd.addVar("THREAD_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
        par.subDbMode = 1;
        cmd.addVar("CREATESUBDB_PAR", par.createParameterString(par.createsubdb));
        par.translate = 1;
        cmd.addVar("ORF_PAR", par.createParameterString(par.extractorfs));
        cmd.addVar("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment));
        cmd.addVar("SEARCH", program);
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
                cmd.addVar("EXTRACTFRAMES","TRUE");
                extract_frames = true;
                break;
            case 1:
                par.forwardFrames= "1";
                par.reverseFrames= "";
                break;
            case 2:
                par.forwardFrames= "1";
                par.reverseFrames= "1";
                cmd.addVar("EXTRACTFRAMES","TRUE");
                extract_frames = true;
                break;
        }
        cmd.addVar("SPLITSEQUENCE_PAR", par.createParameterString(par.splitsequence));
        if(indexStr=="") {
            cmd.addVar("NEEDTARGETSPLIT", "TRUE");
            needTargetSplit = true;
        }
        cmd.addVar("NEEDQUERYSPLIT","TRUE");
        need_query_split = true;
        cmd.addVar("EXTRACT_FRAMES_PAR", par.createParameterString(par.extractframes));
        cmd.addVar("OFFSETALIGNMENT_PAR", par.createParameterString(par.offsetalignment));
        cmd.addVar("SEARCH", program);
        search_prog = program;
        program = std::string(tmpDir + "/blastn.sh");

    }

    std::string query = par.filenames[0];
    std::string target = par.filenames[1];
    std::string result = par.filenames[2];
    std::string searchTmpDir = par.baseTmpPath + parentTmpDir + "/search";
    // mkdir tmpDir
    FileUtil::makeDir(searchTmpDir.c_str());

    if (program == tmpDir + "/blastp.sh") {
        call_blastp(out, par, no_steps, senses, align_module, query, target, result, searchTmpDir, tmpDir + "/blastp.sh", cmd);
        return 0;
    }

    //////////

    std::cout << "step_search A\n" << std::flush;
    if (needTargetSplit) {
        if (!FileUtil::fileExists((tmpDir + "/target_seqs_split.dbtype").c_str())) {
            std::cout << "step_search B\n" << std::flush;
            Parameters splitsequence_par(par);

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
            Parameters extractframes_par(par);
            std::vector<std::string> extractframes_filenames = {
                query,
                tmpDir + "/query_seqs"
            };
            extractframes_par.filenames = extractframes_filenames;
//EXTRACT_FRAMES_PAR: [--forward-frames 1,2,3 --reverse-frames 1,2,3 --create-lookup 0 --threads 1 --compressed 0 -v 3 ]
            extractframes_par.setDBFields(1, query);
            extractframes_par.setDBFields(2, tmpDir + "/query_seqs");
            extractframes_par.forwardFrames = par.forwardFrames;
            extractframes_par.reverseFrames = par.reverseFrames;
            extractframes_par.createLookup = 0;
            extractframes_par.threads = 1;
            extractframes_par.compressed = 0;

            subcall_mmseqs(out, "extractframes", extractframes_par);
            std::cout << "CALL EXTRACT FRAMES {" << query << "} {" << (tmpDir + "/query_seqs") << "}\n" << std::flush;
            out->print();
            //EXIT(EXIT_FAILURE);
        }
        query = tmpDir + "/query_seqs";
        std::cout << "step_search F\n" << std::flush;
    }

    std::cout << "step_search G\n" << std::flush;
    if (need_query_split) {
        if (!FileUtil::fileExists((tmpDir + "/query_seqs_split.dbtype").c_str())) {
            std::cout << "step_search H " << query << " " << (tmpDir + "/query_seqs_split") << "\n" << std::flush;
            Parameters splitsequence_par(par);
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

//            out->print();
//            EXIT(EXIT_FAILURE);
        }
        query = tmpDir + "/query_seqs_split";
        std::cout << "step_search I\n" << std::flush;
    }


    std::cout << "step_search J\n" << std::flush;
    if (!FileUtil::fileExists((tmpDir + "/aln.dbtype").c_str())) {
        call_blastp(out, par, no_steps, senses, align_module, query, target, tmpDir + "/aln", searchTmpDir, "", cmd);
    }

    std::cout << "step_search M\n" << std::flush;
    if (!FileUtil::fileExists((result+".dbtype").c_str())) {


        Parameters offsetalignment_par(par);
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
        offsetalignment_par.searchType = par.searchType;
        offsetalignment_par.baseTmpPath = par.baseTmpPath;
//        offsetalignment_par.chainAlignment = 0;
//        offsetalignment_par.mergeQuery = 1;
//        //offsetalignment_par.searchType = 0;
//        offsetalignment_par.threads = 1;
//        offsetalignment_par.compressed = 0;
//        offsetalignment_par.preloadMode = 0;
        subcall_mmseqs(out, "offsetalignment", offsetalignment_par);
        out->print();
        std::cout << "mmseqs offsetalignment " << par.filenames[0] << " ";
        std::cout << query << " ";
        std::cout << par.filenames[1] << " ";
        std::cout << target << " ";
        std::cout << tmpDir + "/aln" << " ";
        std::cout << result << "\n";
        std::cout << "step_search O\n" << std::flush;
    }

    std::cout << "step_search P\n" << std::flush;
    if (par.removeTmpFiles) {
        std::cout << "step_search Q\n" << std::flush;
        Parameters rmdb_par(par);
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
