#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <biosnake/commons/application.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/log.h>
#include <biosnake/output.h>
#include <biosnake/api/api.h>

namespace py = pybind11;

biosnake_output call_biosnake_proxy(std::string command_name, Parameters args) {
  // pybind11::gil_scoped_release release;
  return call_biosnake(command_name, args);
}

PYBIND11_MODULE(biosnake_native, m) {
  m.doc() = "Main module of the biosnake";  // optional

    // make a new custom exception and use it as a translation target
    static py::exception<FatalException> ex(m, "MMSEQException");
    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const FatalException &e) {
            ex(e.what());
        }
    });

  pybind11::class_<SearchResults>(m, "SearchResults")
      .def(pybind11::init<>())
      .def_readonly("_records", &SearchResults::_records)
      .def_readonly("_headers", &SearchResults::_headers)
      .def_property_readonly("records", &SearchResults::getRecords);

  pybind11::class_<Database>(m, "Database")
      .def(pybind11::init<>())
      .def("remove", &Database::remove)
      .def("to_fasta", &Database::to_fasta,
            py::arg("output_path") = "")
      .def("copy", &Database::copy,
            py::arg("search_input_fasta") = "")
      .def("search", &Database::search,
            py::arg("sequences"),
            py::arg("search_type") = DEFAULT_SEARCH_TYPE,
            py::arg("headers") = std::vector<std::string>(),
            py::arg("sensitivity") = DEFAULT_SENSITIVITY,
            py::arg("max_sequence_length") = DEFAULT_MAX_SEQUENCE_LENGTH,
            py::arg("max_results_count_per_query") = DEFAULT_MAX_RESULTS_COUNT_PER_QUERY,
            py::arg("max_orf_length") = DEFAULT_MAX_ORF_LENGTH,
            py::arg("min_orf_length") = DEFAULT_MIN_ORF_LENGTH,
            py::arg("search_steps") = DEFAULT_SEARCH_STEPS,
            py::arg("start_sensitivity") = DEFAULT_START_SENSITIVITY
            )
      .def("search_file", &Database::search_file,
            py::arg("search_input_fasta"),
            py::arg("search_type") = DEFAULT_SEARCH_TYPE,
            py::arg("headers") = std::vector<std::string>(),
            py::arg("sensitivity") = DEFAULT_SENSITIVITY,
            py::arg("max_sequence_length") = DEFAULT_MAX_SEQUENCE_LENGTH,
            py::arg("max_results_count_per_query") = DEFAULT_MAX_RESULTS_COUNT_PER_QUERY,
            py::arg("max_orf_length") = DEFAULT_MAX_ORF_LENGTH,
            py::arg("min_orf_length") = DEFAULT_MIN_ORF_LENGTH,
            py::arg("search_steps") = DEFAULT_SEARCH_STEPS,
            py::arg("start_sensitivity") = DEFAULT_START_SENSITIVITY
            )
      .def("create_index", &Database::create_index,
            py::arg("search_type") = DEFAULT_SEARCH_TYPE,
            py::arg("sensitivity") = DEFAULT_SENSITIVITY,
            py::arg("max_sequence_length") = DEFAULT_MAX_SEQUENCE_LENGTH,
            py::arg("max_orf_length") = DEFAULT_MAX_ORF_LENGTH,
            py::arg("min_orf_length") = DEFAULT_MIN_ORF_LENGTH,
            py::arg("orf_start_mode") = DEFAULT_ORF_START_MODE
            )
      .def_property("name", &Database::getName, &Database::setName)
      .def_property("description", &Database::getDescription, &Database::setDescription)
      .def_property_readonly("type", &Database::getType)
      .def_property_readonly("columns_data", &Database::getColumnData)
      .def_property_readonly("size", &Database::getSize)
      .def("__iter__", [](Database &db) { return py::make_iterator(db.begin(), db.end()); }, py::keep_alive<0, 1>());

  pybind11::class_<Database::Record>(m, "Record")
      .def(pybind11::init<>())
      .def_readwrite("seq", &Database::Record::_sequence)
      .def_readwrite("id", &Database::Record::_header);

  pybind11::class_<Client>(m, "Client")
      .def(pybind11::init<const std::string&, const std::string&>())
      .def("create", &Client::create,
            py::arg("name"),
            py::arg("description"),
            py::arg("input_fasta"),
            py::arg("mode") = "copy",
            py::arg("database_type") = "auto",
            py::arg("offset") = 0,
            py::arg("shuffle") = false)
      .def("list", &Client::list)
      .def("get", &Client::get)
      .def("__getitem__", &Client::get)
      .def("__iter__", [](Client &c) { return py::make_iterator(c.begin(), c.end()); }, py::keep_alive<0, 1>());

  pybind11::class_<biosnake_call_args>(m, "BiosnakeCallArgs")
      .def(pybind11::init<>())
      .def_readwrite("cli_args", &biosnake_call_args::cli_args);

  pybind11::class_<biosnake_output>(m, "BiosnakeCallOutput")
      .def(pybind11::init<>())
      .def_readwrite("vars_str", &biosnake_output::vars_str)
      .def_readwrite("blast_tab_records", &biosnake_output::blast_tab_records);

  pybind11::class_<MultiParam<char*>>(m, "BiosnakeMultiParamString")
      .def(pybind11::init<>())
      .def_readwrite("aminoacids", &MultiParam<char*>::aminoacids)
      .def_readwrite("nucleotides", &MultiParam<char*>::nucleotides);

  pybind11::class_<MultiParam<int>>(m, "BiosnakeMultiParamInt")
      .def(pybind11::init<>())
      .def_readwrite("aminoacids", &MultiParam<int>::aminoacids)
      .def_readwrite("nucleotides", &MultiParam<int>::nucleotides);

  pybind11::class_<biosnake_blast_tab_record>(m, "BiosnakeSearchRecord")
      .def(pybind11::init<>())
      .def_readwrite("query_sequence_id",
                     &biosnake_blast_tab_record::query_sequence_id)
      .def_readwrite("target_sequence_id",
                     &biosnake_blast_tab_record::target_sequence_id)
      .def_readwrite("sequence_identity",
                     &biosnake_blast_tab_record::sequence_identity)
      .def_readwrite("alignment_length",
                     &biosnake_blast_tab_record::alignment_length)
      .def_readwrite("number_of_mismatches",
                     &biosnake_blast_tab_record::number_of_mismatches)
      .def_readwrite("number_of_gap_openings",
                     &biosnake_blast_tab_record::number_of_gap_openings)
      .def_readwrite("domain_start_index_query",
                     &biosnake_blast_tab_record::domain_start_index_query)
      .def_readwrite("domain_end_index_query",
                     &biosnake_blast_tab_record::domain_end_index_query)
      .def_readwrite("domain_start_index_target",
                     &biosnake_blast_tab_record::domain_start_index_target)
      .def_readwrite("domain_end_index_target",
                     &biosnake_blast_tab_record::domain_end_index_target)
      .def_readwrite("query_sequence_aligned",
                     &biosnake_blast_tab_record::query_sequence_aligned)
      .def_readwrite("target_sequence_aligned",
                     &biosnake_blast_tab_record::target_sequence_aligned)
      .def_readwrite("target_sequence_content",
                     &biosnake_blast_tab_record::target_sequence_content)
      .def_readwrite("query_sequence_content",
                     &biosnake_blast_tab_record::query_sequence_content)
      .def_readwrite("target_sequence_length",
                     &biosnake_blast_tab_record::target_sequence_length)
      .def_readwrite("query_sequence_length",
                     &biosnake_blast_tab_record::query_sequence_length)
      .def_readwrite("e_value", &biosnake_blast_tab_record::e_value)
      .def_readwrite("bit_score", &biosnake_blast_tab_record::bit_score)
      .def("__repr__", &biosnake_blast_tab_record::toString);

  pybind11::class_<Parameters>(m, "BiosnakeCallConfig")
      .def(pybind11::init<>())
      .def_readwrite("logFilePath", &Parameters::logFilePath)
      .def_readwrite("baseTmpPath", &Parameters::baseTmpPath)
      .def_readwrite("db1", &Parameters::db1)
      .def_readwrite("db1Index", &Parameters::db1Index)
      .def_readwrite("db1dbtype", &Parameters::db1dbtype)
      .def_readwrite("hdr1", &Parameters::hdr1)
      .def_readwrite("hdr1Index", &Parameters::hdr1Index)
      .def_readwrite("hdr1dbtype", &Parameters::hdr1dbtype)
      .def_readwrite("db2", &Parameters::db2)
      .def_readwrite("db2Index", &Parameters::db2Index)
      .def_readwrite("db2dbtype", &Parameters::db2dbtype)
      .def_readwrite("hdr2", &Parameters::hdr2)
      .def_readwrite("hdr2Index", &Parameters::hdr2Index)
      .def_readwrite("hdr2dbtype", &Parameters::hdr2dbtype)
      .def_readwrite("db3", &Parameters::db3)
      .def_readwrite("db3Index", &Parameters::db3Index)
      .def_readwrite("db3dbtype", &Parameters::db3dbtype)
      .def_readwrite("hdr3", &Parameters::hdr3)
      .def_readwrite("hdr3Index", &Parameters::hdr3Index)
      .def_readwrite("hdr3dbtype", &Parameters::hdr3dbtype)
      .def_readwrite("db4", &Parameters::db4)
      .def_readwrite("db4Index", &Parameters::db4Index)
      .def_readwrite("db4dbtype", &Parameters::db4dbtype)
      .def_readwrite("hdr4", &Parameters::hdr4)
      .def_readwrite("hdr4Index", &Parameters::hdr4Index)
      .def_readwrite("hdr4dbtype", &Parameters::hdr4dbtype)
      .def_readwrite("db5", &Parameters::db5)
      .def_readwrite("db5Index", &Parameters::db5Index)
      .def_readwrite("db5dbtype", &Parameters::db5dbtype)
      .def_readwrite("hdr5", &Parameters::hdr5)
      .def_readwrite("hdr5Index", &Parameters::hdr5Index)
      .def_readwrite("hdr5dbtype", &Parameters::hdr5dbtype)
      .def_readwrite("db6", &Parameters::db6)
      .def_readwrite("db6Index", &Parameters::db6Index)
      .def_readwrite("db6dbtype", &Parameters::db6dbtype)
      .def_readwrite("hdr6", &Parameters::hdr6)
      .def_readwrite("hdr6Index", &Parameters::hdr6Index)
      .def_readwrite("hdr6dbtype", &Parameters::hdr6dbtype)
      .def_readwrite("filenames", &Parameters::filenames)
      //.def_readwrite("restArgv", &Parameters::restArgv)
      //.def_readwrite("restArgc", &Parameters::restArgc)
      .def_readwrite("scoringMatrixFile", &Parameters::scoringMatrixFile)
      .def_readwrite("seedScoringMatrixFile",
                     &Parameters::seedScoringMatrixFile)
      .def_readwrite("maxSeqLen", &Parameters::maxSeqLen)
      .def_readwrite("maxResListLen", &Parameters::maxResListLen)
      .def_readwrite("verbosity", &Parameters::verbosity)
      .def_readwrite("threads", &Parameters::threads)
      .def_readwrite("compressed", &Parameters::compressed)
      .def_readwrite("removeTmpFiles", &Parameters::removeTmpFiles)
      .def_readwrite("includeIdentity", &Parameters::includeIdentity)
      .def_readwrite("sensitivity", &Parameters::sensitivity)
      .def_readwrite("kmerSize", &Parameters::kmerSize)
      .def_readwrite("kmerScore", &Parameters::kmerScore)
      .def_readwrite("alphabetSize", &Parameters::alphabetSize)
      .def_readwrite("compBiasCorrection", &Parameters::compBiasCorrection)
      .def_readwrite("diagonalScoring", &Parameters::diagonalScoring)
      .def_readwrite("exactKmerMatching", &Parameters::exactKmerMatching)
      .def_readwrite("maskMode", &Parameters::maskMode)
      .def_readwrite("maskLowerCaseMode", &Parameters::maskLowerCaseMode)
      .def_readwrite("minDiagScoreThr", &Parameters::minDiagScoreThr)
      .def_readwrite("spacedKmer", &Parameters::spacedKmer)
      .def_readwrite("split", &Parameters::split)
      .def_readwrite("splitMode", &Parameters::splitMode)
      .def_readwrite("splitMemoryLimit", &Parameters::splitMemoryLimit)
      .def_readwrite("diskSpaceLimit", &Parameters::diskSpaceLimit)
      .def_readwrite("splitAA", &Parameters::splitAA)
      .def_readwrite("preloadMode", &Parameters::preloadMode)
      .def_readwrite("scoreBias", &Parameters::scoreBias)
      .def_readwrite("realignScoreBias", &Parameters::realignScoreBias)
      .def_readwrite("realignMaxSeqs", &Parameters::realignMaxSeqs)
      .def_readwrite("spacedKmerPattern", &Parameters::spacedKmerPattern)
      .def_readwrite("localTmp", &Parameters::localTmp)
      .def_readwrite("alignmentMode", &Parameters::alignmentMode)
      .def_readwrite("evalThr", &Parameters::evalThr)
      .def_readwrite("covThr", &Parameters::covThr)
      .def_readwrite("covMode", &Parameters::covMode)
      .def_readwrite("seqIdMode", &Parameters::seqIdMode)
      .def_readwrite("maxRejected", &Parameters::maxRejected)
      .def_readwrite("maxAccept", &Parameters::maxAccept)
      .def_readwrite("altAlignment", &Parameters::altAlignment)
      .def_readwrite("seqIdThr", &Parameters::seqIdThr)
      .def_readwrite("alnLenThr", &Parameters::alnLenThr)
      .def_readwrite("addBacktrace", &Parameters::addBacktrace)
      .def_readwrite("realign", &Parameters::realign)
      .def_readwrite("gapOpen", &Parameters::gapOpen)
      .def_readwrite("gapExtend", &Parameters::gapExtend)
      .def_readwrite("zdrop", &Parameters::zdrop)
      .def_readwrite("runner", &Parameters::runner)
      .def_readwrite("reuseLatest", &Parameters::reuseLatest)
      .def_readwrite("clusteringMode", &Parameters::clusteringMode)
      .def_readwrite("clusterSteps", &Parameters::clusterSteps)
      .def_readwrite("singleStepClustering", &Parameters::singleStepClustering)
      .def_readwrite("clusterReassignment", &Parameters::clusterReassignment)
      .def_readwrite("numIterations", &Parameters::numIterations)
      .def_readwrite("startSens", &Parameters::startSens)
      .def_readwrite("sensSteps", &Parameters::sensSteps)
      .def_readwrite("sliceSearch", &Parameters::sliceSearch)
      .def_readwrite("strand", &Parameters::strand)
      .def_readwrite("orfFilter", &Parameters::orfFilter)
      .def_readwrite("orfFilterSens", &Parameters::orfFilterSens)
      .def_readwrite("orfFilterEval", &Parameters::orfFilterEval)
      .def_readwrite("lcaSearch", &Parameters::lcaSearch)
      .def_readwrite("greedyBestHits", &Parameters::greedyBestHits)
      .def_readwrite("maxIteration", &Parameters::maxIteration)
      .def_readwrite("similarityScoreType", &Parameters::similarityScoreType)
      .def_readwrite("orfMinLength", &Parameters::orfMinLength)
      .def_readwrite("orfMaxLength", &Parameters::orfMaxLength)
      .def_readwrite("orfMaxGaps", &Parameters::orfMaxGaps)
      .def_readwrite("contigStartMode", &Parameters::contigStartMode)
      .def_readwrite("contigEndMode", &Parameters::contigEndMode)
      .def_readwrite("orfStartMode", &Parameters::orfStartMode)
      .def_readwrite("forwardFrames", &Parameters::forwardFrames)
      .def_readwrite("reverseFrames", &Parameters::reverseFrames)
      .def_readwrite("useAllTableStarts", &Parameters::useAllTableStarts)
      .def_readwrite("translate", &Parameters::translate)
      .def_readwrite("createLookup", &Parameters::createLookup)
      .def_readwrite("formatAlignmentMode", &Parameters::formatAlignmentMode)
      .def_readwrite("outfmt", &Parameters::outfmt)
      .def_readwrite("dbOut", &Parameters::dbOut)
      .def_readwrite("rescoreMode", &Parameters::rescoreMode)
      .def_readwrite("wrappedScoring", &Parameters::wrappedScoring)
      .def_readwrite("filterHits", &Parameters::filterHits)
      .def_readwrite("globalAlignment", &Parameters::globalAlignment)
      .def_readwrite("sortResults", &Parameters::sortResults)
      .def_readwrite("msaFormatMode", &Parameters::msaFormatMode)
      .def_readwrite("allowDeletion", &Parameters::allowDeletion)
      .def_readwrite("summaryPrefix", &Parameters::summaryPrefix)
      .def_readwrite("skipQuery", &Parameters::skipQuery)
      .def_readwrite("identifierField", &Parameters::identifierField)
      .def_readwrite("matchMode", &Parameters::matchMode)
      .def_readwrite("matchRatio", &Parameters::matchRatio)
      .def_readwrite("maskProfile", &Parameters::maskProfile)
      .def_readwrite("filterMaxSeqId", &Parameters::filterMaxSeqId)
      .def_readwrite("evalProfile", &Parameters::evalProfile)
      .def_readwrite("filterMsa", &Parameters::filterMsa)
      .def_readwrite("qsc", &Parameters::qsc)
      .def_readwrite("qid", &Parameters::qid)
      .def_readwrite("covMSAThr", &Parameters::covMSAThr)
      .def_readwrite("Ndiff", &Parameters::Ndiff)
      .def_readwrite("wg", &Parameters::wg)
      .def_readwrite("pca", &Parameters::pca)
      .def_readwrite("pcb", &Parameters::pcb)
      .def_readwrite("neff", &Parameters::neff)
      .def_readwrite("tau", &Parameters::tau)
      .def_readwrite("firstSeqRepr", &Parameters::firstSeqRepr)
      .def_readwrite("idxSeqSrc", &Parameters::idxSeqSrc)
      .def_readwrite("fullHeader", &Parameters::fullHeader)
      .def_readwrite("targetTsvColumn", &Parameters::targetTsvColumn)
      .def_readwrite("stat", &Parameters::stat)
      .def_readwrite("kmersPerSequence", &Parameters::kmersPerSequence)
      .def_readwrite("kmersPerSequenceScale",
                     &Parameters::kmersPerSequenceScale)
      .def_readwrite("includeOnlyExtendable",
                     &Parameters::includeOnlyExtendable)
      .def_readwrite("ignoreMultiKmer", &Parameters::ignoreMultiKmer)
      .def_readwrite("hashShift", &Parameters::hashShift)
      .def_readwrite("pickNbest", &Parameters::pickNbest)
      .def_readwrite("adjustKmerLength", &Parameters::adjustKmerLength)
      .def_readwrite("resultDirection", &Parameters::resultDirection)
      .def_readwrite("checkCompatible", &Parameters::checkCompatible)
      .def_readwrite("searchType", &Parameters::searchType)
      .def_readwrite("identifierOffset", &Parameters::identifierOffset)
      .def_readwrite("dbType", &Parameters::dbType)
      .def_readwrite("createdbMode", &Parameters::createdbMode)
      .def_readwrite("shuffleDatabase", &Parameters::shuffleDatabase)
      .def_readwrite("sequenceOverlap", &Parameters::sequenceOverlap)
      .def_readwrite("sequenceSplitMode", &Parameters::sequenceSplitMode)
      .def_readwrite("headerSplitMode", &Parameters::headerSplitMode)
      .def_readwrite("useHeaderFile", &Parameters::useHeaderFile)
      .def_readwrite("writeLookup", &Parameters::writeLookup)
      .def_readwrite("useHeader", &Parameters::useHeader)
      .def_readwrite("gffType", &Parameters::gffType)
      .def_readwrite("translationTable", &Parameters::translationTable)
      .def_readwrite("addOrfStop", &Parameters::addOrfStop)
      .def_readwrite("minSequences", &Parameters::minSequences)
      .def_readwrite("maxSequences", &Parameters::maxSequences)
      .def_readwrite("hhFormat", &Parameters::hhFormat)
      .def_readwrite("filterColumn", &Parameters::filterColumn)
      .def_readwrite("columnToTake", &Parameters::columnToTake)
      .def_readwrite("filterColumnRegex", &Parameters::filterColumnRegex)
      .def_readwrite("filteringFile", &Parameters::filteringFile)
      .def_readwrite("mappingFile", &Parameters::mappingFile)
      .def_readwrite("filterExpression", &Parameters::filterExpression)
      .def_readwrite("positiveFilter", &Parameters::positiveFilter)
      .def_readwrite("trimToOneColumn", &Parameters::trimToOneColumn)
      .def_readwrite("extractLines", &Parameters::extractLines)
      .def_readwrite("compValue", &Parameters::compValue)
      .def_readwrite("compOperator", &Parameters::compOperator)
      .def_readwrite("sortEntries", &Parameters::sortEntries)
      .def_readwrite("beatsFirst", &Parameters::beatsFirst)
      .def_readwrite("joinDB", &Parameters::joinDB)
      .def_readwrite("simpleBestHit", &Parameters::simpleBestHit)
      .def_readwrite("alpha", &Parameters::alpha)
      .def_readwrite("shortOutput", &Parameters::shortOutput)
      .def_readwrite("aggregationMode", &Parameters::aggregationMode)
      .def_readwrite("mergePrefixes", &Parameters::mergePrefixes)
      .def_readwrite("mergeStopEmpty", &Parameters::mergeStopEmpty)
      .def_readwrite("overlap", &Parameters::overlap)
      .def_readwrite("msaType", &Parameters::msaType)
      .def_readwrite("extractMode", &Parameters::extractMode)
      .def_readwrite("kbColumns", &Parameters::kbColumns)
      .def_readwrite("preserveKeysB", &Parameters::preserveKeysB)
      .def_readwrite("takeLargerEntry", &Parameters::takeLargerEntry)
      .def_readwrite("chainAlignment", &Parameters::chainAlignment)
      .def_readwrite("mergeQuery", &Parameters::mergeQuery)
      .def_readwrite("outputDbType", &Parameters::outputDbType)
      .def_readwrite("useSequenceId", &Parameters::useSequenceId)
      .def_readwrite("prefix", &Parameters::prefix)
      .def_readwrite("tsvOut", &Parameters::tsvOut)
      .def_readwrite("recoverDeleted", &Parameters::recoverDeleted)
      .def_readwrite("headerType", &Parameters::headerType)
      .def_readwrite("taxonList", &Parameters::taxonList)
      .def_readwrite("idList", &Parameters::idList)
      .def_readwrite("idxEntryType", &Parameters::idxEntryType)
      .def_readwrite("pickIdFrom", &Parameters::pickIdFrom)
      .def_readwrite("lcaRanks", &Parameters::lcaRanks)
      .def_readwrite("showTaxLineage", &Parameters::showTaxLineage)
      .def_readwrite("blacklist", &Parameters::blacklist)
      .def_readwrite("majorityThr", &Parameters::majorityThr)
      .def_readwrite("voteMode", &Parameters::voteMode)
      .def_readwrite("reportMode", &Parameters::reportMode)
      .def_readwrite("ncbiTaxDump", &Parameters::ncbiTaxDump)
      .def_readwrite("taxMappingFile", &Parameters::taxMappingFile)
      .def_readwrite("taxMappingMode", &Parameters::taxMappingMode)
      .def_readwrite("taxDbMode", &Parameters::taxDbMode)
      .def_readwrite("expansionMode", &Parameters::expansionMode)
      .def_readwrite("taxonomySearchMode", &Parameters::taxonomySearchMode)
      .def_readwrite("taxonomyOutputMode", &Parameters::taxonomyOutputMode)
      .def_readwrite("subDbMode", &Parameters::subDbMode)
      .def_readwrite("tarInclude", &Parameters::tarInclude)
      .def_readwrite("tarExclude", &Parameters::tarExclude)
      .def_readwrite("help", &Parameters::help)
      .def_readwrite("citations", &Parameters::citations);

  m.def("_call_biosnake", &call_biosnake_proxy, R"doc(
        Run biosnake2
    )doc");
}
