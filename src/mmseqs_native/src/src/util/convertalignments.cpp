#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "TranslateNucl.h"
#include "Sequence.h"
#include "Orf.h"
#include "MemoryMapped.h"
#include "NcbiTaxonomy.h"
#include "output.h"

#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>
#include "result_viz_prelude.html.zst.h"

#include <map>

#ifdef OPENMP
#include <omp.h>
#endif


void printSeqBasedOnAln(mmseqs_output* mmout, std::string &out, const char *seq, unsigned int offset,
                        const std::string &bt, bool reverse, bool isReverseStrand,
                        bool translateSequence, const TranslateNucl &translateNucl) {
    unsigned int seqPos = 0;
    char codon[3];
    for (uint32_t i = 0; i < bt.size(); ++i) {
        char seqChar = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos]) : seq[offset + seqPos];
        if (translateSequence) {
            codon[0] = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos])     : seq[offset + seqPos];
            codon[1] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+1)]) : seq[offset + (seqPos+1)];
            codon[2] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+2)]) : seq[offset + (seqPos+2)];
            seqChar = translateNucl.translateSingleCodon(codon);
        }
        switch (bt[i]) {
            case 'M':
                out.append(1, seqChar);
                seqPos += (translateSequence) ?  3 : 1;
                break;
            case 'I':
                if (reverse == true) {
                    out.append(1, '-');
                } else {
                    out.append(1, seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                }
                break;
            case 'D':
                if (reverse == true) {
                    out.append(1, seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                } else {
                    out.append(1, '-');
                }
                break;
        }
    }
}

/*
query       Query sequence label
target      Target sequenc label
evalue      E-value
gapopen     Number of gap opens
pident      Percentage of identical matches
nident      Number of identical matches
qstart      1-based start position of alignment in query sequence
qend        1-based end position of alignment in query sequence
qlen        Query sequence length
tstart      1-based start position of alignment in target sequence
tend        1-based end position of alignment in target sequence
tlen        Target sequence length
alnlen      Number of alignment columns
raw         Raw alignment score
bits        Bit score
cigar       Alignment as string M=letter pair, D=delete (gap in query), I=insert (gap in target)
qseq        Full-length query sequence
tseq        Full-length target sequence
qheader     Header of Query sequence
theader     Header of Target sequence
qaln        Aligned query sequence with gaps
taln        Aligned target sequence with gaps
qframe      Query frame (-3 to +3)
tframe      Target frame (-3 to +3)
mismatch    Number of mismatches
qcov        Fraction of query sequence covered by alignment
tcov        Fraction of target sequence covered by alignment
qset        Query set
tset        Target set
 */

std::map<unsigned int, unsigned int> readKeyToSet(const std::string& file) {
    std::map<unsigned int, unsigned int> mapping;
    if (file.length() == 0) {
        return mapping;
    }

    MemoryMapped lookup(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char* data = (char *) lookup.getData();
    const char* entry[255];
    while (*data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        mapping.emplace(Util::fast_atoi<unsigned int>(entry[0]), Util::fast_atoi<unsigned int>(entry[2]));
        data = Util::skipLine(data);
    }
    lookup.close();
    return mapping;
}


std::map<unsigned int, std::string> readSetToSource(const std::string& file) {
    std::map<unsigned int, std::string> mapping;
    if (file.length() == 0) {
        return mapping;
    }

    MemoryMapped source(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char* data = (char *) source.getData();
    const char* entry[255];
    while (*data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 2) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        data = Util::skipLine(data);
        std::string source(entry[1], data - entry[1] - 1);
        mapping.emplace(Util::fast_atoi<unsigned int>(entry[0]), source);
    }
    source.close();
    return mapping;
}

static bool compareToFirstInt(const std::pair<unsigned int, unsigned int>& lhs, const std::pair<unsigned int, unsigned int>&  rhs){
    return (lhs.first <= rhs.first);
}

int convertalignments(mmseqs_output* out, Parameters &par) {
//    Parameters &par = Parameters::getInstance();
//    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const int format = par.formatAlignmentMode;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool needSequenceDB = false;
    bool needBacktrace = false;
    bool needFullHeaders = false;
    bool needLookup = false;
    bool needSource = false;
    bool needTaxonomy = false;
    bool needTaxonomyMapping = false;
    const std::vector<int> outcodes = Parameters::getOutputFormat(format, par.outfmt, needSequenceDB, needBacktrace, needFullHeaders,
                                                                  needLookup, needSource, needTaxonomyMapping, needTaxonomy);

    NcbiTaxonomy * t = NULL;
    std::vector<std::pair<unsigned int, unsigned int>> mapping;
    if(needTaxonomy){
        std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
        t = NcbiTaxonomy::openTaxonomy(db2NoIndexName);
    }
    if(needTaxonomy || needTaxonomyMapping){
        std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
        if(FileUtil::fileExists(std::string(db2NoIndexName + "_mapping").c_str()) == false){
            Debug(Debug::ERROR) << db2NoIndexName + "_mapping" << " does not exist. Please create the taxonomy mapping!\n";
            EXIT(EXIT_FAILURE);
        }
        bool isSorted = Util::readMapping( db2NoIndexName + "_mapping", mapping);
        if(isSorted == false){
            std::stable_sort(mapping.begin(), mapping.end(), compareToFirstInt);
        }
    }

    bool isTranslatedSearch = false;

    int dbaccessMode = needSequenceDB ? (DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA) : (DBReader<unsigned int>::USE_INDEX);

    std::map<unsigned int, unsigned int> qKeyToSet;
    std::map<unsigned int, unsigned int> tKeyToSet;
    if (needLookup) {
        std::string file1 = par.db1 + ".lookup";
        std::string file2 = par.db2 + ".lookup";
        qKeyToSet = readKeyToSet(file1);
        tKeyToSet = readKeyToSet(file2);
    }

    std::map<unsigned int, std::string> qSetToSource;
    std::map<unsigned int, std::string> tSetToSource;
    if (needSource) {
        std::string file1 = par.db1 + ".source";
        std::string file2 = par.db2 + ".source";
        qSetToSource = readSetToSource(file1);
        tSetToSource = readSetToSource(file2);
    }

    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    IndexReader *tDbr;
    IndexReader *tDbrHeader;
    if (sameDB) {
        tDbr = &qDbr;
        tDbrHeader= &qDbrHeader;
    } else {
        tDbr = new IndexReader(par.db2, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tDbrHeader = new IndexReader(par.db2, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    bool queryNucs = Parameters::isEqualDbtype(qDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    bool targetNucs = Parameters::isEqualDbtype(tDbr->sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    if (needSequenceDB) {
        // try to figure out if search was translated. This is can not be solved perfectly.
        bool seqtargetAA = false;
        if(Parameters::isEqualDbtype(tDbr->getDbtype(), Parameters::DBTYPE_INDEX_DB)){
            IndexReader tseqDbr(par.db2, par.threads, IndexReader::SEQUENCES, 0, IndexReader::PRELOAD_INDEX);
            seqtargetAA = Parameters::isEqualDbtype(tseqDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_AMINO_ACIDS);
        } else if(targetNucs == true && queryNucs == true && par.searchType == Parameters::SEARCH_TYPE_AUTO){
            Debug(Debug::WARNING) << "It is unclear from the input if a translated or nucleotide search was performed\n "
                                     "Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)\n";
            EXIT(EXIT_FAILURE);
        } else if(par.searchType == Parameters::SEARCH_TYPE_TRANSLATED){
            seqtargetAA = true;
        }

        if((targetNucs == true && queryNucs == false )  || (targetNucs == false && queryNucs == true ) || (targetNucs == true && seqtargetAA == true && queryNucs == true )  ){
            isTranslatedSearch = true;
        }
    }

    int gapOpen, gapExtend;
    SubstitutionMatrix * subMat= NULL;
    if (targetNucs == true && queryNucs == true && isTranslatedSearch == false) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
        gapOpen = par.gapOpen.nucleotides;
        gapExtend = par.gapExtend.nucleotides;
    }else{
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
        gapOpen = par.gapOpen.aminoacids;
        gapExtend = par.gapExtend.aminoacids;
    }
    EvalueComputation *evaluer = NULL;
    bool queryProfile = false;
    bool targetProfile = false;
    if (needSequenceDB) {
        queryProfile = Parameters::isEqualDbtype(qDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_HMM_PROFILE);
        targetProfile = Parameters::isEqualDbtype(tDbr->sequenceReader->getDbtype(), Parameters::DBTYPE_HMM_PROFILE);
        evaluer = new EvalueComputation(tDbr->sequenceReader->getAminoAcidDBSize(), subMat, gapOpen, gapExtend);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    unsigned int localThreads = 1;
#ifdef OPENMP
    localThreads = std::min((unsigned int)par.threads, (unsigned int)alnDbr.getSize());
#endif

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads, shouldCompress, dbType);
    resultWriter.open();

    const bool isDb = par.dbOut;
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));

    if (format == Parameters::FORMAT_ALIGNMENT_SAM) {
        char buffer[1024];
        unsigned int lastKey = tDbr->sequenceReader->getLastKey();
        bool *headerWritten = new bool[lastKey + 1];
        memset(headerWritten, 0, sizeof(bool) * (lastKey + 1));
        resultWriter.writeStart(0);
        std::string header = "@HD\tVN:1.4\tSO:queryname\n";
        resultWriter.writeAdd(header.c_str(), header.size(), 0);

        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            char *data = alnDbr.getData(i, 0);
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                if (headerWritten[dbKey] == false) {
                    headerWritten[dbKey] = true;
                    unsigned int tId = tDbr->sequenceReader->getId(dbKey);
                    unsigned int seqLen = tDbr->sequenceReader->getSeqLen(tId);
                    unsigned int tHeaderId = tDbrHeader->sequenceReader->getId(dbKey);
                    const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, 0);
                    std::string targetId = Util::parseFastaHeader(tHeader);
                    int count = snprintf(buffer, sizeof(buffer), "@SQ\tSN:%s\tLN:%d\n", targetId.c_str(),
                                         (int32_t) seqLen);
                    if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                        Debug(Debug::WARNING) << "Truncated line in header " << i << "!\n";
                        continue;
                    }
                    resultWriter.writeAdd(buffer, count, 0);
                }
                resultWriter.writeEnd(0, 0, false, 0);
                data = Util::skipLine(data);
            }
        }
        delete[] headerWritten;
    } else if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
        size_t dstSize = ZSTD_findDecompressedSize(result_viz_prelude_html_zst, result_viz_prelude_html_zst_len);
        char* dst = (char*)malloc(sizeof(char) * dstSize);
        size_t realSize = ZSTD_decompress(dst, dstSize, result_viz_prelude_html_zst, result_viz_prelude_html_zst_len);
        resultWriter.writeData(dst, realSize, 0, 0, false, false);
        const char* scriptBlock = "<script>render([";
        resultWriter.writeData(scriptBlock, strlen(scriptBlock), 0, 0, false, false);
        free(dst);
    }

    Debug::Progress progress(alnDbr.getSize());
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024];

        std::string result;
        result.reserve(1024*1024);

        std::string queryProfData;
        queryProfData.reserve(1024);

        std::string queryBuffer;
        queryBuffer.reserve(1024);

        std::string queryHeaderBuffer;
        queryHeaderBuffer.reserve(1024);

        std::string targetProfData;
        targetProfData.reserve(1024);

        std::string newBacktrace;
        newBacktrace.reserve(1024);

        const TaxonNode * taxonNode = NULL;

#pragma omp  for schedule(dynamic, 10)
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            progress.updateProgress();

            const unsigned int queryKey = alnDbr.getDbKey(i);
            char *querySeqData = NULL;
            size_t querySeqLen = 0;
            queryProfData.clear();
            if (needSequenceDB) {
                size_t qId = qDbr.sequenceReader->getId(queryKey);
                querySeqData = qDbr.sequenceReader->getData(qId, thread_idx);
                querySeqLen = qDbr.sequenceReader->getSeqLen(qId);
                if(sameDB && qDbr.sequenceReader->isCompressed()){
                    queryBuffer.assign(querySeqData, querySeqLen);
                    querySeqData = (char*) queryBuffer.c_str();
                }
                if (queryProfile) {
                    Sequence::extractProfileConsensus(querySeqData, *subMat, queryProfData);
                }
            }

            std::vector<mmseqs_blast_tab_record> results_row;

            size_t qHeaderId = qDbrHeader.sequenceReader->getId(queryKey);
            const char *qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
            size_t qHeaderLen = qDbrHeader.sequenceReader->getSeqLen(qHeaderId);
            std::string queryId = Util::parseFastaHeader(qHeader);
            if (sameDB && needFullHeaders) {
                queryHeaderBuffer.assign(qHeader, qHeaderLen);
                qHeader = (char*) queryHeaderBuffer.c_str();
            }

            if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
                const char* jsStart = "{\"query\": {\"accession\": \"%s\",\"sequence\": \"";
                int count = snprintf(buffer, sizeof(buffer), jsStart, queryId.c_str(), querySeqData);
                if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                    Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                    continue;
                }
                result.append(buffer, count);
                if (queryProfile) {
                    result.append(queryProfData);
                } else {
                    result.append(querySeqData, querySeqLen);
                }
                result.append("\"}, \"alignments\": [\n");
            }

            char *data = alnDbr.getData(i, thread_idx);
            while (*data != '\0') {
                mmseqs_blast_tab_record record;
                Matcher::result_t res = Matcher::parseAlignmentRecord(data, true);
                data = Util::skipLine(data);

                if (res.backtrace.empty() && needBacktrace == true) {
                    Debug(Debug::ERROR) << "Backtrace cigar is missing in the alignment result. Please recompute the alignment with the -a flag.\n"
                                           "Command: mmseqs align " << par.db1 << " " << par.db2 << " " << par.db3 << " " << "alnNew -a\n";
                    EXIT(EXIT_FAILURE);
                }

                size_t tHeaderId = tDbrHeader->sequenceReader->getId(res.dbKey);
                const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, thread_idx);
                size_t tHeaderLen = tDbrHeader->sequenceReader->getSeqLen(tHeaderId);
                std::string targetId = Util::parseFastaHeader(tHeader);

                unsigned int gapOpenCount = 0;
                unsigned int alnLen = res.alnLength;
                unsigned int missMatchCount = 0;
                unsigned int identical = 0;
                if (res.backtrace.empty() == false) {
                    size_t matchCount = 0;
                    alnLen = 0;
                    for (size_t pos = 0; pos < res.backtrace.size(); pos++) {
                        int cnt = 0;
                        if (isdigit(res.backtrace[pos])) {
                            cnt += Util::fast_atoi<int>(res.backtrace.c_str() + pos);
                            while (isdigit(res.backtrace[pos])) {
                                pos++;
                            }
                        }
                        alnLen += cnt;

                        switch (res.backtrace[pos]) {
                            case 'M':
                                matchCount += cnt;
                                break;
                            case 'D':
                            case 'I':
                                gapOpenCount += 1;
                                break;
                        }
                    }
//                res.seqId = X / alnLen;
                    identical = static_cast<unsigned int>(res.seqId * static_cast<float>(alnLen) + 0.5);
                    //res.alnLength = alnLen;
                    missMatchCount = static_cast<unsigned int>( matchCount - identical);
                } else {
                    const int adjustQstart = (res.qStartPos == -1) ? 0 : res.qStartPos;
                    const int adjustDBstart = (res.dbStartPos == -1) ? 0 : res.dbStartPos;
                    const float bestMatchEstimate = static_cast<float>(std::min(abs(res.qEndPos - adjustQstart), abs(res.dbEndPos - adjustDBstart)));
                    missMatchCount = static_cast<unsigned int>(bestMatchEstimate * (1.0f - res.seqId) + 0.5);
                }

                switch (format) {
                    case Parameters::FORMAT_ALIGNMENT_BLAST_TAB: {
                        if (outcodes.empty()) {
                            int count = snprintf(buffer, sizeof(buffer),
                                                 "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                                                 queryId.c_str(), targetId.c_str(), res.seqId, alnLen,
                                                 missMatchCount, gapOpenCount,
                                                 res.qStartPos + 1, res.qEndPos + 1,
                                                 res.dbStartPos + 1, res.dbEndPos + 1,
                                                 res.eval, res.score);
                            if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                                Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                                continue;
                            }
                            //results_row.push_back(std::string(buffer, count));
                        } else {
                            char *targetSeqData = NULL;
                            targetProfData.clear();
                            unsigned int taxon = 0;

                            if(needTaxonomy || needTaxonomyMapping) {
                                std::pair<unsigned int, unsigned int> val;
                                val.first = res.dbKey;
                                std::vector<std::pair<unsigned int, unsigned int>>::iterator mappingIt;
                                mappingIt = std::upper_bound(mapping.begin(), mapping.end(), val, compareToFirstInt);
                                if (mappingIt == mapping.end() || mappingIt->first != val.first) {
                                    taxon = 0;
                                    taxonNode = NULL;
                                }else{
                                    taxon = mappingIt->second;
                                    if(needTaxonomy){
                                        taxonNode = t->taxonNode(taxon, false);
                                    }
                                }

                            }

                            if (needSequenceDB) {
                                size_t tId = tDbr->sequenceReader->getId(res.dbKey);
                                targetSeqData = tDbr->sequenceReader->getData(tId, thread_idx);
                                if (targetProfile) {
                                    Sequence::extractProfileConsensus(targetSeqData, *subMat, targetProfData);
                                }
                            }
                            for(size_t i = 0; i < outcodes.size(); i++) {
                                switch (outcodes[i]) {
                                    case Parameters::OUTFMT_QUERY:
                                        //results_row.push_back(queryId);
                                        record.query_sequence_id = queryId;
                                        break;
                                    case Parameters::OUTFMT_TARGET:
                                        record.target_sequence_id = targetId;
                                        //results_row.push_back(targetId);
                                        break;
                                    case Parameters::OUTFMT_EVALUE:
                                        record.e_value = res.eval;
                                        //results_row.push_back(SSTR(res.eval));
                                        break;
                                    case Parameters::OUTFMT_GAPOPEN:
                                        record.number_of_gap_openings = gapOpenCount;
                                        //results_row.push_back(SSTR(gapOpenCount));
                                        break;
                                    case Parameters::OUTFMT_FIDENT:
                                        record.sequence_identity = res.seqId;
                                        //results_row.push_back(SSTR(res.seqId));
                                        break;
                                    case Parameters::OUTFMT_PIDENT:
                                        //results_row.push_back(SSTR(res.seqId*100));
                                        break;
                                    case Parameters::OUTFMT_NIDENT:
                                        //results_row.push_back(SSTR(identical));
                                        break;
                                    case Parameters::OUTFMT_QSTART:
                                        record.domain_start_index_query = res.qStartPos + 1;
                                        //results_row.push_back(SSTR(res.qStartPos + 1));
                                        break;
                                    case Parameters::OUTFMT_QEND:
                                        record.domain_end_index_query = res.qEndPos + 1;
                                        //results_row.push_back(SSTR(res.qEndPos + 1));
                                        break;
                                    case Parameters::OUTFMT_QLEN:
                                        //results_row.push_back(SSTR(res.qLen));
                                        break;
                                    case Parameters::OUTFMT_TSTART:
                                        record.domain_start_index_target = res.dbStartPos + 1;
                                        //results_row.push_back(SSTR(res.dbStartPos + 1));
                                        break;
                                    case Parameters::OUTFMT_TEND:
                                        record.domain_end_index_target = res.dbEndPos + 1;
                                        //results_row.push_back(SSTR(res.dbEndPos + 1));
                                        break;
                                    case Parameters::OUTFMT_TLEN:
                                        //results_row.push_back(SSTR(res.dbLen));
                                        break;
                                    case Parameters::OUTFMT_ALNLEN:
                                        record.alignment_length = alnLen;
                                        //results_row.push_back(SSTR(alnLen));
                                        break;
                                    case Parameters::OUTFMT_RAW:
                                        //results_row.push_back(SSTR(static_cast<int>(evaluer->computeRawScoreFromBitScore(res.score) + 0.5)));
                                        break;
                                    case Parameters::OUTFMT_BITS:
                                        record.bit_score = res.score;
                                        //results_row.push_back(SSTR(res.score));
                                        break;
                                    case Parameters::OUTFMT_CIGAR:
                                        if(isTranslatedSearch == true && targetNucs == true && queryNucs == true ){
                                            Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                                            res.backtrace = newBacktrace;
                                        }
                                        //results_row.push_back(SSTR(res.backtrace));
                                        newBacktrace.clear();
                                        break;
                                    case Parameters::OUTFMT_QSEQ:
                                        if (queryProfile) {
                                            //results_row.push_back(std::string(queryProfData.c_str(), res.qLen));
                                        } else {
                                            //results_row.push_back(std::string(querySeqData, res.qLen));
                                        }
                                        break;
                                    case Parameters::OUTFMT_TSEQ:
                                        if (targetProfile) {
                                            //results_row.push_back(std::string(targetProfData.c_str(), res.dbLen));
                                        } else {
                                            //results_row.push_back(std::string(targetSeqData, res.dbLen));
                                        }
                                        break;
                                    case Parameters::OUTFMT_QHEADER:
                                        //results_row.push_back(std::string(qHeader, qHeaderLen));
                                        break;
                                    case Parameters::OUTFMT_THEADER:
                                        //results_row.push_back(std::string(tHeader, tHeaderLen));
                                        break;
                                    case Parameters::OUTFMT_QALN:
                                        if (queryProfile) {
                                            printSeqBasedOnAln(out, result, queryProfData.c_str(), res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                        } else {
                                            printSeqBasedOnAln(out, result, querySeqData, res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                        }
                                        break;
                                    case Parameters::OUTFMT_TALN: {
                                        if (targetProfile) {
                                            printSeqBasedOnAln(out, result, targetProfData.c_str(), res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                        } else {
                                            printSeqBasedOnAln(out, result, targetSeqData, res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                        }
                                        break;
                                    }
                                    case Parameters::OUTFMT_MISMATCH:
                                        record.number_of_mismatches = missMatchCount;
                                        //results_row.push_back(SSTR(missMatchCount));
                                        break;
                                    case Parameters::OUTFMT_QCOV:
                                        //results_row.push_back(SSTR(res.qcov));
                                        break;
                                    case Parameters::OUTFMT_TCOV:
                                        //results_row.push_back(SSTR(res.dbcov));
                                        break;
                                    case Parameters::OUTFMT_QSET:
                                        //results_row.push_back(SSTR(qSetToSource[qKeyToSet[queryKey]]));
                                        break;
                                    case Parameters::OUTFMT_QSETID:
                                        //results_row.push_back(SSTR(qKeyToSet[queryKey]));
                                        break;
                                    case Parameters::OUTFMT_TSET:
                                        //results_row.push_back(SSTR(tSetToSource[tKeyToSet[res.dbKey]]));
                                        break;
                                    case Parameters::OUTFMT_TSETID:
                                        //results_row.push_back(SSTR(tKeyToSet[res.dbKey]));
                                        break;
                                    case Parameters::OUTFMT_TAXID:
                                        //results_row.push_back(SSTR(taxon));
                                        break;
                                    case Parameters::OUTFMT_TAXNAME:
                                        //results_row.push_back((taxonNode != NULL) ? t->getString(taxonNode->nameIdx) : "unclassified");
                                        break;
                                    case Parameters::OUTFMT_TAXLIN:
                                        //results_row.push_back((taxonNode != NULL) ? t->taxLineage(taxonNode, true) : "unclassified");
                                        break;
                                    case Parameters::OUTFMT_EMPTY:
                                        //results_row.push_back("-");
                                        break;
                                    case Parameters::OUTFMT_QORFSTART:
                                        //results_row.push_back(SSTR(res.queryOrfStartPos));
                                        break;
                                    case Parameters::OUTFMT_QORFEND:
                                        //results_row.push_back(SSTR(res.queryOrfEndPos));
                                        break;
                                    case Parameters::OUTFMT_TORFSTART:
                                        //results_row.push_back(SSTR(res.dbOrfStartPos));
                                        break;
                                    case Parameters::OUTFMT_TORFEND:
                                        //results_row.push_back(SSTR(res.dbOrfEndPos));
                                        break;
                                }
                                if (i < outcodes.size() - 1) {
                                    result.push_back('\t');
                                }
                            }
                            result.push_back('\n');
                        }
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_BLAST_WITH_LEN: {
                        int count = snprintf(buffer, sizeof(buffer),
                                             "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\t%d\t%d\n",
                                             queryId.c_str(), targetId.c_str(), res.seqId, alnLen,
                                             missMatchCount, gapOpenCount,
                                             res.qStartPos + 1, res.qEndPos + 1,
                                             res.dbStartPos + 1, res.dbEndPos + 1,
                                             res.eval, res.score,
                                             res.qLen, res.dbLen);

                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }

                        result.append(buffer, count);
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_SAM: {
                        bool strand = res.qEndPos > res.qStartPos;
                        int rawScore = static_cast<int>(evaluer->computeRawScoreFromBitScore(res.score) + 0.5);
                        uint32_t mapq = -4.343 * log(exp(static_cast<double>(-rawScore)));
                        mapq = (uint32_t) (mapq + 4.99);
                        mapq = mapq < 254 ? mapq : 254;
                        int count = snprintf(buffer, sizeof(buffer), "%s\t%d\t%s\t%d\t%d\t",  queryId.c_str(), (strand) ? 16: 0, targetId.c_str(), res.dbStartPos + 1, mapq);
                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        if (isTranslatedSearch == true && targetNucs == true && queryNucs == true) {
                            Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                            result.append(newBacktrace);
                            newBacktrace.clear();

                        } else {
                            result.append(res.backtrace);
                        }
                        result.append("\t*\t0\t0\t");
                        int start = std::min(res.qStartPos, res.qEndPos);
                        int end   = std::max(res.qStartPos, res.qEndPos);
                        if (queryProfile) {
                            result.append(queryProfData.c_str() + start, (end + 1) - start);
                        } else {
                            result.append(querySeqData + start, (end + 1) - start);
                        }
                        count = snprintf(buffer, sizeof(buffer), "\t*\tAS:i:%d\tNM:i:%d\n", rawScore, missMatchCount);
                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_HTML: {
                        const char* jsAln = "{\"target\": \"%s\", \"seqId\": %1.3f, \"alnLen\": %d, \"mismatch\": %d, \"gapopen\": %d, \"qStartPos\": %d, \"qEndPos\": %d, \"dbStartPos\": %d, \"dbEndPos\": %d, \"eval\": %.2E, \"score\": %d, \"qLen\": %d, \"dbLen\": %d, \"qAln\": \"";
                        int count = snprintf(buffer, sizeof(buffer), jsAln,
                                             targetId.c_str(), res.seqId, alnLen,
                                             missMatchCount, gapOpenCount,
                                             res.qStartPos + 1, res.qEndPos + 1,
                                             res.dbStartPos + 1, res.dbEndPos + 1,
                                             res.eval, res.score,
                                             res.qLen, res.dbLen);

                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        if (queryProfile) {
                            printSeqBasedOnAln(out, result, queryProfData.c_str(), res.qStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                        } else {
                            printSeqBasedOnAln(out, result, querySeqData, res.qStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                        }
                        result.append("\", \"dbAln\": \"");
                        size_t tId = tDbr->sequenceReader->getId(res.dbKey);
                        char* targetSeqData = tDbr->sequenceReader->getData(tId, thread_idx);
                        if (targetProfile) {
                            Sequence::extractProfileConsensus(targetSeqData, *subMat, targetProfData);
                            printSeqBasedOnAln(out, result, targetProfData.c_str(), res.dbStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), true,
                                               (res.dbStartPos > res.dbEndPos),
                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                        } else {
                            printSeqBasedOnAln(out, result, targetSeqData, res.dbStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), true,
                                               (res.dbStartPos > res.dbEndPos),
                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                        }
                        result.append("\" },\n");
                        break;
                    }

//                    case Parameters::FORMAT_ALIGNMENT_GFF:{
//                        // for TBLASTX
//                        bool strand = res.qEndPos > res.qStartPos;
//                        int currStart = std::min(res.qStartPos, res.qEndPos);
//                        int currEnd = std::max(res.qStartPos, res.qEndPos);
//                        int currLen = currEnd - currStart;
//                        result.append(queryId);
//                        result.append("\tconserve\tprotein_match\t");
//                        result.append(SSTR(currStart+1));
//                        result.push_back('\t');
//                        result.append(SSTR(currEnd+1));
//                        result.push_back('\t');
//                        result.append(SSTR(currLen));
//                        result.push_back('\t');
//                        result.push_back((strand) ? '-' : '+');
//                        result.append("\t.\t");
//                        result.append("ID=");
//                        result.append(queryId);
//                        result.append(":hsp:");
//                        result.append(SSTR(counter));
//                        result.append(";");
//                        break;
//                    }
                    default:
                        Debug(Debug::ERROR) << "Not implemented yet";
                        EXIT(EXIT_FAILURE);
                }
                results_row.push_back(record);
            }

            if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
                result.append("]},\n");
            }
            out->blast_tab_records.push_back(results_row);
            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx, isDb);
            result.clear();
        }
    }
    if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
        const char* endBlock = "]);</script>";
        resultWriter.writeData(endBlock, strlen(endBlock), 0, localThreads - 1, false, false);
    }
    // tsv output
    resultWriter.close(true);
    if (isDb == false) {
        FileUtil::remove(par.db4Index.c_str());
    }
    if(needTaxonomy){
        delete t;
    }
    alnDbr.close();
    if (sameDB == false) {
        delete tDbr;
        delete tDbrHeader;
    }
    if (needSequenceDB) {
        delete evaluer;
    }
    delete subMat;

    return EXIT_SUCCESS;
}

