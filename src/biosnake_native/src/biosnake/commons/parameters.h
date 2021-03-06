// Written by Martin Steinegger martin.steinegger@snu.ac.kr
//
// Represents parameters of Biosnake2
//
#ifndef BIOSNAKE_PARAMETERS
#define BIOSNAKE_PARAMETERS
#include <cstddef>
#include <map>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include <biosnake/commons/command.h>
#include <biosnake/commons/multiParam.h>

static const char *version = "1.0";

#define PARAMETER(x)                     \
  const static int x##_ID = __COUNTER__; \
  BiosnakeParameter x;

struct BiosnakeParameter {
  const char *name;
  const char *display;
  const char *description;
  const std::type_info &type;
  void *value;
  const char *regex;
  const int uniqid;
  unsigned int category;
  bool wasSet;

  static const unsigned int COMMAND_PREFILTER = 1;
  static const unsigned int COMMAND_ALIGN = 2;
  static const unsigned int COMMAND_CLUST = 4;
  static const unsigned int COMMAND_COMMON = 8;
  static const unsigned int COMMAND_PROFILE = 16;
  static const unsigned int COMMAND_MISC = 32;
  static const unsigned int COMMAND_CLUSTLINEAR = 64;
  static const unsigned int COMMAND_EXPERT = 128;
  static const unsigned int COMMAND_HIDDEN = 256;

  BiosnakeParameter(int uid, const char *n, const char *display, const char *d,
                  const std::type_info &hash, void *value, const char *regex,
                  unsigned int category = COMMAND_MISC)
      : name(n),
        display(display),
        description(d),
        type(hash),
        value(value),
        regex(regex),
        uniqid(uid),
        category(category),
        wasSet(false) {}

  void addCategory(unsigned int cat) { category |= cat; }

  void removeCategory(unsigned int cat) { category &= ~cat; }

  void replaceCategory(unsigned int cat) { category = cat; }
};

class Parameters {
 public:
  static constexpr int DBTYPE_AMINO_ACIDS = 0;
  static constexpr int DBTYPE_NUCLEOTIDES = 1;
  static constexpr int DBTYPE_HMM_PROFILE = 2;
  static constexpr int DBTYPE_PROFILE_STATE_SEQ = 3;
  static constexpr int DBTYPE_PROFILE_STATE_PROFILE = 4;
  static constexpr int DBTYPE_ALIGNMENT_RES = 5;
  static constexpr int DBTYPE_CLUSTER_RES = 6;
  static constexpr int DBTYPE_PREFILTER_RES = 7;
  static constexpr int DBTYPE_TAXONOMICAL_RESULT = 8;
  static constexpr int DBTYPE_INDEX_DB = 9;
  static constexpr int DBTYPE_CA3M_DB = 10;
  static constexpr int DBTYPE_MSA_DB = 11;
  static constexpr int DBTYPE_GENERIC_DB = 12;
  static constexpr int DBTYPE_OMIT_FILE = 13;
  static constexpr int DBTYPE_PREFILTER_REV_RES = 14;
  static constexpr int DBTYPE_OFFSETDB = 15;
  static constexpr int DBTYPE_DIRECTORY = 16;  // needed for verification
  static constexpr int DBTYPE_FLATFILE = 17;   // needed for verification
  static constexpr int DBTYPE_SEQTAXDB = 18;   // needed for verification
  static constexpr int DBTYPE_STDIN = 19;      // needed for verification

  // don't forget to add new database types to DBReader::getDbTypeName and
  // Parameters::PARAM_OUTPUT_DBTYPE

  static constexpr int SEARCH_TYPE_AUTO = 0;
  static constexpr int SEARCH_TYPE_PROTEIN = 1;
  static constexpr int SEARCH_TYPE_TRANSLATED = 2;
  static constexpr int SEARCH_TYPE_NUCLEOTIDES = 3;
  static constexpr int SEARCH_TYPE_TRANS_NUCL_ALN = 4;
  // flag
  static constexpr int SEARCH_MODE_FLAG_QUERY_AMINOACID = 1;
  static constexpr int SEARCH_MODE_FLAG_TARGET_AMINOACID = 2;
  static constexpr int SEARCH_MODE_FLAG_QUERY_TRANSLATED = 4;
  static constexpr int SEARCH_MODE_FLAG_TARGET_TRANSLATED = 8;
  static constexpr int SEARCH_MODE_FLAG_QUERY_PROFILE = 16;
  static constexpr int SEARCH_MODE_FLAG_TARGET_PROFILE = 32;
  static constexpr int SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE = 64;
  static constexpr int SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE = 128;

  static constexpr unsigned int ALIGNMENT_MODE_FAST_AUTO = 0;
  static constexpr unsigned int ALIGNMENT_MODE_SCORE_ONLY = 1;
  static constexpr unsigned int ALIGNMENT_MODE_SCORE_COV = 2;
  static constexpr unsigned int ALIGNMENT_MODE_SCORE_COV_SEQID = 3;
  static constexpr unsigned int ALIGNMENT_MODE_UNGAPPED = 4;

  static constexpr unsigned int EXPAND_TRANSFER_EVALUE = 0;
  static constexpr unsigned int EXPAND_RESCORE_BACKTRACE = 1;

  static constexpr unsigned int WRITER_ASCII_MODE = 0;
  static constexpr unsigned int WRITER_COMPRESSED_MODE = 1;
  static constexpr unsigned int WRITER_LEXICOGRAPHIC_MODE = 2;

  // convertalis alignment
  static constexpr int FORMAT_ALIGNMENT_BLAST_TAB = 0;
  static constexpr int FORMAT_ALIGNMENT_SAM = 1;
  static constexpr int FORMAT_ALIGNMENT_BLAST_WITH_LEN = 2;
  static constexpr int FORMAT_ALIGNMENT_HTML = 3;

  // result2msa
  static constexpr int FORMAT_MSA_CA3M = 0;
  static constexpr int FORMAT_MSA_CA3M_CONSENSUS = 1;
  static constexpr int FORMAT_MSA_FASTADB = 2;
  static constexpr int FORMAT_MSA_FASTADB_SUMMARY = 3;
  static constexpr int FORMAT_MSA_STOCKHOLM_FLAT = 4;

  // outfmt
  static constexpr int OUTFMT_QUERY = 0;
  static constexpr int OUTFMT_TARGET = 1;
  static constexpr int OUTFMT_EVALUE = 2;
  static constexpr int OUTFMT_GAPOPEN = 3;
  static constexpr int OUTFMT_PIDENT = 4;
  static constexpr int OUTFMT_NIDENT = 5;
  static constexpr int OUTFMT_QSTART = 6;
  static constexpr int OUTFMT_QEND = 7;
  static constexpr int OUTFMT_QLEN = 8;
  static constexpr int OUTFMT_TSTART = 9;
  static constexpr int OUTFMT_TEND = 10;
  static constexpr int OUTFMT_TLEN = 11;
  static constexpr int OUTFMT_ALNLEN = 12;
  static constexpr int OUTFMT_RAW = 13;
  static constexpr int OUTFMT_BITS = 14;
  static constexpr int OUTFMT_CIGAR = 15;
  static constexpr int OUTFMT_QSEQ = 16;
  static constexpr int OUTFMT_TSEQ = 17;
  static constexpr int OUTFMT_QHEADER = 18;
  static constexpr int OUTFMT_THEADER = 19;
  static constexpr int OUTFMT_QALN = 20;
  static constexpr int OUTFMT_TALN = 21;
  static constexpr int OUTFMT_QFRAME = 22;
  static constexpr int OUTFMT_TFRAME = 23;
  static constexpr int OUTFMT_MISMATCH = 24;
  static constexpr int OUTFMT_QCOV = 25;
  static constexpr int OUTFMT_TCOV = 26;
  static constexpr int OUTFMT_EMPTY = 27;
  static constexpr int OUTFMT_QSET = 28;
  static constexpr int OUTFMT_QSETID = 29;
  static constexpr int OUTFMT_TSET = 30;
  static constexpr int OUTFMT_TSETID = 31;
  static constexpr int OUTFMT_TAXID = 32;
  static constexpr int OUTFMT_TAXNAME = 33;
  static constexpr int OUTFMT_TAXLIN = 34;
  static constexpr int OUTFMT_QORFSTART = 35;
  static constexpr int OUTFMT_QORFEND = 36;
  static constexpr int OUTFMT_TORFSTART = 37;
  static constexpr int OUTFMT_TORFEND = 38;
  static constexpr int OUTFMT_FIDENT = 39;

  static std::vector<int> getOutputFormat(
      biosnake_output* out,
      int formatMode, const std::string &outformat, bool &needSequences,
      bool &needBacktrace, bool &needFullHeaders, bool &needLookup,
      bool &needSource, bool &needTaxonomyMapping, bool &needTaxonomy);

  // clustering
  static constexpr int SET_COVER = 0;
  static constexpr int CONNECTED_COMPONENT = 1;
  static constexpr int GREEDY = 2;
  static constexpr int GREEDY_MEM = 3;

  // clustering
  static constexpr int APC_ALIGNMENTSCORE = 1;
  static constexpr int APC_SEQID = 2;
  // split mode
  static constexpr int TARGET_DB_SPLIT = 0;
  static constexpr int QUERY_DB_SPLIT = 1;
  static constexpr int DETECT_BEST_DB_SPLIT = 2;

  // taxonomy output
  static constexpr int TAXONOMY_OUTPUT_LCA = 0;
  static constexpr int TAXONOMY_OUTPUT_ALIGNMENT = 1;
  static constexpr int TAXONOMY_OUTPUT_BOTH = 2;

  // aggregate taxonomy
  static constexpr int AGG_TAX_UNIFORM = 0;
  static constexpr int AGG_TAX_MINUS_LOG_EVAL = 1;
  static constexpr int AGG_TAX_SCORE = 2;

  // taxonomy search strategy
  static constexpr int TAXONOMY_SINGLE_SEARCH = 1;
  static constexpr int TAXONOMY_2BLCA = 2;
  static constexpr int TAXONOMY_ACCEL_2BLCA = 3;
  static constexpr int TAXONOMY_TOP_HIT = 4;

  static constexpr int PARSE_VARIADIC = 1;
  static constexpr int PARSE_REST = 2;
  static constexpr int PARSE_ALLOW_EMPTY = 4;

  // preload mode
  static constexpr int PRELOAD_MODE_AUTO = 0;
  static constexpr int PRELOAD_MODE_FREAD = 1;
  static constexpr int PRELOAD_MODE_MMAP = 2;
  static constexpr int PRELOAD_MODE_MMAP_TOUCH = 3;

  static std::string getSplitModeName(int splitMode) {
    switch (splitMode) {
      case 0:
        return "Target";
      case 1:
        return "Query";
      case 2:
        return "Auto";
      default:
        return "Error";
    }
  };

  // split
  static constexpr int AUTO_SPLIT_DETECTION = 0;

  static constexpr int MAX_SEQ_LEN = 65535;

  // extractalignedregion
  static constexpr int EXTRACT_QUERY = 1;
  static constexpr int EXTRACT_TARGET = 2;

  static constexpr int CLUST_HASH_DEFAULT_ALPH_SIZE = 3;
  static constexpr int CLUST_HASH_DEFAULT_MIN_SEQ_ID = 99;
  static constexpr int CLUST_LINEAR_DEFAULT_ALPH_SIZE = 13;
  static constexpr int CLUST_LINEAR_DEFAULT_K = 0;
  static constexpr int CLUST_LINEAR_KMER_PER_SEQ = 0;

  // cov mode
  static constexpr int COV_MODE_BIDIRECTIONAL = 0;
  static constexpr int COV_MODE_TARGET = 1;
  static constexpr int COV_MODE_QUERY = 2;
  static constexpr int COV_MODE_LENGTH_QUERY = 3;
  static constexpr int COV_MODE_LENGTH_TARGET = 4;
  static constexpr int COV_MODE_LENGTH_SHORTER = 5;

  // seq. id mode
  static constexpr int SEQ_ID_ALN_LEN = 0;
  static constexpr int SEQ_ID_SHORT = 1;
  static constexpr int SEQ_ID_LONG = 2;

  // seq. split mode
  static constexpr int SEQUENCE_SPLIT_MODE_HARD = 0;
  static constexpr int SEQUENCE_SPLIT_MODE_SOFT = 1;

  // rescorediagonal
  static constexpr int RESCORE_MODE_HAMMING = 0;
  static constexpr int RESCORE_MODE_SUBSTITUTION = 1;
  static constexpr int RESCORE_MODE_ALIGNMENT = 2;
  static constexpr int RESCORE_MODE_GLOBAL_ALIGNMENT = 3;
  static constexpr int RESCORE_MODE_WINDOW_QUALITY_ALIGNMENT = 4;

  // combinepvalperset
  static constexpr int AGGREGATION_MODE_MULTIHIT = 0;
  static constexpr int AGGREGATION_MODE_MIN_PVAL = 1;
  static constexpr int AGGREGATION_MODE_PRODUCT = 2;
  static constexpr int AGGREGATION_MODE_TRUNCATED_PRODUCT = 3;

  // header type
  static constexpr int HEADER_TYPE_UNICLUST = 1;
  static constexpr int HEADER_TYPE_METACLUST = 2;

  // createsubdb, filtertaxseqdb type
  static constexpr int SUBDB_MODE_HARD = 0;
  static constexpr int SUBDB_MODE_SOFT = 1;

  // result direction
  static constexpr int PARAM_RESULT_DIRECTION_QUERY = 0;
  static constexpr int PARAM_RESULT_DIRECTION_TARGET = 1;

  std::string baseTmpPath = "";

  std::string logFilePath = "";

  // path to databases
  std::string db1;
  std::string db1Index;
  std::string db1dbtype;

  std::string hdr1;
  std::string hdr1Index;
  std::string hdr1dbtype;

  std::string db2;
  std::string db2Index;
  std::string db2dbtype;

  std::string hdr2;
  std::string hdr2Index;
  std::string hdr2dbtype;

  std::string db3;
  std::string db3Index;
  std::string db3dbtype;

  std::string hdr3;
  std::string hdr3Index;
  std::string hdr3dbtype;

  std::string db4;
  std::string db4Index;
  std::string db4dbtype;

  std::string hdr4;
  std::string hdr4Index;
  std::string hdr4dbtype;

  std::string db5;
  std::string db5Index;
  std::string db5dbtype;

  std::string hdr5;
  std::string hdr5Index;
  std::string hdr5dbtype;

  std::string db6;
  std::string db6Index;
  std::string db6dbtype;

  std::string hdr6;
  std::string hdr6Index;
  std::string hdr6dbtype;

  std::vector<std::string> filenames;

  const char **restArgv;
  int restArgc;

  MultiParam<char *> scoringMatrixFile;      // path to scoring matrix
  MultiParam<char *> seedScoringMatrixFile;  // seed sub. matrix
  size_t maxSeqLen;                          // sequence length
  size_t maxResListLen;  // Maximal result list length per query
  int verbosity;         // log level
  int threads;           // Amounts of threads
  int compressed;        // compressed writer
  bool removeTmpFiles;   // Do not delete temp files
  bool includeIdentity;  // include identical ids as hit

  // PREFILTER
  float sensitivity;             // target sens
  int kmerSize;                  // kmer size for the prefilter
  int kmerScore;                 // kmer score for the prefilter
  MultiParam<int> alphabetSize;  // alphabet size for the prefilter
  int compBiasCorrection;        // Aminoacid composiont correction
  bool diagonalScoring;          // switch diagonal scoring
  int exactKmerMatching;         // only exact k-mer matching
  int maskMode;                  // mask low complex areas
  int maskLowerCaseMode;         // mask lowercase letters in prefilter and
                                 // kmermatchers

  int minDiagScoreThr;      // min diagonal score
  int spacedKmer;           // Spaced Kmers
  int split;                // Split database in n equal chunks
  int splitMode;            // Split by query or target DB
  size_t splitMemoryLimit;  // Maximum memory in bytes a split can use
  size_t diskSpaceLimit;    // Maximum disk space in bytes for sliced reverse
                            // profile search
  bool splitAA;             // Split database by amino acid count instead
  int preloadMode;          // Preload mode of database
  float scoreBias;  // Add this bias to the score when computing the alignements
  float realignScoreBias;         // Add this bias additionally when realigning
  int realignMaxSeqs;             // Max alignments to realign
  std::string spacedKmerPattern;  // User-specified kmer pattern
  std::string localTmp;           // Local temporary path

  // ALIGNMENT
  int alignmentMode;  // alignment mode 0=fastest on parameters,
                      // 1=score only, 2=score, cov, start/end pos, 3=score,
                      // cov, start/end pos, seq.id,
  double evalThr;     // e-value threshold for acceptance
  float covThr;       // coverage query&target threshold for acceptance
  int covMode;        // coverage target threshold for acceptance
  int seqIdMode;      // seq. id. normalize mode

  int maxRejected;            // after n sequences that are above eval stop
  int maxAccept;              // after n accepted sequences stop
  int altAlignment;           // show up to this many alternative alignments
  float seqIdThr;             // sequence identity threshold for acceptance
  int alnLenThr;              // min. alignment length
  bool addBacktrace;          // store backtrace string (M=Match, D=deletion,
                              // I=insertion)
  bool realign;               // realign hit with more conservative score
  MultiParam<int> gapOpen;    // gap open cost
  MultiParam<int> gapExtend;  // gap extension cost
  int zdrop;                  // zdrop

  // workflow
  std::string runner;
  bool reuseLatest;

  // CLUSTERING
  int clusteringMode;
  int clusterSteps;
  bool singleStepClustering;
  int clusterReassignment;

  // SEARCH WORKFLOW
  int numIterations;
  float startSens;
  int sensSteps;
  bool sliceSearch;
  int strand;
  int orfFilter;
  float orfFilterSens;
  double orfFilterEval;
  bool lcaSearch;

  // easysearch
  bool greedyBestHits;

  // CLUSTERING
  int maxIteration;  // Maximum depth of breadth first search in connected
                     // component
  int similarityScoreType;  // Type of score to use for reassignment 1=alignment
                            // score. 2=coverage 3=sequence identity 4=E-value
                            // 5= Score per Column

  // extractorfs
  int orfMinLength;
  int orfMaxLength;
  int orfMaxGaps;
  int contigStartMode;
  int contigEndMode;
  int orfStartMode;
  std::string forwardFrames;
  std::string reverseFrames;
  bool useAllTableStarts;
  int translate;
  int createLookup;

  // convertalis
  int formatAlignmentMode;
  std::string outfmt;
  bool dbOut;

  // rescorediagonal
  int rescoreMode;
  bool wrappedScoring;
  bool filterHits;
  bool globalAlignment;
  int sortResults;

  // result2msa
  int msaFormatMode;
  bool allowDeletion;
  std::string summaryPrefix;
  bool skipQuery;

  // convertmsa
  int identifierField;

  // msa2profile
  int matchMode;
  float matchRatio;

  // result2profile
  int maskProfile;
  float filterMaxSeqId;
  double evalProfile;
  int filterMsa;
  float qsc;
  float qid;
  float covMSAThr;
  int Ndiff;
  bool wg;
  float pca;
  float pcb;

  // sequence2profile
  float neff;
  float tau;

  // createtsv
  bool firstSeqRepr;
  int idxSeqSrc;
  bool fullHeader;
  size_t targetTsvColumn;

  // result2stats
  std::string stat;

  // linearcluster
  int kmersPerSequence;
  MultiParam<float> kmersPerSequenceScale;
  bool includeOnlyExtendable;
  bool ignoreMultiKmer;
  int hashShift;
  int pickNbest;
  int adjustKmerLength;
  int resultDirection;

  // indexdb
  int checkCompatible;
  int searchType;

  // createdb
  int identifierOffset;
  int dbType;
  int createdbMode;
  bool shuffleDatabase;

  // splitsequence
  int sequenceOverlap;
  int sequenceSplitMode;
  int headerSplitMode;

  // convert2fasta
  bool useHeaderFile;
  int writeLookup;

  // result2flat
  bool useHeader;

  // gff2db
  std::string gffType;

  // translate nucleotide
  int translationTable;
  bool addOrfStop;

  // createseqfiledb
  int minSequences;
  int maxSequences;
  bool hhFormat;

  // filterDb
  int filterColumn;
  int columnToTake;
  std::string filterColumnRegex;
  std::string filteringFile;
  std::string mappingFile;
  std::string filterExpression;
  bool positiveFilter;
  bool trimToOneColumn;
  int extractLines;
  double compValue;
  std::string compOperator;
  int sortEntries;
  bool beatsFirst;
  std::string joinDB;

  // besthitperset
  bool simpleBestHit;
  float alpha;
  bool shortOutput;
  int aggregationMode;

  // mergedbs
  std::string mergePrefixes;
  bool mergeStopEmpty;

  // summarizetabs
  float overlap;
  int msaType;

  // extractalignedregion
  int extractMode;

  // convertkb
  std::string kbColumns;

  // concatdbs
  bool preserveKeysB;
  bool takeLargerEntry;

  // offsetalignments
  int chainAlignment;
  int mergeQuery;

  // tsv2db
  int outputDbType;

  // diff
  bool useSequenceId;

  // prefixid
  std::string prefix;
  bool tsvOut;

  // clusterUpdate;
  bool recoverDeleted;

  // summarize headers
  int headerType;

  // filtertaxdb, filtertaxseqdb
  std::string taxonList;

  // view
  std::string idList;
  int idxEntryType;

  // lca
  int pickIdFrom;
  std::string lcaRanks;
  int showTaxLineage;
  std::string blacklist;

  // aggregatetax
  float majorityThr;
  int voteMode;

  // taxonomyreport
  int reportMode;

  // createtaxdb
  std::string ncbiTaxDump;
  std::string taxMappingFile;
  int taxMappingMode;
  int taxDbMode;

  // exapandaln
  int expansionMode;

  // taxonomy
  int taxonomySearchMode;
  int taxonomyOutputMode;

  // createsubdb
  int subDbMode;

  // tar2db
  std::string tarInclude;
  std::string tarExclude;

  // for modules that should handle -h themselves
  bool help;

  // tool citations
  std::map<unsigned int, const char *> citations;

  static Parameters &getInstance() {
    if (instance == NULL) {
      initInstance();
    }
    return *instance;
  }
  static void initInstance() { new Parameters; }

  void setDefaults();
  //    void parseParameters(int argc, const char *pargv[], const Command
  //    &command, bool printPar, int parseFlags,
  //                         int outputFlags);
  //    void printUsageMessage(const Command& command, unsigned int outputFlag,
  //    const char* extraText = NULL); void printParameters(const std::string
  //    &module, int argc, const char* pargv[],
  //                         const std::vector<BiosnakeParameter*> &par);
  //
  //    void checkIfDatabaseIsValid(const Command& command, int argc, const
  //    char** argv, bool isStartVar, bool isMiddleVar, bool isEndVar);

  std::vector<BiosnakeParameter *> removeParameter(
      const std::vector<BiosnakeParameter *> &par, const BiosnakeParameter &x);

  PARAMETER(PARAM_S)
  PARAMETER(PARAM_K)
  PARAMETER(PARAM_THREADS)
  PARAMETER(PARAM_COMPRESSED)
  PARAMETER(PARAM_ALPH_SIZE)
  PARAMETER(PARAM_MAX_SEQ_LEN)
  PARAMETER(PARAM_DIAGONAL_SCORING)
  PARAMETER(PARAM_EXACT_KMER_MATCHING)
  PARAMETER(PARAM_MASK_RESIDUES)
  PARAMETER(PARAM_MASK_LOWER_CASE)

  PARAMETER(PARAM_MIN_DIAG_SCORE)
  PARAMETER(PARAM_K_SCORE)
  PARAMETER(PARAM_MAX_SEQS)
  PARAMETER(PARAM_SPLIT)
  PARAMETER(PARAM_SPLIT_MODE)
  PARAMETER(PARAM_SPLIT_MEMORY_LIMIT)
  PARAMETER(PARAM_DISK_SPACE_LIMIT)
  PARAMETER(PARAM_SPLIT_AMINOACID)
  PARAMETER(PARAM_SUB_MAT)
  PARAMETER(PARAM_SEED_SUB_MAT)
  PARAMETER(PARAM_NO_COMP_BIAS_CORR)
  PARAMETER(PARAM_SPACED_KMER_MODE)
  PARAMETER(PARAM_REMOVE_TMP_FILES)
  PARAMETER(PARAM_INCLUDE_IDENTITY)
  PARAMETER(PARAM_PRELOAD_MODE)
  PARAMETER(PARAM_SPACED_KMER_PATTERN)
  PARAMETER(PARAM_LOCAL_TMP)
  std::vector<BiosnakeParameter *> prefilter;
  std::vector<BiosnakeParameter *> ungappedprefilter;

  // alignment
  PARAMETER(PARAM_ALIGNMENT_MODE)
  PARAMETER(PARAM_E)
  PARAMETER(PARAM_C)
  PARAMETER(PARAM_COV_MODE)
  PARAMETER(PARAM_SEQ_ID_MODE)
  PARAMETER(PARAM_MAX_REJECTED)
  PARAMETER(PARAM_MAX_ACCEPT)
  PARAMETER(PARAM_ADD_BACKTRACE)
  PARAMETER(PARAM_REALIGN)
  PARAMETER(PARAM_MIN_SEQ_ID)
  PARAMETER(PARAM_MIN_ALN_LEN)
  PARAMETER(PARAM_SCORE_BIAS)
  PARAMETER(PARAM_REALIGN_SCORE_BIAS)
  PARAMETER(PARAM_REALIGN_MAX_SEQS)
  PARAMETER(PARAM_ALT_ALIGNMENT)
  PARAMETER(PARAM_GAP_OPEN)
  PARAMETER(PARAM_GAP_EXTEND)
  PARAMETER(PARAM_ZDROP)

  // clustering
  PARAMETER(PARAM_CLUSTER_MODE)
  PARAMETER(PARAM_CLUSTER_STEPS)
  PARAMETER(PARAM_CASCADED)
  PARAMETER(PARAM_CLUSTER_REASSIGN)

  // affinity clustering
  PARAMETER(PARAM_MAXITERATIONS)
  PARAMETER(PARAM_SIMILARITYSCORE)

  // logging
  PARAMETER(PARAM_V)
  std::vector<BiosnakeParameter *> clust;

  // format alignment
  PARAMETER(PARAM_FORMAT_MODE)
  PARAMETER(PARAM_FORMAT_OUTPUT)
  PARAMETER(PARAM_DB_OUTPUT)

  // rescoremode
  PARAMETER(PARAM_RESCORE_MODE)
  PARAMETER(PARAM_WRAPPED_SCORING)
  PARAMETER(PARAM_FILTER_HITS)
  PARAMETER(PARAM_SORT_RESULTS)

  // result2msa
  PARAMETER(PARAM_MSA_FORMAT_MODE)
  PARAMETER(PARAM_ALLOW_DELETION)
  PARAMETER(PARAM_SUMMARY_PREFIX)
  PARAMETER(PARAM_SKIP_QUERY)

  // convertmsa
  PARAMETER(PARAM_IDENTIFIER_FIELD)

  // msa2profile
  PARAMETER(PARAM_MATCH_MODE)
  PARAMETER(PARAM_MATCH_RATIO)

  // result2profile
  PARAMETER(PARAM_MASK_PROFILE)
  PARAMETER(PARAM_E_PROFILE)
  PARAMETER(PARAM_FILTER_MSA)
  PARAMETER(PARAM_FILTER_MAX_SEQ_ID)
  PARAMETER(PARAM_FILTER_QSC)
  PARAMETER(PARAM_FILTER_QID)
  PARAMETER(PARAM_FILTER_COV)
  PARAMETER(PARAM_FILTER_NDIFF)
  PARAMETER(PARAM_WG)
  PARAMETER(PARAM_PCA)
  PARAMETER(PARAM_PCB)

  // sequence2profile
  PARAMETER(PARAM_NEFF)
  PARAMETER(PARAM_TAU)

  // createtsv
  PARAMETER(PARAM_TARGET_COLUMN)
  PARAMETER(PARAM_FIRST_SEQ_REP_SEQ)
  PARAMETER(PARAM_FULL_HEADER)
  PARAMETER(PARAM_IDX_SEQ_SRC)

  // result2stat
  PARAMETER(PARAM_STAT)

  // linearcluster
  PARAMETER(PARAM_KMER_PER_SEQ)
  PARAMETER(PARAM_KMER_PER_SEQ_SCALE)
  PARAMETER(PARAM_INCLUDE_ONLY_EXTENDABLE)
  PARAMETER(PARAM_IGNORE_MULTI_KMER)
  PARAMETER(PARAM_HASH_SHIFT)
  PARAMETER(PARAM_PICK_N_SIMILAR)
  PARAMETER(PARAM_ADJUST_KMER_LEN)
  PARAMETER(PARAM_RESULT_DIRECTION)
  // workflow
  PARAMETER(PARAM_RUNNER)
  PARAMETER(PARAM_REUSELATEST)

  // search workflow
  PARAMETER(PARAM_NUM_ITERATIONS)
  PARAMETER(PARAM_START_SENS)
  PARAMETER(PARAM_SENS_STEPS)
  PARAMETER(PARAM_SLICE_SEARCH)
  PARAMETER(PARAM_STRAND)
  PARAMETER(PARAM_ORF_FILTER)
  PARAMETER(PARAM_ORF_FILTER_S)
  PARAMETER(PARAM_ORF_FILTER_E)
  PARAMETER(PARAM_LCA_SEARCH)

  // easysearch
  PARAMETER(PARAM_GREEDY_BEST_HITS)

  // extractorfs
  PARAMETER(PARAM_ORF_MIN_LENGTH)
  PARAMETER(PARAM_ORF_MAX_LENGTH)
  PARAMETER(PARAM_ORF_MAX_GAP)
  PARAMETER(PARAM_CONTIG_START_MODE)
  PARAMETER(PARAM_CONTIG_END_MODE)
  PARAMETER(PARAM_ORF_START_MODE)
  PARAMETER(PARAM_ORF_FORWARD_FRAMES)
  PARAMETER(PARAM_ORF_REVERSE_FRAMES)
  PARAMETER(PARAM_USE_ALL_TABLE_STARTS)
  PARAMETER(PARAM_TRANSLATE)
  PARAMETER(PARAM_CREATE_LOOKUP)

  // indexdb
  PARAMETER(PARAM_CHECK_COMPATIBLE)
  PARAMETER(PARAM_SEARCH_TYPE)

  // createdb
  PARAMETER(PARAM_USE_HEADER)  // also used by extractorfs
  PARAMETER(PARAM_ID_OFFSET)   // same
  PARAMETER(PARAM_DB_TYPE)
  PARAMETER(PARAM_CREATEDB_MODE)
  PARAMETER(PARAM_SHUFFLE)
  PARAMETER(PARAM_WRITE_LOOKUP)

  // convert2fasta
  PARAMETER(PARAM_USE_HEADER_FILE)

  // split sequence
  PARAMETER(PARAM_SEQUENCE_OVERLAP)
  PARAMETER(PARAM_SEQUENCE_SPLIT_MODE)
  PARAMETER(PARAM_HEADER_SPLIT_MODE)

  // gff2db
  PARAMETER(PARAM_GFF_TYPE)

  // translate_nucleotide
  PARAMETER(PARAM_TRANSLATION_TABLE)
  PARAMETER(PARAM_ADD_ORF_STOP)

  // createseqfiledb
  PARAMETER(PARAM_MIN_SEQUENCES)
  PARAMETER(PARAM_MAX_SEQUENCES)
  PARAMETER(PARAM_HH_FORMAT)

  // filterDb
  PARAMETER(PARAM_FILTER_COL)
  PARAMETER(PARAM_COLUMN_TO_TAKE)
  PARAMETER(PARAM_FILTER_REGEX)
  PARAMETER(PARAM_FILTER_POS)
  PARAMETER(PARAM_FILTER_FILE)
  PARAMETER(PARAM_FILTER_EXPRESSION)
  PARAMETER(PARAM_MAPPING_FILE)
  PARAMETER(PARAM_TRIM_TO_ONE_COL)
  PARAMETER(PARAM_EXTRACT_LINES)
  PARAMETER(PARAM_COMP_OPERATOR)
  PARAMETER(PARAM_COMP_VALUE)
  PARAMETER(PARAM_SORT_ENTRIES)
  PARAMETER(PARAM_BEATS_FIRST)
  PARAMETER(PARAM_JOIN_DB)

  // besthitperset
  PARAMETER(PARAM_SIMPLE_BEST_HIT)
  PARAMETER(PARAM_ALPHA)
  PARAMETER(PARAM_SHORT_OUTPUT)
  PARAMETER(PARAM_AGGREGATION_MODE)

  // concatdb
  PARAMETER(PARAM_PRESERVEKEYS)
  PARAMETER(PARAM_TAKE_LARGER_ENTRY)

  // offsetalignment
  PARAMETER(PARAM_CHAIN_ALIGNMENT)
  PARAMETER(PARAM_MERGE_QUERY)

  // tsv2db
  PARAMETER(PARAM_OUTPUT_DBTYPE)

  // diff
  PARAMETER(PARAM_USESEQID)

  // prefixid
  PARAMETER(PARAM_PREFIX)
  PARAMETER(PARAM_TSV)

  // summarize headers
  PARAMETER(PARAM_HEADER_TYPE)

  // mergedbs
  PARAMETER(PARAM_MERGE_PREFIXES)
  PARAMETER(PARAM_MERGE_STOP_EMPTY)

  // summarizetabs
  PARAMETER(PARAM_OVERLAP)

  // extractdomains
  PARAMETER(PARAM_MSA_TYPE)

  // extract aligned region
  PARAMETER(PARAM_EXTRACT_MODE)

  // convertkb
  PARAMETER(PARAM_KB_COLUMNS)

  // clusterupdate
  PARAMETER(PARAM_RECOVER_DELETED)

  // filtertaxdb, filtertaxseqdb
  PARAMETER(PARAM_TAXON_LIST)

  // view
  PARAMETER(PARAM_ID_LIST)
  PARAMETER(PARAM_IDX_ENTRY_TYPE)

  // lca, addtaxonomy and aggregatetax
  PARAMETER(PARAM_PICK_ID_FROM)
  PARAMETER(PARAM_LCA_RANKS)
  PARAMETER(PARAM_BLACKLIST)
  PARAMETER(PARAM_TAXON_ADD_LINEAGE)

  // aggregatetax
  PARAMETER(PARAM_MAJORITY)
  PARAMETER(PARAM_VOTE_MODE)

  // taxonomyreport
  PARAMETER(PARAM_REPORT_MODE)

  // createtaxdb
  PARAMETER(PARAM_NCBI_TAX_DUMP)
  PARAMETER(PARAM_TAX_MAPPING_FILE)
  PARAMETER(PARAM_TAX_MAPPING_MODE)
  PARAMETER(PARAM_TAX_DB_MODE)

  // exapandaln
  PARAMETER(PARAM_EXPANSION_MODE)

  // taxonomy
  PARAMETER(PARAM_LCA_MODE)
  PARAMETER(PARAM_TAX_OUTPUT_MODE)

  // createsubdb
  PARAMETER(PARAM_SUBDB_MODE)

  // tar2db
  PARAMETER(PARAM_TAR_INCLUDE)
  PARAMETER(PARAM_TAR_EXCLUDE)

  // for modules that should handle -h themselves
  PARAMETER(PARAM_HELP)
  PARAMETER(PARAM_HELP_LONG)

  struct PredefinedSubstitutionMatrix {
    std::string name;
    const unsigned char *subMatData;
    unsigned int subMatDataLen;
    PredefinedSubstitutionMatrix(const char *name,
                                 const unsigned char *subMatData,
                                 const unsigned int subMatDataLen)
        : name(name), subMatData(subMatData), subMatDataLen(subMatDataLen) {}
  };
  std::vector<PredefinedSubstitutionMatrix> substitutionMatrices;

  std::vector<BiosnakeParameter *> empty;
  std::vector<BiosnakeParameter *> onlyverbosity;
  std::vector<BiosnakeParameter *> view;
  std::vector<BiosnakeParameter *> verbandcompression;
  std::vector<BiosnakeParameter *> onlythreads;
  std::vector<BiosnakeParameter *> threadsandcompression;

  std::vector<BiosnakeParameter *> alignall;
  std::vector<BiosnakeParameter *> align;
  std::vector<BiosnakeParameter *> rescorediagonal;
  std::vector<BiosnakeParameter *> alignbykmer;
  std::vector<BiosnakeParameter *> createFasta;
  std::vector<BiosnakeParameter *> convertprofiledb;
  std::vector<BiosnakeParameter *> sequence2profile;
  std::vector<BiosnakeParameter *> result2profile;
  std::vector<BiosnakeParameter *> result2pp;
  std::vector<BiosnakeParameter *> result2msa;
  std::vector<BiosnakeParameter *> result2dnamsa;
  std::vector<BiosnakeParameter *> filterresult;
  std::vector<BiosnakeParameter *> convertmsa;
  std::vector<BiosnakeParameter *> msa2profile;
  std::vector<BiosnakeParameter *> createtsv;
  std::vector<BiosnakeParameter *> result2stats;
  std::vector<BiosnakeParameter *> extractorfs;
  std::vector<BiosnakeParameter *> extractframes;
  std::vector<BiosnakeParameter *> orftocontig;
  std::vector<BiosnakeParameter *> reverseseq;
  std::vector<BiosnakeParameter *> splitdb;
  std::vector<BiosnakeParameter *> splitsequence;
  std::vector<BiosnakeParameter *> indexdb;
  std::vector<BiosnakeParameter *> kmerindexdb;
  std::vector<BiosnakeParameter *> createindex;
  std::vector<BiosnakeParameter *> createlinindex;
  std::vector<BiosnakeParameter *> convertalignments;
  std::vector<BiosnakeParameter *> createdb;
  std::vector<BiosnakeParameter *> convert2fasta;
  std::vector<BiosnakeParameter *> result2flat;
  std::vector<BiosnakeParameter *> result2repseq;
  std::vector<BiosnakeParameter *> gff2db;
  std::vector<BiosnakeParameter *> clusthash;
  std::vector<BiosnakeParameter *> kmermatcher;
  std::vector<BiosnakeParameter *> kmersearch;
  std::vector<BiosnakeParameter *> countkmer;
  std::vector<BiosnakeParameter *> easylinclustworkflow;
  std::vector<BiosnakeParameter *> linclustworkflow;
  std::vector<BiosnakeParameter *> easysearchworkflow;
  std::vector<BiosnakeParameter *> searchworkflow;
  std::vector<BiosnakeParameter *> linsearchworkflow;
  std::vector<BiosnakeParameter *> easylinsearchworkflow;
  std::vector<BiosnakeParameter *> mapworkflow;
  std::vector<BiosnakeParameter *> easyclusterworkflow;
  std::vector<BiosnakeParameter *> clusterworkflow;
  std::vector<BiosnakeParameter *> clusterUpdateSearch;
  std::vector<BiosnakeParameter *> clusterUpdateClust;
  std::vector<BiosnakeParameter *> mergeclusters;
  std::vector<BiosnakeParameter *> clusterUpdate;
  std::vector<BiosnakeParameter *> translatenucs;
  std::vector<BiosnakeParameter *> swapresult;
  std::vector<BiosnakeParameter *> swapdb;
  std::vector<BiosnakeParameter *> createseqfiledb;
  std::vector<BiosnakeParameter *> filterDb;
  std::vector<BiosnakeParameter *> offsetalignment;
  std::vector<BiosnakeParameter *> proteinaln2nucl;
  std::vector<BiosnakeParameter *> subtractdbs;
  std::vector<BiosnakeParameter *> diff;
  std::vector<BiosnakeParameter *> concatdbs;
  std::vector<BiosnakeParameter *> mergedbs;
  std::vector<BiosnakeParameter *> summarizeheaders;
  std::vector<BiosnakeParameter *> prefixid;
  std::vector<BiosnakeParameter *> summarizeresult;
  std::vector<BiosnakeParameter *> summarizetabs;
  std::vector<BiosnakeParameter *> extractdomains;
  std::vector<BiosnakeParameter *> extractalignedregion;
  std::vector<BiosnakeParameter *> convertkb;
  std::vector<BiosnakeParameter *> tsv2db;
  std::vector<BiosnakeParameter *> lca;
  std::vector<BiosnakeParameter *> majoritylca;
  std::vector<BiosnakeParameter *> addtaxonomy;
  std::vector<BiosnakeParameter *> taxonomyreport;
  std::vector<BiosnakeParameter *> filtertaxdb;
  std::vector<BiosnakeParameter *> filtertaxseqdb;
  std::vector<BiosnakeParameter *> aggregatetax;
  std::vector<BiosnakeParameter *> aggregatetaxweights;
  std::vector<BiosnakeParameter *> taxonomy;
  std::vector<BiosnakeParameter *> taxpercontig;
  std::vector<BiosnakeParameter *> easytaxonomy;
  std::vector<BiosnakeParameter *> createsubdb;
  std::vector<BiosnakeParameter *> renamedbkeys;
  std::vector<BiosnakeParameter *> createtaxdb;
  std::vector<BiosnakeParameter *> profile2pssm;
  std::vector<BiosnakeParameter *> profile2seq;
  std::vector<BiosnakeParameter *> profile2cs;
  std::vector<BiosnakeParameter *> besthitbyset;
  std::vector<BiosnakeParameter *> combinepvalbyset;
  std::vector<BiosnakeParameter *> multihitdb;
  std::vector<BiosnakeParameter *> multihitsearch;
  std::vector<BiosnakeParameter *> expandaln;
  std::vector<BiosnakeParameter *> expand2profile;
  std::vector<BiosnakeParameter *> sortresult;
  std::vector<BiosnakeParameter *> enrichworkflow;
  std::vector<BiosnakeParameter *> databases;
  std::vector<BiosnakeParameter *> tar2db;

  std::vector<DbType> databases_types;

  std::vector<BiosnakeParameter *> combineList(
      const std::vector<BiosnakeParameter *> &par1,
      const std::vector<BiosnakeParameter *> &par2);

  size_t hashParameter(biosnake_output* out, const std::vector<DbType> &dbtypes,
                       const std::vector<std::string> &filenames,
                       const std::vector<BiosnakeParameter *> &par);

  std::string createParameterString(
      biosnake_output* out,
      const std::vector<BiosnakeParameter *> &vector, bool wasSet = false);

  void overrideParameterDescription(BiosnakeParameter &par,
                                    const char *description,
                                    const char *regex = NULL, int category = 0);

  static std::vector<std::string> findMissingTaxDbFiles(
      biosnake_output* out,
      const std::string &filename);
  static void printTaxDbError(biosnake_output* out, const std::string &filename,
                              const std::vector<std::string> &missingFiles);

  static bool isEqualDbtype(const int type1, const int type2) {
    return ((type1 & 0x3FFFFFFF) == (type2 & 0x3FFFFFFF));
  }

  void setSubstitutionMatrices(std::string aminoacids, std::string nucleotides);
  void setSeedSubstitutionMatrices(std::string aminoacids,
                                   std::string nucleotides);
  void setDBFields(int no, std::string path);

  static const char *getDbTypeName(int dbtype) {
    switch (dbtype & 0x7FFFFFFF) {
      case DBTYPE_AMINO_ACIDS:
        return "Aminoacid";
      case DBTYPE_NUCLEOTIDES:
        return "Nucleotide";
      case DBTYPE_HMM_PROFILE:
        return "Profile";
      case DBTYPE_PROFILE_STATE_SEQ:
        return "Profile state";
      case DBTYPE_PROFILE_STATE_PROFILE:
        return "Profile profile";
      case DBTYPE_ALIGNMENT_RES:
        return "Alignment";
      case DBTYPE_CLUSTER_RES:
        return "Clustering";
      case DBTYPE_PREFILTER_RES:
        return "Prefilter";
      case DBTYPE_TAXONOMICAL_RESULT:
        return "Taxonomy";
      case DBTYPE_INDEX_DB:
        return "Index";
      case DBTYPE_CA3M_DB:
        return "CA3M";
      case DBTYPE_MSA_DB:
        return "MSA";
      case DBTYPE_GENERIC_DB:
        return "Generic";
      case DBTYPE_PREFILTER_REV_RES:
        return "Bi-directional prefilter";
      case DBTYPE_OFFSETDB:
        return "Offsetted headers";
      case DBTYPE_DIRECTORY:
        return "Directory";
      case DBTYPE_FLATFILE:
        return "Flatfile";
      case DBTYPE_STDIN:
        return "stdin";

      default:
        return "Unknown";
    }
  }

 protected:
  static Parameters *instance;

 public:
  Parameters();
  void operator=(Parameters const &);
  virtual ~Parameters(){};
};

struct Command {
  const char *cmd;
  int (*commandFunction)(biosnake_output *out, Parameters &par);
};

extern int align(biosnake_output *out, Parameters &par);
extern int alignall(biosnake_output *out, Parameters &par);
extern int alignbykmer(biosnake_output *out, Parameters &par);
extern int apply(biosnake_output *out, Parameters &par);
extern int besthitperset(biosnake_output *out, Parameters &par);
extern int transitivealign(biosnake_output *out, Parameters &par);
extern int clust(biosnake_output *out, Parameters &par);
extern int clusteringworkflow(biosnake_output *out, Parameters &par);
extern int clusterupdate(biosnake_output *out, Parameters &par);
extern int clusthash(biosnake_output *out, Parameters &par);
extern int combinepvalperset(biosnake_output *out, Parameters &par);
extern int compress(biosnake_output *out, Parameters &par);
extern int concatdbs(biosnake_output *out, Parameters &par);
extern int convert2fasta(biosnake_output *out, Parameters &par);
extern int convertalignments(biosnake_output *out, Parameters &par);
extern int convertca3m(biosnake_output *out, Parameters &par);
extern int convertkb(biosnake_output *out, Parameters &par);
extern int convertmsa(biosnake_output *out, Parameters &par);
extern int convertprofiledb(biosnake_output *out, Parameters &par);
extern int createdb(biosnake_output *out, Parameters &par);
extern int createindex(biosnake_output *out, Parameters &par);
extern int createlinindex(biosnake_output *out, Parameters &par);
extern int createseqfiledb(biosnake_output *out, Parameters &par);
extern int createsubdb(biosnake_output *out, Parameters &par);
extern int view(biosnake_output *out, Parameters &par);
extern int rmdb(biosnake_output *out, Parameters &par);
extern int mvdb(biosnake_output *out, Parameters &par);
extern int cpdb(biosnake_output *out, Parameters &par);
extern int lndb(biosnake_output *out, Parameters &par);
extern int createtsv(biosnake_output *out, Parameters &par);
extern int databases(biosnake_output *out, Parameters &par);
extern int dbtype(biosnake_output *out, Parameters &par);
extern int decompress(biosnake_output *out, Parameters &par);
extern int diffseqdbs(biosnake_output *out, Parameters &par);
extern int easycluster(biosnake_output *out, Parameters &par);
extern int easyrbh(biosnake_output *out, Parameters &par);
extern int easylinclust(biosnake_output *out, Parameters &par);
extern int easysearch(biosnake_output *out, Parameters &par);
extern int easylinsearch(biosnake_output *out, Parameters &par);
extern int enrich(biosnake_output *out, Parameters &par);
extern int expandaln(biosnake_output *out, Parameters &par);
extern int expand2profile(biosnake_output *out, Parameters &par);
extern int countkmer(biosnake_output *out, Parameters &par);
extern int extractalignedregion(biosnake_output *out, Parameters &par);
extern int extractdomains(biosnake_output *out, Parameters &par);
extern int extractorfs(biosnake_output *out, Parameters &par);
extern int extractframes(biosnake_output *out, Parameters &par);
extern int filterdb(biosnake_output *out, Parameters &par);
extern int filterresult(biosnake_output *out, Parameters &par);
extern int gff2db(biosnake_output *out, Parameters &par);
extern int masksequence(biosnake_output *out, Parameters &par);
extern int indexdb(biosnake_output *out, Parameters &par);
extern int kmermatcher(biosnake_output *out, Parameters &par);
extern int kmersearch(biosnake_output *out, Parameters &par);
extern int kmerindexdb(biosnake_output *out, Parameters &par);
extern int lca(biosnake_output *out, Parameters &par);
extern int lcaalign(biosnake_output *out, Parameters &par);
extern int taxonomyreport(biosnake_output *out, Parameters &par);
extern int linclust(biosnake_output *out, Parameters &par);
extern int map(biosnake_output *out, Parameters &par);
extern int renamedbkeys(biosnake_output *out, Parameters &par);
extern int majoritylca(biosnake_output *out, Parameters &par);
extern int maskbygff(biosnake_output *out, Parameters &par);
extern int mergeclusters(biosnake_output *out, Parameters &par);
extern int mergedbs(biosnake_output *out, Parameters &par);
extern int mergeresultsbyset(biosnake_output *out, Parameters &par);
extern int msa2profile(biosnake_output *out, Parameters &par);
extern int msa2result(biosnake_output *out, Parameters &par);
extern int multihitdb(biosnake_output *out, Parameters &par);
extern int multihitsearch(biosnake_output *out, Parameters &par);
extern int nrtotaxmapping(biosnake_output *out, Parameters &par);
extern int offsetalignment(biosnake_output *out, Parameters &par);
extern int orftocontig(biosnake_output *out, Parameters &par);
extern int touchdb(biosnake_output *out, Parameters &par);
extern int prefilter(biosnake_output *out, Parameters &par);
extern int prefixid(biosnake_output *out, Parameters &par);
extern int profile2cs(biosnake_output *out, Parameters &par);
extern int profile2pssm(biosnake_output *out, Parameters &par);
extern int profile2consensus(biosnake_output *out, Parameters &par);
extern int profile2repseq(biosnake_output *out, Parameters &par);
extern int proteinaln2nucl(biosnake_output *out, Parameters &par);
extern int rescorediagonal(biosnake_output *out, Parameters &par);
extern int ungappedprefilter(biosnake_output *out, Parameters &par);
extern int rbh(biosnake_output *out, Parameters &par);
extern int result2flat(biosnake_output *out, Parameters &par);
extern int result2msa(biosnake_output *out, Parameters &par);
extern int result2dnamsa(biosnake_output *out, Parameters &par);
extern int result2pp(biosnake_output *out, Parameters &par);
extern int result2profile(biosnake_output *out, Parameters &par);
extern int result2rbh(biosnake_output *out, Parameters &par);
extern int result2repseq(biosnake_output *out, Parameters &par);
extern int result2stats(biosnake_output *out, Parameters &par);
extern int reverseseq(biosnake_output *out, Parameters &par);
extern int search(biosnake_output *out, Parameters &par);
extern int linsearch(biosnake_output *out, Parameters &par);
extern int sortresult(biosnake_output *out, Parameters &par);
extern int splitdb(biosnake_output *out, Parameters &par);
extern int splitsequence(biosnake_output *out, Parameters &par);
extern int subtractdbs(biosnake_output *out, Parameters &par);
extern int suffixid(biosnake_output *out, Parameters &par);
extern int summarizeheaders(biosnake_output *out, Parameters &par);
extern int summarizeresult(biosnake_output *out, Parameters &par);
extern int summarizealis(biosnake_output *out, Parameters &par);
extern int summarizetabs(biosnake_output *out, Parameters &par);
extern int swapdb(biosnake_output *out, Parameters &par);
extern int swapresults(biosnake_output *out, Parameters &par);
extern int taxonomy(biosnake_output *out, Parameters &par);
extern int taxpercontig(biosnake_output *out, Parameters &par);
extern int easytaxonomy(biosnake_output *out, Parameters &par);
extern int createtaxdb(biosnake_output *out, Parameters &par);
extern int createbintaxonomy(biosnake_output *out, Parameters &par);
extern int translateaa(biosnake_output *out, Parameters &par);
extern int translatenucs(biosnake_output *out, Parameters &par);
extern int tsv2db(biosnake_output *out, Parameters &par);
extern int tar2db(biosnake_output *out, Parameters &par);
extern int versionstring(biosnake_output *out, Parameters &par);
extern int addtaxonomy(biosnake_output *out, Parameters &par);
extern int filtertaxdb(biosnake_output *out, Parameters &par);
extern int filtertaxseqdb(biosnake_output *out, Parameters &par);
extern int aggregatetax(biosnake_output *out, Parameters &par);
extern int aggregatetaxweights(biosnake_output *out, Parameters &par);
extern int diskspaceavail(biosnake_output *out, Parameters &par);

static std::vector<Command> commands = {
    {"convert2fasta", convert2fasta},
    {"summarizeresult", summarizeresult},
    {"convertalis", convertalignments},
    {"align", align},
    {"prefilter", prefilter},
    {"search", search},
    {"easy-search", easysearch},
    {"easy-linsearch", easylinsearch},
    {"easy-cluster", easycluster},
    {"easy-linclust", easylinclust},
    {"easy-taxonomy", easytaxonomy},
    {"easy-rbh", easyrbh},
    {"databases", databases},
    {"createdb", createdb},
    {"indexdb", indexdb},
    {"createindex", createindex},
    {"createlinindex", createlinindex},
    {"convertmsa", convertmsa},
    {"linsearch", linsearch},
    {"map", map},
    {"rbh", rbh},
    {"linclust", linclust},
    {"cluster", clusteringworkflow},
    {"clusterupdate", clusterupdate},
    {"taxonomy", taxonomy},
    {"createtsv", createtsv},
    {"createseqfiledb", createseqfiledb},
    {"createbintaxonomy", createbintaxonomy},
    {"addtaxonomy", addtaxonomy},
    {"taxonomyreport", taxonomyreport},
    {"filtertaxdb", filtertaxdb},
    {"filtertaxseqdb", filtertaxseqdb},
    {"aggregatetax", aggregatetax},
    {"aggregatetaxweights", aggregatetaxweights},
    {"lcaalign", lcaalign},
    {"lca", lca},
    {"majoritylca", majoritylca},
    {"multihitsearch", multihitsearch},
    {"besthitperset", besthitperset},
    {"combinepvalperset", combinepvalperset},
    {"mergeresultsbyset", mergeresultsbyset},
    {"kmermatcher", kmermatcher},
    {"kmersearch", kmersearch},
    {"kmerindexdb", kmerindexdb},
    {"alignall", alignall},
    {"transitivealign", transitivealign},
    {"rescorediagonal", rescorediagonal},
    {"alignbykmer", alignbykmer},
    {"clusthash", clusthash},
    {"mergeclusters", mergeclusters},
    {"decompress", decompress},
    {"rmdb", rmdb},
    {"mvdb", mvdb},
    {"cpdb", cpdb},
    {"lndb", lndb},
    {"touchdb", touchdb},
    {"concatdbs", concatdbs},
    {"splitdb", splitdb},
    {"mergedbs", mergedbs},
    {"subtractdbs", subtractdbs},
    {"apply", apply},
    {"filterdb", filterdb},
    {"swapdb", swapdb},
    {"prefixid", prefixid},
    {"suffixid", suffixid},
    {"renamedbkeys", renamedbkeys},
    {"extractframes", extractframes},
    {"orftocontig", orftocontig},
    {"reverseseq", reverseseq},
    {"translatenucs", translatenucs},
    {"translateaa", translateaa},
    {"splitsequence", splitsequence},
    {"masksequence", masksequence},
    {"extractalignedregion", extractalignedregion},
    {"filterresult", filterresult},
    {"offsetalignment", offsetalignment},
    {"sortresult", sortresult},
    {"summarizealis", summarizealis},
    {"summarizeresult", summarizeresult},
    {"convertprofiledb", convertprofiledb},
    {"expandaln", expandaln},
    {"summarizetabs", summarizetabs},
    {"maskbygff", maskbygff},
    {"convertkb", convertkb},
    {"summarizeheaders", summarizeheaders},
    {"nrtotaxmapping", nrtotaxmapping},
    {"extractdomains", extractdomains},
    {"countkmer", countkmer},
    {"version", versionstring},
    {"diskspaceavail", diskspaceavail},
};

#endif
