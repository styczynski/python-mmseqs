#ifndef MMSEQS_INDEXBUILDER_H
#define MMSEQS_INDEXBUILDER_H

#include <mmseqs/prefiltering/indexTable.h>
#include <mmseqs/output.h>

class IndexBuilder {
 public:
  static void fillDatabase(mmseqs_output* out, IndexTable *indexTable,
                           SequenceLookup **maskedLookup,
                           SequenceLookup **unmaskedLookup, BaseMatrix &subMat,
                           Sequence *seq, DBReader<unsigned int> *dbr,
                           size_t dbFrom, size_t dbTo, int kmerThr, bool mask,
                           bool maskLowerCaseMode);
};

#endif
