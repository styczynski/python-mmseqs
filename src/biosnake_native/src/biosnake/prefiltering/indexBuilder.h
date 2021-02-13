#ifndef BIOSNAKE_INDEXBUILDER_H
#define BIOSNAKE_INDEXBUILDER_H

#include <biosnake/prefiltering/indexTable.h>
#include <biosnake/output.h>

class IndexBuilder {
 public:
  static void fillDatabase(biosnake_output* out, IndexTable *indexTable,
                           SequenceLookup **maskedLookup,
                           SequenceLookup **unmaskedLookup, BaseMatrix &subMat,
                           Sequence *seq, DBReader<unsigned int> *dbr,
                           size_t dbFrom, size_t dbTo, int kmerThr, bool mask,
                           bool maskLowerCaseMode);
};

#endif
