//
// Written by Martin Steinegger
//
// Wrapper for KSW2 aligner.
// Local banded nucleotide aligner
//
#include <biosnake/commons/nucleotideMatrix.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/alignment/stripedSmithWaterman.h>
#include <biosnake/output.h>
#include <biosnake/commons/substitutionMatrix.h>
#include <biosnake/commons/util.h>

class BandedNucleotideAligner {
 public:
  BandedNucleotideAligner(biosnake_output* output, BaseMatrix *subMat, size_t maxSequenceLength,
                          int gapo, int gape, int zdrop);

  ~BandedNucleotideAligner();

  void initQuery(Sequence *q);

  s_align align(Sequence *targetSeqObj, int diagonal, bool reverse,
                std::string &backtrace, int &aaIds, EvalueComputation *evaluer,
                bool wrappedScoring = false);

 private:
  biosnake_output* out;
  SubstitutionMatrix::FastMatrix fastMatrix;
  uint8_t *targetSeqRev;
  int targetSeqRevDataLen;
  uint8_t *querySeq;
  uint8_t *querySeqRev;
  int querySeqRevDataLen;
  uint8_t *queryRevCompSeq;
  char *queryRevCompCharSeq;
  uint8_t *queryRevCompSeqRev;
  Sequence *querySeqObj;
  int8_t *mat;
  NucleotideMatrix *subMat;
  //    uint32_t * cigar;
  int gapo;
  int gape;
  int zdrop;
};
