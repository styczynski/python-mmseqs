//
// Written by Martin Steinegger
//
// Wrapper for KSW2 aligner.
// Local banded nucleotide aligner
//
#include <mmseqs/commons/nucleotideMatrix.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/alignment/stripedSmithWaterman.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/substitutionMatrix.h>
#include <mmseqs/commons/util.h>

class BandedNucleotideAligner {
 public:
  BandedNucleotideAligner(mmseqs_output* output, BaseMatrix *subMat, size_t maxSequenceLength,
                          int gapo, int gape, int zdrop);

  ~BandedNucleotideAligner();

  void initQuery(Sequence *q);

  s_align align(Sequence *targetSeqObj, int diagonal, bool reverse,
                std::string &backtrace, int &aaIds, EvalueComputation *evaluer,
                bool wrappedScoring = false);

 private:
  mmseqs_output* out;
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
