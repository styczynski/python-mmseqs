#include <biosnake/alignment/multipleAlignment.h>

#include <biosnake/output.h>
#include <biosnake/commons/sequence.h>
#include <biosnake/commons/substitutionMatrix.h>
#include <biosnake/commons/util.h>

MultipleAlignment::MultipleAlignment(biosnake_output* output, size_t maxSeqLen,
                                     SubstitutionMatrix *subMat): out(output) {
  this->maxSeqLen = maxSeqLen;
  this->maxMsaSeqLen = maxSeqLen * 2;

  this->subMat = subMat;
  this->queryGaps = new unsigned int[maxMsaSeqLen];
}

char *MultipleAlignment::initX(int len) {
  int seqSimdLength = (len) / (VECSIZE_INT * 4) + 2;
  seqSimdLength *= (VECSIZE_INT * 4);
  char *ptr = (char *)malloc_simd_int(seqSimdLength);
  std::fill(ptr, ptr + seqSimdLength, MultipleAlignment::GAP);
  return ptr;
}

void MultipleAlignment::deleteMSA(MultipleAlignment::MSAResult *res) {
  for (size_t i = 0; i < res->setSize; i++) {
    free(res->msaSequence[i]);
  }
  delete[] res->msaSequence;
}

MultipleAlignment::~MultipleAlignment() { delete[] queryGaps; }

void MultipleAlignment::print(MSAResult msaResult, SubstitutionMatrix *subMat) {
  for (size_t i = 0; i < msaResult.setSize; i++) {
    for (size_t pos = 0; pos < msaResult.msaSequenceLength; pos++) {
      char aa = msaResult.msaSequence[i][pos];
      printf("%c", (aa < NAA) ? subMat->num2aa[(int)aa] : '-');
    }
    printf("\n");
  }
}

void MultipleAlignment::computeQueryGaps(
    unsigned int *queryGaps, Sequence *centerSeq, size_t edges,
    const std::vector<Matcher::result_t> &alignmentResults) {
  // init query gaps
  memset(queryGaps, 0, sizeof(unsigned int) * centerSeq->L);
  for (size_t i = 0; i < edges; i++) {
    const Matcher::result_t &alignment = alignmentResults[i];
    const std::string &bt = alignment.backtrace;
    size_t queryPos = 0;
    size_t targetPos = 0;
    size_t currentQueryGapSize = 0;
    queryPos = alignment.qStartPos;
    targetPos = alignment.dbStartPos;
    // compute query gaps (deletions)
    for (size_t pos = 0; pos < bt.size(); ++pos) {
      char bt_letter = bt.at(pos);
      if (bt_letter == 'M') {  // match state
        ++queryPos;
        ++targetPos;
        currentQueryGapSize = 0;
      } else {
        if (bt_letter == 'I') {  // insertion
          ++queryPos;
          currentQueryGapSize = 0;
        } else {  // deletion
          ++targetPos;
          currentQueryGapSize += 1;
          size_t gapCount = queryGaps[queryPos];
          queryGaps[queryPos] = std::max(gapCount, currentQueryGapSize);
        }
      }
    }
  }
}

size_t MultipleAlignment::updateGapsInCenterSequence(char **msaSequence,
                                                     Sequence *centerSeq,
                                                     bool noDeletionMSA) {
  size_t centerSeqPos = 0;
  for (int queryPos = 0; queryPos < centerSeq->L; queryPos++) {
    if (centerSeqPos >= maxMsaSeqLen) {
      out->failure("queryMSASize ({}) is >= maxMsaSeqLen ({})", centerSeqPos, maxMsaSeqLen);
    }
    for (size_t gapIdx = 0; gapIdx < queryGaps[queryPos]; gapIdx++) {
      if (noDeletionMSA == false) {
        msaSequence[0][centerSeqPos] = '-';
        centerSeqPos++;
      }
    }
    msaSequence[0][centerSeqPos] =
        subMat->num2aa[centerSeq->numSequence[queryPos]];
    centerSeqPos++;
  }
  return centerSeqPos;
}

void MultipleAlignment::updateGapsInSequenceSet(
    char **msaSequence, size_t centerSeqSize,
    const std::vector<std::vector<unsigned char>> &seqs,
    const std::vector<Matcher::result_t> &alignmentResults,
    unsigned int *queryGaps, bool noDeletionMSA) {
  for (size_t i = 0; i < seqs.size(); i++) {
    const Matcher::result_t &result = alignmentResults[i];
    const std::string &bt = result.backtrace;
    char *edgeSeqMSA = msaSequence[i + 1];
    const std::vector<unsigned char> &edgeSeq = seqs[i];
    unsigned int queryPos = result.qStartPos;
    unsigned int targetPos = result.dbStartPos;
    // HACK: score was 0 and sequence was rejected, so we fill in an empty gap
    // sequence
    if (targetPos == UINT_MAX) {
      out->warn("Edge sequence {} was not aligned.", i);
      // fill up with gaps
      for (size_t pos = 0; pos < centerSeqSize; pos++) {
        edgeSeqMSA[pos] = '-';
      }
      continue;
    }
    size_t bufferPos = 0;
    // fill initial positions with gaps (local alignment)
    for (int pos = 0; pos < result.qStartPos; pos++) {
      edgeSeqMSA[bufferPos] = '-';
      bufferPos++;
    }
    for (size_t alnPos = 0; alnPos < bt.size(); alnPos++) {
      if (bufferPos >= maxMsaSeqLen) {
        out->failure("BufferPos ({}) is >= maxMsaSeqLen ({})", bufferPos, maxMsaSeqLen);
      }
      if (bt.at(alnPos) == 'I') {
        edgeSeqMSA[bufferPos] = '-';
        bufferPos++;
        queryPos++;
      } else {
        // D state in target Sequence
        if (bt.at(alnPos) == 'D') {
          while (bt.at(alnPos) == 'D' && alnPos < bt.size()) {
            if (noDeletionMSA == false) {
              edgeSeqMSA[bufferPos] = subMat->num2aa[edgeSeq[targetPos]];
              bufferPos++;
            }
            targetPos++;
            alnPos++;
          }
          if (alnPos >= bt.size()) {
            break;
          } else if (bt.at(alnPos) == 'I') {
            edgeSeqMSA[bufferPos] = '-';
            bufferPos++;
            queryPos++;
          } else if (bt.at(alnPos) == 'M') {
            edgeSeqMSA[bufferPos] = subMat->num2aa[edgeSeq[targetPos]];
            bufferPos++;
            queryPos++;
            targetPos++;
          }
          continue;
        } else if (bt.at(alnPos) == 'M') {
          // add query deletion gaps
          for (size_t gapIdx = 0; gapIdx < queryGaps[queryPos]; gapIdx++) {
            if (noDeletionMSA == false) {
              edgeSeqMSA[bufferPos] = '-';
              bufferPos++;
            }
          }
          // M state
          edgeSeqMSA[bufferPos] = subMat->num2aa[edgeSeq[targetPos]];

          bufferPos++;
          queryPos++;
          targetPos++;
        }
      }
    }
    // fill up rest with gaps
    for (size_t pos = bufferPos; pos < centerSeqSize; pos++) {
      edgeSeqMSA[bufferPos] = '-';
      bufferPos++;
    }
  }
}

MultipleAlignment::MSAResult MultipleAlignment::computeMSA(
    Sequence *centerSeq,
    const std::vector<std::vector<unsigned char>> &edgeSeqs,
    const std::vector<Matcher::result_t> &alignmentResults,
    bool noDeletionMSA) {
  if (edgeSeqs.empty()) {
    return singleSequenceMSA(centerSeq);
  }

  if (edgeSeqs.size() != alignmentResults.size()) {
    out->failure("edgeSeqs.size ({}) is != alignmentResults.size ({})", edgeSeqs.size(), alignmentResults.size());
  }

  char **msaSequence = new char *[edgeSeqs.size() + 1];
  for (size_t i = 0; i <= edgeSeqs.size(); i++) {
    // FIXME: in deletion case, the msa could become even larger than maxSeqLen
    msaSequence[i] = initX(noDeletionMSA ? centerSeq->L + 1 : maxSeqLen + 1);
  }

  computeQueryGaps(queryGaps, centerSeq, edgeSeqs.size(), alignmentResults);
  // process gaps in Query (update sequences)
  // and write query Alignment at position 0

  size_t centerSeqSize =
      updateGapsInCenterSequence(msaSequence, centerSeq, noDeletionMSA);

  // compute the MSA alignment
  updateGapsInSequenceSet(msaSequence, centerSeqSize, edgeSeqs,
                          alignmentResults, queryGaps, noDeletionMSA);

  // clean vector
  // alignmentResults.clear();
  // map to int
  for (size_t k = 0; k < edgeSeqs.size() + 1; ++k) {
    for (size_t pos = 0; pos < centerSeqSize; ++pos) {
      msaSequence[k][pos] =
          (msaSequence[k][pos] == '-')
              ? GAP
              : static_cast<int>(
                    subMat->aa2num[static_cast<int>(msaSequence[k][pos])]);
    }
    int len = std::min(maxMsaSeqLen, (centerSeqSize + VECSIZE_INT * 4));
    int startPos = std::min(centerSeqSize, maxMsaSeqLen - 1);
    for (int pos = startPos; pos < len; pos++) {
      msaSequence[k][pos] = GAP;
    }
  }

  // +1 for the query
  return MSAResult(centerSeqSize, centerSeq->L, edgeSeqs.size() + 1,
                   msaSequence);
}

MultipleAlignment::MSAResult MultipleAlignment::singleSequenceMSA(
    Sequence *centerSeq) {
  size_t queryMSASize = 0;
  char **msaSequence = new char *[1];
  msaSequence[0] = initX(centerSeq->L);
  for (int queryPos = 0; queryPos < centerSeq->L; queryPos++) {
    if (queryMSASize >= maxMsaSeqLen) {
      out->failure("queryMSASize ({}) is >= maxMsaSeqLen ({})", queryMSASize, maxMsaSeqLen);
    }
    msaSequence[0][queryMSASize] = (char)centerSeq->numSequence[queryPos];
    queryMSASize++;
  }
  return MSAResult(queryMSASize, centerSeq->L, 1, msaSequence);
}
