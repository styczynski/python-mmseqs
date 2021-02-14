//
// Created by lars on 10.06.15.
//

#ifndef BIOSNAKE_ALIGNMENTSYMMETRY_H
#define BIOSNAKE_ALIGNMENTSYMMETRY_H
#include <biosnake/output.h>
#include <biosnake/commons/util.h>
#include <list>
#include <set>

#include <biosnake/commons/dBReader.h>

class AlignmentSymmetry {
 public:
  static void readInData(biosnake_output* out, DBReader<unsigned int> *pReader,
                         DBReader<unsigned int> *pDBReader, unsigned int **pInt,
                         unsigned short **elementScoreTable, int scoretype,
                         size_t *offsets);
  template <typename T>
  static void computeOffsetFromCounts(biosnake_output* out, T *elementSizes, size_t dbSize) {
    size_t prevElementLength = elementSizes[0];
    elementSizes[0] = 0;
    for (size_t i = 0; i < dbSize; i++) {
      const size_t currElementLength = elementSizes[i + 1];
      elementSizes[i + 1] = elementSizes[i] + prevElementLength;
      prevElementLength = currElementLength;
    }
  }
  static size_t findMissingLinks(biosnake_output* out, unsigned int **elementLookupTable,
                                 size_t *offsetTable, size_t dbSize,
                                 int threads);
  static void addMissingLinks(biosnake_output* out, unsigned int **elementLookupTable,
                              size_t *offsetTable, size_t *newOffset,
                              size_t dbSize,
                              unsigned short **elementScoreTable);
  static void sortElements(biosnake_output* out, unsigned int **elementLookupTable, size_t *offsets,
                           size_t dbSize);

  template <typename T>
  static void setupPointers(biosnake_output* out, T *elements, T **elementLookupTable,
                            size_t *elementOffset, unsigned int dbSize,
                            size_t totalElementCount) {
    for (size_t i = 0; i < dbSize; i++) {
      if (totalElementCount < elementOffset[i]) {
        out->failure("Error in setupPointers. totalElementCount ({}) < elementOffset[{}] ({})", totalElementCount, i, elementOffset[i]);
      }
      elementLookupTable[i] = elements + elementOffset[i];
    }
  }
};
#endif  // BIOSNAKE_ALIGNMENTSYMMETRY_H
