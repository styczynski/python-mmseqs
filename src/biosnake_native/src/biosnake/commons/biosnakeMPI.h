#ifndef BIOSNAKE_MPI_H
#define BIOSNAKE_MPI_H

#ifdef HAVE_MPI
#include <mpi.h>
#endif

class BiosnakeMPI {
 public:
  static const int MASTER = 0;

  static bool active;
  static int rank;
  static int numProc;

  static void init(int argc, const char **argv);
  static inline bool isMaster() {
#ifdef HAVE_MPI
    return rank == MASTER;
#else
    return true;
#endif
  };
};

#endif