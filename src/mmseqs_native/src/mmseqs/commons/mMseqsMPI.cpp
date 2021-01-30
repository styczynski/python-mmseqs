#include <mmseqs/commons/mMseqsMPI.h>
#include <mmseqs/commons/debug.h>
#include <mmseqs/commons/parameters.h>

bool MMseqsMPI::active = false;
int MMseqsMPI::rank = -1;
int MMseqsMPI::numProc = -1;

#ifdef HAVE_MPI
void MMseqsMPI::init(int argc, const char **argv) {
  MPI_Init(&argc, const_cast<char ***>(&argv));
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  active = true;

  if (!isMaster()) {
    Parameters &par = Parameters::getInstance();
    par.verbosity = Debug::ERROR;
    Debug::setDebugLevel(Debug::ERROR);
  }

  out->info("MPI Init");
  out->info("Rank: {}, Size: {}", rank, numProc);
}
#else
void MMseqsMPI::init(int, const char **) { rank = 0; }
#endif
