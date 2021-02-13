#include <biosnake/commons/biosnakeMPI.h>
#include <biosnake/output.h>
#include <biosnake/commons/parameters.h>

bool BiosnakeMPI::active = false;
int BiosnakeMPI::rank = -1;
int BiosnakeMPI::numProc = -1;

#ifdef HAVE_MPI
void BiosnakeMPI::init(int argc, const char **argv) {
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
void BiosnakeMPI::init(int, const char **) { rank = 0; }
#endif
