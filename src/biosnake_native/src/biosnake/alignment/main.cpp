#include <biosnake/alignment/alignment.h>
#include <biosnake/output.h>
#include <biosnake/commons/biosnakeMPI.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

#ifdef OPENMP
#include <omp.h>
#endif

int align(biosnake_output* out, Parameters& par) {
  //    BiosnakeMPI::init(argc, argv);
  //
  //    Parameters& par = Parameters::getInstance();
  //    par.overrideParameterDescription(par.PARAM_ALIGNMENT_MODE, "How to
  //    compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also
  //    start_pos and cov\n3: also seq.id", NULL, 0); par.parseParameters(argc,
  //    argv, command, true, 0, BiosnakeParameter::COMMAND_ALIGN);

  Alignment aln(out, par.db1, par.db2, par.db3, par.db3Index, par.db4, par.db4Index,
                par, false);

  out->info("Calculation of alignments");

#ifdef HAVE_MPI
  aln.run(BiosnakeMPI::rank, BiosnakeMPI::numProc);
#else
  aln.run();
#endif

  return EXIT_SUCCESS;
}

int lcaalign(biosnake_output* out, Parameters& par) {
  //    BiosnakeMPI::init(argc, argv);
  //
  //    Parameters& par = Parameters::getInstance();
  //    par.overrideParameterDescription(par.PARAM_ALIGNMENT_MODE, "How to
  //    compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also
  //    start_pos and cov\n3: also seq.id", NULL, 0); par.parseParameters(argc,
  //    argv, command, true, 0, BiosnakeParameter::COMMAND_ALIGN);

  Alignment aln(out, par.db1, par.db2, par.db3, par.db3Index, par.db4, par.db4Index,
                par, true);

#ifdef HAVE_MPI
  aln.run(BiosnakeMPI::rank, BiosnakeMPI::numProc);
#else
  aln.run();
#endif

  return EXIT_SUCCESS;
}
