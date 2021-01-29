#include "alignment.h"
#include "debug.h"
#include "mMseqsMPI.h"
#include "parameters.h"
#include "util.h"
#include "output.h"

#ifdef OPENMP
#include <omp.h>
#endif

int align(mmseqs_output* out, Parameters& par) {
  //    MMseqsMPI::init(argc, argv);
  //
  //    Parameters& par = Parameters::getInstance();
  //    par.overrideParameterDescription(par.PARAM_ALIGNMENT_MODE, "How to
  //    compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also
  //    start_pos and cov\n3: also seq.id", NULL, 0); par.parseParameters(argc,
  //    argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

  Alignment aln(par.db1, par.db2, par.db3, par.db3Index, par.db4, par.db4Index,
                par, false);

  Debug(Debug::INFO) << "Calculation of alignments\n";

#ifdef HAVE_MPI
  aln.run(MMseqsMPI::rank, MMseqsMPI::numProc);
#else
  aln.run();
#endif

  return EXIT_SUCCESS;
}

int lcaalign(mmseqs_output* out, Parameters& par) {
  //    MMseqsMPI::init(argc, argv);
  //
  //    Parameters& par = Parameters::getInstance();
  //    par.overrideParameterDescription(par.PARAM_ALIGNMENT_MODE, "How to
  //    compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also
  //    start_pos and cov\n3: also seq.id", NULL, 0); par.parseParameters(argc,
  //    argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

  Alignment aln(par.db1, par.db2, par.db3, par.db3Index, par.db4, par.db4Index,
                par, true);

#ifdef HAVE_MPI
  aln.run(MMseqsMPI::rank, MMseqsMPI::numProc);
#else
  aln.run();
#endif

  return EXIT_SUCCESS;
}
