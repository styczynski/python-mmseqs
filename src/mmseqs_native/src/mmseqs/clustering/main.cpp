#include <mmseqs/clustering/clustering.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/output.h>

int clust(mmseqs_output* out, Parameters& par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  Clustering clu(par.db1, par.db1Index, par.db2, par.db2Index, par.db3,
                 par.db3Index, par.maxIteration, par.similarityScoreType,
                 par.threads, par.compressed);
  clu.run(par.clusteringMode);
  return EXIT_SUCCESS;
}
