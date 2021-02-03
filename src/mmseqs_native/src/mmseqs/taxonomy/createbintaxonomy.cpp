#include <mmseqs/output.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/taxonomy/ncbiTaxonomy.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/output.h>

int createbintaxonomy(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, false, 0, 0);
  NcbiTaxonomy taxonomy(out, par.db1, par.db2, par.db3);
  std::pair<char*, size_t> serialized = NcbiTaxonomy::serialize(taxonomy);
  FILE* handle = fopen(par.db4.c_str(), "w");
  if (handle == NULL) {
    out->failure("Could not open {} for writing\n", par.db4 );
  }
  fwrite(serialized.first, serialized.second, sizeof(char), handle);
  fclose(handle);
  free(serialized.first);
  EXIT(EXIT_SUCCESS);
}
