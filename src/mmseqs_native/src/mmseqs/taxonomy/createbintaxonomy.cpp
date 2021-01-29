#include <mmseqs/commons/debug.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/taxonomy/ncbiTaxonomy.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/output.h>

int createbintaxonomy(mmseqs_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, false, 0, 0);
  NcbiTaxonomy taxonomy(par.db1, par.db2, par.db3);
  std::pair<char*, size_t> serialized = NcbiTaxonomy::serialize(taxonomy);
  FILE* handle = fopen(par.db4.c_str(), "w");
  if (handle == NULL) {
    Debug(Debug::ERROR) << "Could not open " << par.db4 << " for writing\n";
    EXIT(EXIT_FAILURE);
  }
  fwrite(serialized.first, serialized.second, sizeof(char), handle);
  fclose(handle);
  free(serialized.first);
  EXIT(EXIT_SUCCESS);
}
