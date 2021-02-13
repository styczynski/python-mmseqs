#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/taxonomy/ncbiTaxonomy.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/output.h>

int createbintaxonomy(biosnake_output* out, Parameters& par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, false, 0, 0);
  NcbiTaxonomy taxonomy(out, par.db1, par.db2, par.db3);
  std::pair<char*, size_t> serialized = NcbiTaxonomy::serialize(taxonomy);
  FILE* handle = fopen(par.db4.c_str(), "w");
  if (handle == NULL) {
    out->failure("Could not open {} for writing", par.db4 );
  }
  fwrite(serialized.first, serialized.second, sizeof(char), handle);
  fclose(handle);
  free(serialized.first);
}
