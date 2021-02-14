#include <biosnake/output.h>
#include <biosnake/commons/indexReader.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

int view(biosnake_output* out, Parameters& par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, false, 0, 0);
  std::vector<std::string> ids = Util::split(par.idList, ",");
  int indexSrcType = IndexReader::SEQUENCES;
  switch (par.idxEntryType) {
    case 0:
      indexSrcType = IndexReader::SEQUENCES;
      break;
    case 1:
      indexSrcType = IndexReader::SRC_SEQUENCES;
      break;
    case 2:
      indexSrcType = IndexReader::HEADERS;
      break;
    case 3:
      indexSrcType = IndexReader::SRC_HEADERS;
      break;
  }
  IndexReader reader(out, par.db1, par.threads, indexSrcType, 0);
  for (size_t i = 0; i < ids.size(); ++i) {
    const unsigned int key = Util::fast_atoi<unsigned int>(ids[i].c_str());
    const size_t id = reader.sequenceReader->getId(key);
    if (id >= UINT_MAX) {
      out->error("Key {} not found in database", ids[i]);
      continue;
    }
    char* data = reader.sequenceReader->getData(id, 0);
    size_t size = reader.sequenceReader->getEntryLen(id) - 1;
    fwrite(data, sizeof(char), size, stdout);
  }
}
