#include <fstream>
#include <sstream>

#include <biosnake/commons/dBWriter.h>
#include <biosnake/output.h>
#include <biosnake/commons/parameters.h>
#include <biosnake/commons/util.h>
#include <biosnake/output.h>

int tsv2db(biosnake_output *out, Parameters &par) {
  //    Parameters &par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  out->info("Output database type: {}", Parameters::getDbTypeName(par.outputDbType));
  if (par.PARAM_OUTPUT_DBTYPE.wasSet == false) {
    out->info("Consider setting --output-dbtype.");
  }

  DBWriter writer(out, par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed,
                  par.outputDbType);
  writer.open();

  std::ifstream tsv(par.db1);
  if (tsv.fail()) {
    out->failure("File {} not found", par.db1);
  }

  std::ostringstream ss;
  char keyData[255];
  bool skippedFirst = false;
  std::string lastKey;
  std::string line;
  while (std::getline(tsv, line)) {
    char *current = (char *)line.c_str();
    Util::parseKey(current, keyData);
    const std::string key(keyData);

    if (key != lastKey && skippedFirst == true) {
      if (par.includeIdentity) {
        const std::string temp = ss.str();
        ss.seekp(0);
        ss << lastKey << "\n";
        ss << temp;
      }
      const std::string result = ss.str();
      unsigned int keyId = strtoull(lastKey.c_str(), NULL, 10);
      writer.writeData(result.c_str(), result.length(), keyId);
      ss.str("");
      ss.clear();
    }

    char *restStart = current + key.length();
    restStart = restStart + Util::skipWhitespace(restStart);
    char *restEnd = restStart;
    restEnd = Util::seekToNextEntry(restEnd) - 1;

    const std::string rest(restStart, restEnd - restStart);

    skippedFirst = true;
    ss << rest << "\n";
    lastKey = key;
  }

  if (par.includeIdentity) {
    const std::string temp = ss.str();
    ss.seekp(0);
    ss << lastKey << "\n";
    ss << temp;
  }
  const std::string result = ss.str();
  unsigned int keyId = strtoull(lastKey.c_str(), NULL, 10);
  writer.writeData(result.c_str(), result.length(), keyId);

  writer.close();

  return EXIT_SUCCESS;
}
