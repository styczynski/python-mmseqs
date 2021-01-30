/*
 * convert2fasta
 * written by Milot Mirdita <milot@mirdita.de>
 */

#include <cstdio>
#include <cstring>

#include <mmseqs/commons/dBReader.h>
#include <mmseqs/commons/debug.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

const char headerStart[] = {'>'};
const char newline[] = {'\n'};

int convert2fasta(mmseqs_output* out, Parameters& par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  DBReader<unsigned int> db(
      par.db1.c_str(), par.db1Index.c_str(), 1,
      DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
  db.open(DBReader<unsigned int>::NOSORT);

  DBReader<unsigned int> db_header(
      par.hdr1.c_str(), par.hdr1Index.c_str(), 1,
      DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
  db_header.open(DBReader<unsigned int>::NOSORT);

  FILE* fastaFP = fopen(par.db2.c_str(), "w");
  if (fastaFP == NULL) {
    perror(par.db2.c_str());
    EXIT(EXIT_FAILURE);
  }

  DBReader<unsigned int>* from = &db;
  if (par.useHeaderFile) {
    from = &db_header;
  }

  out->info("Start writing file to {}\n", par.db2);
  for (size_t i = 0; i < from->getSize(); i++) {
    unsigned int key = from->getDbKey(i);
    unsigned int headerKey = db_header.getId(key);
    const char* headerData = db_header.getData(headerKey, 0);
    const size_t headerLen = db_header.getEntryLen(headerKey);

    fwrite(headerStart, sizeof(char), 1, fastaFP);
    fwrite(headerData, sizeof(char), headerLen - 2, fastaFP);
    fwrite(newline, sizeof(char), 1, fastaFP);

    unsigned int bodyKey = db.getId(key);
    const char* bodyData = db.getData(bodyKey, 0);
    const size_t bodyLen = db.getEntryLen(bodyKey);
    fwrite(bodyData, sizeof(char), bodyLen - 2, fastaFP);
    fwrite(newline, sizeof(char), 1, fastaFP);
  }
  if (fclose(fastaFP) != 0) {
    Debug(Debug::ERROR) << "Cannot close file " << par.db2 << "\n";
    EXIT(EXIT_FAILURE);
  }
  db_header.close();
  db.close();

  return EXIT_SUCCESS;
}
