#include <mmseqs/commons/dBReader.h>
#include <mmseqs/commons/dBWriter.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/commons/parameters.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

#include <climits>

int createsubdb(mmseqs_output *out, Parameters &par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  FILE *orderFile = NULL;
  if (FileUtil::fileExists(par.db1Index.c_str())) {
    orderFile = fopen(par.db1Index.c_str(), "r");
  } else {
    if (FileUtil::fileExists(par.db1.c_str())) {
      orderFile = fopen(par.db1.c_str(), "r");
    } else {
      out->failure("File {} does not exist", par.db1 );
    }
  }

  DBReader<unsigned int> reader(
      par.db2.c_str(), par.db2Index.c_str(), 1,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  reader.open(DBReader<unsigned int>::NOSORT);
  const bool isCompressed = reader.isCompressed();

  DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, 0,
                  Parameters::DBTYPE_OMIT_FILE);
  writer.open();
  // getline reallocs automatic
  char *line = NULL;
  size_t len = 0;
  char dbKey[256];
  unsigned int prevKey = 0;
  bool isOrdered = true;
  while (getline(&line, &len, orderFile) != -1) {
    Util::parseKey(line, dbKey);
    const unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
    isOrdered &= (prevKey <= key);
    prevKey = key;
    const size_t id = reader.getId(key);
    if (id >= UINT_MAX) {
      out->warn("Key {} not found in database", dbKey);
      continue;
    }
    if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
      writer.writeIndexEntry(key, reader.getOffset(id), reader.getEntryLen(id),
                             0);
    } else {
      char *data = reader.getDataUncompressed(id);
      size_t originalLength = reader.getEntryLen(id);
      size_t entryLength = std::max(originalLength, static_cast<size_t>(1)) - 1;

      if (isCompressed) {
        // copy also the null byte since it contains the information if
        // compressed or not
        entryLength = *(reinterpret_cast<unsigned int *>(data)) +
                      sizeof(unsigned int) + 1;
        writer.writeData(data, entryLength, key, 0, false, false);
      } else {
        writer.writeData(data, entryLength, key, 0, true, false);
      }
      // do not write null byte since
      writer.writeIndexEntry(key, writer.getStart(0), originalLength, 0);
    }
  }
  // merge any kind of sequence database
  const bool shouldMerge =
      Parameters::isEqualDbtype(reader.getDbtype(),
                                Parameters::DBTYPE_HMM_PROFILE) ||
      Parameters::isEqualDbtype(reader.getDbtype(),
                                Parameters::DBTYPE_AMINO_ACIDS) ||
      Parameters::isEqualDbtype(reader.getDbtype(),
                                Parameters::DBTYPE_NUCLEOTIDES) ||
      Parameters::isEqualDbtype(reader.getDbtype(),
                                Parameters::DBTYPE_PROFILE_STATE_PROFILE) ||
      Parameters::isEqualDbtype(reader.getDbtype(),
                                Parameters::DBTYPE_PROFILE_STATE_SEQ);
  writer.close(shouldMerge, !isOrdered);
  if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
    DBReader<unsigned int>::softlinkDb(par.db2, par.db3, DBFiles::DATA);
  }
  DBWriter::writeDbtypeFile(par.db3.c_str(), reader.getDbtype(), isCompressed);
  DBReader<unsigned int>::softlinkDb(par.db2, par.db3,
                                     DBFiles::SEQUENCE_ANCILLARY);

  free(line);
  reader.close();
  if (fclose(orderFile) != 0) {
    out->error("Cannot close file {}", par.db1);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
