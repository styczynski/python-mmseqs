#include <climits>
#include <fstream>
#include <string>

#include <mmseqs/commons/dBReader.h>
#include <mmseqs/commons/dBWriter.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/util.h>
#include <mmseqs/output.h>

int maskbygff(mmseqs_output* out, Parameters& par) {
  //    Parameters& par = Parameters::getInstance();
  //    par.parseParameters(argc, argv, command, true, 0, 0);

  DBReader<std::string> reader(
      out, par.db2.c_str(), par.db2Index.c_str(), par.threads,
      DBReader<std::string>::USE_DATA | DBReader<std::string>::USE_WRITABLE);
  reader.open(DBReader<std::string>::NOSORT);

  bool shouldCompareType = par.gffType.length() > 0;

  size_t entries_num = 0;

  std::ifstream gffFile(par.db1);
  std::string gffLine;
  while (std::getline(gffFile, gffLine)) {
    entries_num++;

    // line is a comment
    if (gffLine.find_first_of("#", 0, 1) != std::string::npos) {
      continue;
    }

    std::vector<std::string> fields = Util::split(gffLine, "\t");

    // gff has always 9 fields
    if (fields.size() != 9) {
      out->warn("Invalid GFF format in line {}", entries_num);
      continue;
    }

    std::string name(fields[0]);

    std::string type(fields[2]);

    if (shouldCompareType && type.compare(par.gffType) != 0) {
      continue;
    }

    char* rest;
    size_t start = strtoull(fields[3].c_str(), &rest, 10);
    if ((rest != fields[3].c_str() && *rest != '\0') || errno == ERANGE) {
      out->warn("Invalid start position format in line {}", entries_num);
      continue;
    }
    size_t end = strtoull(fields[4].c_str(), &rest, 10);
    if ((rest != fields[4].c_str() && *rest != '\0') || errno == ERANGE) {
      out->warn("Invalid end position format in line {}", entries_num);
      continue;
    }

    // gff start and end are supposed to be 1-indexed
    if (end <= start || end == 0 || start == 0) {
      out->warn("Invalid sequence length in line {}", entries_num);
      continue;
    }

    // so 0-index them now
    start -= 1;
    end -= 1;

    size_t id = reader.getId(name);
    if (id == UINT_MAX) {
      out->error("GFF entry not found in input database: {}", name);

      return EXIT_FAILURE;
    }

    char* body = reader.getData(id, 0);
    for (char* i = body + start; i <= body + end; ++i) {
      *i = 'X';
    }
  }
  gffFile.close();

  DBWriter writer(out, par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed,
                  reader.getDbtype());
  writer.open();

  DBReader<std::string> headerReader(
      out, par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads,
      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
  headerReader.open(DBReader<std::string>::NOSORT);

  DBWriter headerWriter(out, par.hdr3.c_str(), par.hdr3Index.c_str(), 1,
                        par.compressed, Parameters::DBTYPE_GENERIC_DB);
  headerWriter.open();

  for (size_t i = 0; i < reader.getSize(); ++i) {
    unsigned int id = par.identifierOffset + i;

    // ignore nulls
    writer.writeData(reader.getData(i, 0), reader.getEntryLen(i) - 1, id);
    headerWriter.writeData(headerReader.getData(i, 0),
                           headerReader.getEntryLen(i) - 1, id);
  }
  headerWriter.close();
  writer.close();
  headerReader.close();
  reader.close();

  return EXIT_SUCCESS;
}
