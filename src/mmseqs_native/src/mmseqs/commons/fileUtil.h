#ifndef MMSEQS_FILEUTIL_H
#define MMSEQS_FILEUTIL_H

#include <cstddef>
#include <cstdio>
#include <list>
#include <string>
#include <utility>
#include <vector>
#include <mmseqs/commons/parameters.h>

class FileUtil {
 public:
  static bool fileExists(mmseqs_output* out, const char *fileName);

  static bool fileExistsAndIsNotEmpty(mmseqs_output* out, const char *fileName);

  static bool directoryExists(mmseqs_output* out, const char *directoryName);

  static FILE *openFileOrDie(mmseqs_output* out, const char *fileName, const char *mode,
                             bool shouldExist);

  static size_t countLines(mmseqs_output* out, const char *name);

  static bool makeDir(mmseqs_output* out, const char *dirName, const int mode = 0777);

  static void deleteTempFiles(mmseqs_output* out, const std::list<std::string> &tmpFiles);

  static std::string getRealPathFromSymLink(mmseqs_output* out, const std::string path);

  static std::string getHashFromSymLink(mmseqs_output* out, const std::string path);

  static void *mmapFile(mmseqs_output* out, FILE *file, size_t *dataSize);

  static void munmapData(mmseqs_output* out, void *ptr, size_t dataSize);

  static void writeFile(mmseqs_output* out, const std::string &pathToFile, const unsigned char *sh,
                        size_t len);

  static std::string dirName(mmseqs_output* out, const std::string &file);

  static std::string baseName(mmseqs_output* out, const std::string &file);

  static size_t getFreeSpace(mmseqs_output* out, const char *dir);

  static std::string getCurrentWorkingDirectory(mmseqs_output* out);

  static void symlinkAlias(mmseqs_output* out, const std::string &file, const std::string &alias);
  static void symlinkAbs(mmseqs_output* out, const std::string &target, const std::string &link);

  static size_t getFileSize(mmseqs_output* out, const std::string &fileName);

  static bool symlinkExists(mmseqs_output* out, const std::string &path);

  static void copyFile(mmseqs_output* out, const char *src, const char *dst);

  static FILE *openAndDelete(mmseqs_output* out, const char *fileName, const char *mode);

  static std::vector<std::string> findDatafiles(mmseqs_output* out, const char *datafiles);

  static void remove(mmseqs_output* out, const char *file);

  static void move(mmseqs_output* out, const char *src, const char *dst);

  static int parseDbType(mmseqs_output* out, const char *name);

  static std::string createTemporaryDirectory(mmseqs_output* out, std::string baseTmpPath,
                                              const std::string &basePath,
                                              const std::string &subDirectory);

  static void fixRlimitNoFile(mmseqs_output* out);
};

#endif  // MMSEQS_FILEUTIL_H
