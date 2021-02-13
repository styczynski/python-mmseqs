#ifndef BIOSNAKE_FILEUTIL_H
#define BIOSNAKE_FILEUTIL_H

#include <cstddef>
#include <cstdio>
#include <list>
#include <string>
#include <utility>
#include <vector>
#include <biosnake/commons/parameters.h>

class FileUtil {
 public:
  static bool fileExists(biosnake_output* out, const char *fileName);

  static bool fileExistsAndIsNotEmpty(biosnake_output* out, const char *fileName);

  static bool directoryExists(biosnake_output* out, const char *directoryName);

  static FILE *openFileOrDie(biosnake_output* out, const char *fileName, const char *mode,
                             bool shouldExist);

  static size_t countLines(biosnake_output* out, const char *name);

  static bool makeDir(biosnake_output* out, const char *dirName, const int mode = 0777);

  static void deleteTempFiles(biosnake_output* out, const std::list<std::string> &tmpFiles);

  static std::string getRealPathFromSymLink(biosnake_output* out, const std::string path);

  static std::string getHashFromSymLink(biosnake_output* out, const std::string path);

  static void *mmapFile(biosnake_output* out, FILE *file, size_t *dataSize);

  static void munmapData(biosnake_output* out, void *ptr, size_t dataSize);

  static void writeFile(biosnake_output* out, const std::string &pathToFile, const unsigned char *sh,
                        size_t len);

  static std::string dirName(biosnake_output* out, const std::string &file);

  static std::string baseName(biosnake_output* out, const std::string &file);

  static size_t getFreeSpace(biosnake_output* out, const char *dir);

  static std::string getCurrentWorkingDirectory(biosnake_output* out);

  static void symlinkAlias(biosnake_output* out, const std::string &file, const std::string &alias);
  static void symlinkAbs(biosnake_output* out, const std::string &target, const std::string &link);

  static size_t getFileSize(biosnake_output* out, const std::string &fileName);

  static bool symlinkExists(biosnake_output* out, const std::string &path);

  static void copyFile(biosnake_output* out, const char *src, const char *dst);

  static FILE *openAndDelete(biosnake_output* out, const char *fileName, const char *mode);

  static std::vector<std::string> findDatafiles(biosnake_output* out, const char *datafiles);

  static void remove(biosnake_output* out, const char *file);

  static void move(biosnake_output* out, const char *src, const char *dst);

  static int parseDbType(biosnake_output* out, const char *name);

  static std::string createTemporaryDirectory(biosnake_output* out, std::string baseTmpPath,
                                              const std::string &basePath,
                                              const std::string &subDirectory);

  static void fixRlimitNoFile(biosnake_output* out);
};

#endif  // BIOSNAKE_FILEUTIL_H
