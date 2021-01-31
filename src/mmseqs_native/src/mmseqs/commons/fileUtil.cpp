#include <mmseqs/commons/fileUtil.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/util.h>

#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstring>
#include <fstream>

#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <sys/types.h>

bool FileUtil::fileExists(mmseqs_output* out. const char *fileName) {
  struct stat st;
  return stat(fileName, &st) == 0;
}

bool FileUtil::fileExistsAndIsNotEmpty(mmseqs_output* out, const char *fileName) {
  if (!fileExists(fileName)) {
    return false;
  }
  std::ifstream file(fileName);
  if (!file) {
    return false;
  }
  return file.peek() != std::ifstream::traits_type::eof();
}

bool FileUtil::directoryExists(mmseqs_output* out, const char *directoryName) {
  struct stat st;
  return stat(directoryName, &st) == 0 && S_ISDIR(st.st_mode);
}

bool FileUtil::makeDir(mmseqs_output* out, const char *directoryName, const int mode) {
  return mkdir(directoryName, mode) == 0;
}

void *FileUtil::mmapFile(mmseqs_output* out, FILE *file, size_t *dataSize) {
  struct stat sb;
  if (fstat(fileno(file), &sb) < 0) {
    int errsv = errno;
    out->failure("Failed to fstat. Error {}", errsv);
  }
  *dataSize = sb.st_size;
  int mode = PROT_READ;
  int fd = fileno(file);

  void *ret = mmap(NULL, *dataSize, mode, MAP_PRIVATE, fd, 0);
  if (ret == MAP_FAILED) {
    int errsv = errno;
    out->failure("Failed to mmap memory dataSize={}. Error {}.", *dataSize, errsv);
  }
  return ret;
}

void FileUtil::munmapData(mmseqs_output* out, void *ptr, size_t dataSize) {
  if (munmap(ptr, dataSize) < 0) {
    out->failure("Failed to munmap memory");
  }
}

FILE *FileUtil::openFileOrDie(mmseqs_output* out, const char *fileName, const char *mode,
                              bool shouldExist) {
  bool exists = FileUtil::fileExists(fileName);
  if (exists && !shouldExist) {
    errno = EEXIST;
    perror(fileName);
    EXIT(EXIT_FAILURE);
  }
  if (!exists && shouldExist) {
    errno = ENOENT;
    perror(fileName);
    EXIT(EXIT_FAILURE);
  }

  FILE *file;
  file = fopen(fileName, mode);
  if (file == NULL) {
    perror(fileName);
    EXIT(EXIT_FAILURE);
  }
  return file;
}
size_t FileUtil::countLines(mmseqs_output* out, const char *name) {
  FILE *fp = FileUtil::openFileOrDie(name, "r", true);
  size_t cnt = 0;
  while (!feof(fp)) {
    char ch = fgetc(fp);
    cnt += (ch == '\n') ? 1 : 0;
  }
  if (fclose(fp) != 0) {
    out->failure("Cannot close file {}", name);
  }
  return cnt;
}

void FileUtil::deleteTempFiles(mmseqs_output* out, const std::list<std::string> &tmpFiles) {
  for (std::list<std::string>::const_iterator it = tmpFiles.begin();
       it != tmpFiles.end(); it++) {
    out->debug("Deleting {}", *it);
    std::string file = *it;
    FileUtil::remove(file.c_str());
  }
}

void FileUtil::writeFile(mmseqs_output* out, const std::string &pathToFile,
                         const unsigned char *data, size_t len) {
  int fd = open(pathToFile.c_str(), O_WRONLY | O_CREAT | O_TRUNC,
                S_IRUSR | S_IWUSR | S_IXUSR);
  if (fd == -1) {
    out->failure("Could not write file {}", pathToFile);
  }

  ssize_t res = write(fd, data, len);
  if (res == -1) {
    out->failure("Error writing file {}", pathToFile);
  }

  if (close(fd) != 0) {
    out->failure("Error closing file {}", pathToFile);
  }
}

std::string FileUtil::dirName(mmseqs_output* out, const std::string &file) {
  size_t pos = file.find_last_of("\\/");
  return (std::string::npos == pos) ? "." : file.substr(0, pos);
}

std::string FileUtil::baseName(mmseqs_output* out, const std::string &file) {
  size_t pos = file.find_last_of("\\/");
  return (std::string::npos == pos) ? file
                                    : file.substr(pos + 1, file.length());
}

size_t FileUtil::getFreeSpace(mmseqs_output* out, const char *path) {
  struct statvfs stat;
  if (statvfs(path, &stat) != 0) {
    // error happens, just quits here
    return SIZE_MAX;
  }

  // the available size is f_bsize * f_bavail
  return stat.f_bfree * stat.f_frsize;
}

std::string FileUtil::getRealPathFromSymLink(mmseqs_output* out, const std::string path) {
  char *p = realpath(path.c_str(), NULL);
  if (p == NULL) {
    out->failure("Could not get path of {}", path);
  }

  std::string name(p);
  free(p);
  return name;
}

std::string FileUtil::getHashFromSymLink(mmseqs_output* out, const std::string path) {
  char *p = realpath(path.c_str(), NULL);
  if (p == NULL) {
    out->failure("Could not get path of {}", path);
  }

  std::string base = baseName(p);
  free(p);

  return base;
}

void FileUtil::symlinkAlias(mmseqs_output* out, const std::string &file, const std::string &alias) {
  char *p = realpath(file.c_str(), NULL);
  if (p == NULL) {
    out->failure("Could not get path of {}", file);
  }

  std::string path = dirName(p);
  std::string base = baseName(p);
  free(p);

  DIR *dir = opendir(path.c_str());
  if (dir == NULL) {
    out->failure("Error opening directory {}", path);
  }

  std::string pathToAlias = (path + "/" + alias);
  if (symlinkExists(pathToAlias) == true) {
    FileUtil::remove(pathToAlias.c_str());
  }
  // symlinkat is not available in Conda macOS
  // Conda uses the macOS 10.9 SDK, and symlinkat was introduced in 10.10
  // We emulate symlinkat by manipulating the CWD instead
  std::string oldWd = FileUtil::getCurrentWorkingDirectory();
  if (chdir(path.c_str()) != 0) {
    out->failure("Could not change working directory to {}", path);
  }
  if (symlink(base.c_str(), alias.c_str()) != 0) {
    out->failure("Could not create symlink of {}", file);
  }
  if (chdir(oldWd.c_str()) != 0) {
    out->failure("Could not change working directory to {}", oldWd);
  }
  if (closedir(dir) != 0) {
    out->failure("Error closing directory {}", path);
  }
}

std::string FileUtil::getCurrentWorkingDirectory(mmseqs_output* out) {
  // CWD can be larger than PATH_MAX and allocating enough memory is somewhat
  // tricky
  char *wd = NULL;
  size_t bufferSize = PATH_MAX;
  do {
    if (wd != NULL) {
      free(wd);
      bufferSize *= 2;
    }
    wd = getcwd(NULL, bufferSize);
    if (wd == NULL && errno != ERANGE && errno != 0) {
      out->failure("Could not get current working directory");
    }
  } while (wd == NULL && errno == ERANGE);
  std::string cwd(wd);
  free(wd);
  return cwd;
}

void FileUtil::symlinkAbs(mmseqs_output* out, const std::string &target, const std::string &link) {
  if (FileUtil::fileExists(link.c_str())) {
    FileUtil::remove(link.c_str());
  }
  char *t = realpath(target.c_str(), NULL);
  if (t == NULL) {
    out->failure("Could not get realpath of {}", target);
  }

  std::string realLink;
  char *l = realpath(link.c_str(), NULL);
  if (l == NULL) {
    std::string path = dirName(link);
    std::string base = baseName(link);
    l = realpath(path.c_str(), NULL);
    if (l == NULL) {
      out->failure("Could not get realpath of {}", link);
    } else {
      realLink = (std::string(l) + "/" + base);
    }
  } else {
    realLink = l;
    if (symlinkExists(realLink) == true) {
      FileUtil::remove(realLink.c_str());
    }
  }

  if (symlink(t, realLink.c_str()) != 0) {
    out->failure("Could not create symlink of {}", target);
  }

  free(t);
  free(l);
}

size_t FileUtil::getFileSize(mmseqs_output* out, const std::string &fileName) {
  struct stat stat_buf;
  int rc = stat(fileName.c_str(), &stat_buf);
  return rc == 0 ? stat_buf.st_size : -1;
}

bool FileUtil::symlinkExists(mmseqs_output* out, const std::string &path) {
  struct stat buf;
  int result = lstat(path.c_str(), &buf);
  return (result == 0);
}

void FileUtil::copyFile(mmseqs_output* out, const char *src, const char *dst) {
  // https://stackoverflow.com/questions/10195343/copy-a-file-in-a-sane-safe-and-efficient-way
  char buf[BUFSIZ];
  size_t size;

  int source = open(src, O_RDONLY, 0);
  if (source == -1) {
    out->failure("Could not open file {}", src);
  }
  int dest = open(dst, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (dest == -1) {
    out->failure("Could not open file {}", dst);
  }
  while ((size = read(source, buf, BUFSIZ)) > 0) {
    size_t res = write(dest, buf, size);
    if (res != size) {
      out->failure("Error writing file {}", dst);
    }
  }
  close(source);
  close(dest);
}

FILE *FileUtil::openAndDelete(mmseqs_output* out, const char *fileName, const char *mode) {
  if (FileUtil::fileExists(fileName) == true) {
    if (FileUtil::directoryExists(fileName)) {
      out->failure("Can not open {} for writing. It is a directory.", fileName);
    } else {
      FileUtil::remove(fileName);
    }
  }
  FILE *file = fopen(fileName, mode);
  if (file == NULL) {
    out->failure("Can not open {} for writing", fileName);
  }
  return file;
}

std::vector<std::string> FileUtil::findDatafiles(mmseqs_output* out, const char *datafiles) {
  std::string baseName = std::string(datafiles);
  std::string checkName = baseName + ".0";
  std::vector<std::string> filenames;
  size_t cnt = 0;
  while (FileUtil::fileExists(checkName.c_str()) == true) {
    filenames.push_back(checkName);
    cnt++;
    checkName = baseName + "." + SSTR(cnt);
  }
  if (cnt == 0) {
    if (FileUtil::fileExists(baseName.c_str())) {
      filenames.push_back(baseName);
    }
  }
  return filenames;
}

void FileUtil::remove(mmseqs_output* out, const char *file) {
  if (std::remove(file) != 0) {
    out->failure("Could not delete {}", file);
  }
}

void FileUtil::move(mmseqs_output* out, const char *src, const char *dst) {
  struct stat srcFileInfo;
  FILE *srcFile = FileUtil::openFileOrDie(src, "rw", true);
  if (fstat(fileno(srcFile), &srcFileInfo) < 0) {
    int errsv = errno;
    out->failure("Failed to fstat File={}. Error {}", src, errsv);
  }
  struct stat srcDirInfo;
  std::string dirName = FileUtil::dirName(dst);
  FILE *dstDir = FileUtil::openFileOrDie(dirName.c_str(), "r", true);
  if (fstat(fileno(dstDir), &srcDirInfo) < 0) {
    int errsv = errno;
    out->failure("Failed to fstat File={}. Error {}", dirName, errsv);
  }
  bool sameFileSystem = (srcDirInfo.st_dev == srcFileInfo.st_dev);
  if (fclose(srcFile) != 0) {
    out->failure("Cannot close file {}", src);
  }
  if (fclose(dstDir) != 0) {
    out->failure("Cannot close directory {}", dirName);
  }
  if (sameFileSystem) {
    if (std::rename(src, dst) != 0) {
      out->failure("Cannot copy file {} to {}", src, dst);
    }
  } else {
    FileUtil::copyFile(src, dst);
    FileUtil::remove(src);
  }
}

int FileUtil::parseDbType(mmseqs_output* out, const char *name) {
  std::string dbTypeFile = std::string(name) + ".dbtype";
  if (FileUtil::fileExists(dbTypeFile.c_str()) == false) {
    return Parameters::DBTYPE_GENERIC_DB;
  }

  size_t fileSize = FileUtil::getFileSize(dbTypeFile);
  if (fileSize != sizeof(int)) {
    out->failure("File size of {} seems to be wrong. It should have {} bytes but it has {} bytes.", dbTypeFile, sizeof(int), fileSize);
  }
  FILE *file = fopen(dbTypeFile.c_str(), "r");
  if (file == NULL) {
    out->failure("Could not open data file {}", dbTypeFile);
  }
  int dbtype;
  size_t result = fread(&dbtype, 1, fileSize, file);
  if (result != fileSize) {
    out->failure("Could not read {}", dbTypeFile);
  }
  if (fclose(file) != 0) {
    out->failure("Cannot close file {}", dbTypeFile);
  }
  return dbtype;
}

std::string FileUtil::createTemporaryDirectory(
    mmseqs_output* out, std::string baseTmpPath, const std::string &tmpPath,
    const std::string &subDirectory) {
  std::string basePath = baseTmpPath + tmpPath;
  std::string tmpDir(basePath);
  if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
    out->info("Temporary path {} does not exist or is not a directory. It will be created.", tmpDir);
    if (FileUtil::makeDir(tmpDir.c_str()) == false) {
      out->failure("Cannot create temporary folder {}", tmpDir)
    } else {
      out->info("Created temporary directory {}", tmpDir);
    }
  }
  tmpDir += "/" + subDirectory;
  if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
    if (FileUtil::makeDir(tmpDir.c_str()) == false) {
      out->failure("Cannot create temporary subfolder {}", tmpDir);
    }
  }
  FileUtil::symlinkAlias(tmpDir, "latest");
  return tmpDir;
}

void FileUtil::fixRlimitNoFile(mmseqs_output* out) {
  static bool increasedRlimitNoFile(false);
  if (increasedRlimitNoFile == false) {
    increasedRlimitNoFile = true;
    struct rlimit limit;
    if (getrlimit(RLIMIT_NOFILE, &limit) != 0) {
      out->warn("Could not increase maximum number of open files (getrlimit {}). Use ulimit manually.", errno);
      return;
    }
    limit.rlim_cur =
        std::min(std::max((rlim_t)8192, limit.rlim_cur), limit.rlim_max);
    limit.rlim_max = std::min(RLIM_INFINITY, limit.rlim_max);
    if (setrlimit(RLIMIT_NOFILE, &limit) != 0) {
      out->warn("Could not increase maximum number of open files (setrlimit {}). Use ulimit manually.", errno);
    }
  }
}
