//
// Created by mad on 3/24/16.
//

#ifndef BIOSNAKE_CONCAT_H
#define BIOSNAKE_CONCAT_H

#include <_simd/simd.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>

#include <biosnake/output.h>
#include <biosnake/commons/util.h>
#define SAFE_READ_ERROR ((size_t)-1)

class Concat {
 public:
  /* Read(write) up to COUNT bytes at BUF from(to) descriptor FD, retrying if
     interrupted.  Return the actual number of bytes read(written), zero for
     EOF, or SAFE_READ_ERROR(SAFE_WRITE_ERROR) upon error.  */
  static size_t safe_write(int fd, void const *buf, size_t count) {
    ssize_t result;
    /* POSIX limits COUNT to SSIZE_MAX, but we limit it further, requiring
       that COUNT <= INT_MAX, to avoid triggering a bug in Tru64 5.1.
       When decreasing COUNT, keep the file pointer block-aligned.
       Note that in any case, read(write) may succeed, yet read(write)
       fewer than COUNT bytes, so the caller must be prepared to handle
       partial results.  */
    if (count > INT_MAX) count = INT_MAX & ~8191;

    do {
      result = write(fd, buf, count);
    } while (result < 0 && ((errno) == EINTR));
    return (size_t)result;
  }

  /*
   * Write all of the supplied buffer out to a file.
   * This does multiple writes as necessary.
   * Returns the amount written, or -1 on an error.
   */
  static ssize_t full_write(int fd, const void *buf, size_t len) {
    ssize_t cc;
    ssize_t total;

    total = 0;

    while (len) {
      cc = safe_write(fd, buf, len);

      if (cc < 0) {
        if (total) {
          /* we already wrote some! */
          /* user can do another write to know the error code */
          return total;
        }
        return cc; /* write() returns -1 on failure. */
      }

      total += cc;
      buf = ((const char *)buf) + cc;
      len -= cc;
    }

    return total;
  }

  /* Read(write) up to COUNT bytes at BUF from(to) descriptor FD, retrying if
     interrupted.  Return the actual number of bytes read(written), zero for
     EOF, or SAFE_READ_ERROR(SAFE_WRITE_ERROR) upon error.  */
  static size_t safe_read(int fd, void const *buf, size_t count) {
    ssize_t result;

    /* POSIX limits COUNT to SSIZE_MAX, but we limit it further, requiring
       that COUNT <= INT_MAX, to avoid triggering a bug in Tru64 5.1.
       When decreasing COUNT, keep the file pointer block-aligned.
       Note that in any case, read(write) may succeed, yet read(write)
       fewer than COUNT bytes, so the caller must be prepared to handle
       partial results.  */
    if (count > INT_MAX) count = INT_MAX & ~8191;

    do {
      result = read(fd, (void *)buf, count);
    } while (result < 0 && ((errno) == EINTR));

    return (size_t)result;
  }

  static void *ptr_align(void const *ptr, size_t alignment) {
    char const *p0 = static_cast<const char *>(ptr);
    char const *p1 = p0 + alignment - 1;
    std::cout << "(size_t) p1 % alignment=" << (size_t)p1 % alignment
              << std::endl;
    return (void *)(p1 - (size_t)p1 % alignment);
  }

  static size_t io_blksize(const struct stat &sb) {
    const size_t IO_BUFSIZE = 64 * 1024;
    return std::max(IO_BUFSIZE, static_cast<size_t>(sb.st_blksize));
  }

  static void concatFiles(biosnake_output* out, const std::vector<FILE *> &files, FILE *outFile) {
    int output_desc = fileno(outFile);
    struct stat stat_buf;
    if (fstat(output_desc, &stat_buf) < 0) {
      out->failure("Error with output file");
    }
    size_t outsize = io_blksize(stat_buf);
    /* Actual number of characters read, and therefore written.  */
    for (size_t fileIdx = 0; fileIdx < files.size(); fileIdx++) {
      int input_desc = fileno(files[fileIdx]);
      if (fstat(input_desc, &stat_buf) < 0) {
        out->failure("Error with input descriptor");
      }

      size_t insize = io_blksize(stat_buf);
      insize = std::max(insize, outsize);

      size_t page_size = getpagesize();
      char *inbuf = (char *)mem_align(page_size, insize);
#if HAVE_POSIX_FADVISE
      if (posix_fadvise(input_desc, 0, 0, POSIX_FADV_SEQUENTIAL) != 0) {
        out->failure("posix_fadvise returned an error");
      }
#endif
      doConcat(out, input_desc, output_desc, inbuf, insize);
      free(inbuf);
    }
  }

  static bool doConcat(biosnake_output* out, int input_desc, int out_desc, const char *buf,
                       size_t bufsize) {
    while (true) {
      /* Read a block of input.  */

      size_t n_read = safe_read(input_desc, buf, bufsize);
      if (n_read == SAFE_READ_ERROR) {
        out->error("Read error nr: {}", errno);
        return false;
      }

      /* End of this file?  */

      if (n_read == 0) return true;

      /* Write this block out.  */

      {
        /* The following is ok, since we know that 0 < n_read.  */
        if (full_write(out_desc, buf, n_read) != (ssize_t)n_read) {
          out->failure("Read error.");
        }
      }
    }
  }
};

#undef SAFE_READ_ERROR

#endif  // BIOSNAKE_CONCAT_H
