#include <biosnake/commons/kSeqWrapper.h>
#include <unistd.h>
#include <biosnake/output.h>
#include <biosnake/commons/fileUtil.h>
#include <biosnake/commons/util.h>
#include "kseq.h"

namespace KSEQFILE {
KSEQ_INIT(int, read)
}

KSeqFile::KSeqFile(biosnake_output* output, const char* fileName): out(output) {
  file = FileUtil::openFileOrDie(out, fileName, "r", true);
  seq = (void*)KSEQFILE::kseq_init(fileno(file));
  type = KSEQ_FILE;
}

bool KSeqFile::ReadEntry() {
  KSEQFILE::kseq_t* s = (KSEQFILE::kseq_t*)seq;
  int result = KSEQFILE::kseq_read(s);
  if (result < 0) return false;
  entry.headerOffset = s->headerOffset;
  entry.sequenceOffset = s->sequenceOffset;
  entry.multiline = s->multiline;
  entry.name = s->name;
  entry.comment = s->comment;
  entry.sequence = s->seq;
  entry.qual = s->qual;

  return true;
}

KSeqFile::~KSeqFile() {
  kseq_destroy((KSEQFILE::kseq_t*)seq);
  if (fclose(file) != 0) {
    out->failure("Cannot close KSeq input file");
  }
}

namespace KSEQSTREAM {
KSEQ_INIT(int, read)
}

KSeqStream::KSeqStream(biosnake_output* out) {
  seq = (void*)KSEQSTREAM::kseq_init(STDIN_FILENO);
  type = KSEQ_STREAM;
}

bool KSeqStream::ReadEntry() {
  KSEQSTREAM::kseq_t* s = (KSEQSTREAM::kseq_t*)seq;
  int result = KSEQSTREAM::kseq_read(s);
  if (result < 0) return false;

  entry.name = s->name;
  entry.comment = s->comment;
  entry.sequence = s->seq;
  entry.qual = s->qual;

  return true;
}

KSeqStream::~KSeqStream() { kseq_destroy((KSEQSTREAM::kseq_t*)seq); }

#ifdef HAVE_ZLIB
namespace KSEQGZIP {
KSEQ_INIT(gzFile, gzread)
}

KSeqGzip::KSeqGzip(biosnake_output* out, const char* fileName) {
  if (FileUtil::fileExists(out, fileName) == false) {
    errno = ENOENT;
    perror(fileName);
    out->failure("KSeqGzip: File cannot be loaded: {}", fileName);
  }

  file = gzopen(fileName, "r");
  if (file == NULL) {
    perror(fileName);
    out->failure("KSeqGzip: File cannot be loaded: {}", fileName);
  }

  seq = (void*)KSEQGZIP::kseq_init(file);
  type = KSEQ_GZIP;
}

bool KSeqGzip::ReadEntry() {
  KSEQGZIP::kseq_t* s = (KSEQGZIP::kseq_t*)seq;
  int result = KSEQGZIP::kseq_read(s);
  if (result < 0) return false;

  entry.name = s->name;
  entry.comment = s->comment;
  entry.sequence = s->seq;
  entry.qual = s->qual;
  entry.headerOffset = 0;
  entry.sequenceOffset = 0;
  entry.multiline = s->multiline;

  return true;
}

KSeqGzip::~KSeqGzip() {
  kseq_destroy((KSEQGZIP::kseq_t*)seq);
  gzclose(file);
}
#endif

#ifdef HAVE_BZLIB
namespace KSEQBZIP {
KSEQ_INIT(BZFILE*, BZ2_bzread)
}

KSeqBzip::KSeqBzip(biosnake_output* out, const char* fileName) {
  if (FileUtil::fileExists(out, fileName) == false) {
    errno = ENOENT;
    perror(fileName);
    out->failure("KSeqBzip: Cannot open file {}", fileName);
  }
  FILE* fp = FileUtil::openFileOrDie(out, fileName, "r+b", true);
  int bzError;
  file = BZ2_bzReadOpen(&bzError, fp, 0, 0, NULL, 0);
  if (bzError != 0) {
    perror(fileName);
    out->failure("KSeqBzip: Cannot open file {}", fileName);
  }
  seq = (void*)KSEQBZIP::kseq_init(file);
  type = KSEQ_BZIP;
}

bool KSeqBzip::ReadEntry() {
  KSEQBZIP::kseq_t* s = (KSEQBZIP::kseq_t*)seq;
  int result = KSEQBZIP::kseq_read(s);
  if (result < 0) return false;

  entry.name = s->name;
  entry.comment = s->comment;
  entry.sequence = s->seq;
  entry.qual = s->qual;
  entry.headerOffset = 0;
  entry.sequenceOffset = 0;
  entry.multiline = s->multiline;

  return true;
}

KSeqBzip::~KSeqBzip() {
  kseq_destroy((KSEQBZIP::kseq_t*)seq);
  int bzError;
  BZ2_bzReadClose(&bzError, file);
}
#endif

KSeqWrapper* KSeqFactory(biosnake_output* out, const char* file) {
  KSeqWrapper* kseq = NULL;
  if (strcmp(file, "stdin") == 0) {
    kseq = new KSeqStream(out);
    return kseq;
  }

  if (Util::endsWith(".gz", file) == false &&
      Util::endsWith(".bz2", file) == false) {
    kseq = new KSeqFile(out, file);
    return kseq;
  }
#ifdef HAVE_ZLIB
  else if (Util::endsWith(".gz", file) == true) {
    kseq = new KSeqGzip(out, file);
    return kseq;
  }
#else
  else if (Util::endsWith(".gz", file) == true) {
    out->failure("Biosnake was not compiled with zlib support. Can not read compressed input");
  }
#endif

#ifdef HAVE_BZLIB
  else if (Util::endsWith(".bz2", file) == true) {
    kseq = new KSeqBzip(out, file);
    return kseq;
  }
#else
  else if (Util::endsWith(".bz2", file) == true) {
    out->failure("Biosnake was not compiled with bz2lib support. Can not read compressed input");
  }
#endif

  return kseq;
}

namespace KSEQBUFFER {
KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)
}

KSeqBuffer::KSeqBuffer(biosnake_output* out, const char* buffer, size_t length) {
  d.buffer = (char*)buffer;
  d.length = length;
  d.position = 0;
  seq = (void*)KSEQBUFFER::kseq_init(&d);
  type = KSEQ_BUFFER;
}

bool KSeqBuffer::ReadEntry() {
  KSEQBUFFER::kseq_t* s = (KSEQBUFFER::kseq_t*)seq;
  int result = KSEQBUFFER::kseq_read(s);
  if (result < 0) return false;
  entry.headerOffset = s->headerOffset;
  entry.sequenceOffset = s->sequenceOffset;
  entry.multiline = s->multiline;
  entry.name = s->name;
  entry.comment = s->comment;
  entry.sequence = s->seq;
  entry.qual = s->qual;

  return true;
}

KSeqBuffer::~KSeqBuffer() { kseq_destroy((KSEQBUFFER::kseq_t*)seq); }
