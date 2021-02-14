#ifndef BIOSNAKE_KSEQWRAPPER_H
#define BIOSNAKE_KSEQWRAPPER_H

#include <biosnake/output.h>
#include <biosnake/commons/kSeqBufferReader.h>
#include "kseq.h"

#include <string>

class KSeqWrapper {
 public:
  struct KSeqEntry {
    kstring_t name;
    kstring_t sequence;
    kstring_t comment;
    kstring_t qual;
    size_t headerOffset;
    size_t sequenceOffset;
    bool multiline;
  } entry;

  enum kseq_type { KSEQ_FILE, KSEQ_STREAM, KSEQ_GZIP, KSEQ_BZIP, KSEQ_BUFFER };
  kseq_type type;

  virtual bool ReadEntry() = 0;
  virtual ~KSeqWrapper(){};

 protected:
  void* seq;
};

class KSeqFile : public KSeqWrapper {
 public:
  KSeqFile(biosnake_output* out, const char* file);
  bool ReadEntry();
  ~KSeqFile();

 private:
  biosnake_output* out;
  FILE* file;
};

class KSeqStream : public KSeqWrapper {
 public:
  KSeqStream(biosnake_output* out);
  bool ReadEntry();
  ~KSeqStream();
};

#ifdef HAVE_ZLIB
#include <zlib.h>

class KSeqGzip : public KSeqWrapper {
 public:
  KSeqGzip(biosnake_output* out, const char* file);
  bool ReadEntry();
  ~KSeqGzip();

 private:
  gzFile file;
};
#endif

#ifdef HAVE_BZLIB
#include <bzlib.h>

class KSeqBzip : public KSeqWrapper {
 public:
  KSeqBzip(biosnake_output* out, const char* file);
  bool ReadEntry();
  ~KSeqBzip();

 private:
  BZFILE* file;
};
#endif

class KSeqBuffer : public KSeqWrapper {
 public:
  KSeqBuffer(biosnake_output* out, const char* buffer, size_t length);
  bool ReadEntry();
  ~KSeqBuffer();

 private:
  kseq_buffer_t d;
};

KSeqWrapper* KSeqFactory(biosnake_output* out, const char* file);

#endif  // BIOSNAKE_KSEQWRAPPER_H
