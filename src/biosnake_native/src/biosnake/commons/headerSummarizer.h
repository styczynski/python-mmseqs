#ifndef BIOSNAKE_HEADERSUMMARIZER_H
#define BIOSNAKE_HEADERSUMMARIZER_H

#include <biosnake/output.h>
#include <string>
#include <vector>

class HeaderSummarizer {
 public:
  virtual std::string summarize(biosnake_output* out, const std::vector<std::string>& headers) = 0;
  virtual ~HeaderSummarizer(){};
};

class UniprotHeaderSummarizer : public HeaderSummarizer {
 public:
  std::string summarize(biosnake_output* out, const std::vector<std::string>& headers);
  ~UniprotHeaderSummarizer(){};
};

class MetaclustHeaderSummarizer : public HeaderSummarizer {
 public:
  std::string summarize(biosnake_output* out, const std::vector<std::string>& headers);
  ~MetaclustHeaderSummarizer(){};
};

#endif  // BIOSNAKE_HEADERSUMMARIZER_H
