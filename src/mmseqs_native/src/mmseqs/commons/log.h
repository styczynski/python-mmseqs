#ifndef DEBUG_H
#define DEBUG_H


#include <spdlog/spdlog.h>

#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstddef>
#include <iostream>
#include <mmseqs/commons/mathUtil.h>
#include <mmseqs/commons/timer.h>
#include <mmseqs/commons/util.h>

class Log {
 public:
  static const int NOTHING = 0;
  static const int ERROR = 1;
  static const int WARNING = 2;
  static const int INFO = 3;

  static int debugLevel;
  std::shared_ptr<logger> logger_instance;

  explicit Log() {
     logger_instance = spdlog::stdout_color_mt("console");
     logger_instance->set_pattern("[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v");
     level = Debug::INFO;
     setLogLevel(level);
  };

  ~Log() {

  }

  template<typename FormatString, typename... Args>
  void error(const FormatString &fmt, Args&&...args) {
      spdlog::get("console")->spdlog::error(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void info(const FormatString &fmt, Args&&...args) {
      spdlog::get("console")->spdlog::info(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void debug(const FormatString &fmt, Args&&...args) {
      spdlog::get("console")->spdlog::debug(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void failure(const FormatString &fmt, Args&&...args) {
      spdlog::get("console")->spdlog::error(fmt, std::forward<Args>(args)...);
      EXIT(EXIT_FAILURE);
  }

  static void setLogLevel(int i) {
    debugLevel = i;
    if (i == Log::NOTHING) {
        spdlog::get("console")->set_level(spdlog::level::critical);
    } else if (i == Log::ERROR) {
        spdlog::get("console")->set_level(spdlog::level::err);
    } else if (i == Log::WARNING) {
        spdlog::get("console")->set_level(spdlog::level::warn);
    } else if (i == Log::INFO) {
        spdlog::get("console")->set_level(spdlog::level::debug);
    } else {
        // Level not supported
        spdlog::get("console")->set_level(spdlog::level::critical);
    }
  }

  class Progress {
   private:
    size_t currentPos;
    size_t prevPrintedId;
    size_t totalEntries;
    bool interactive;
    Timer timer;

    const static int BARWIDTH = 65;

    std::string buildItemString(size_t id) {}

   public:
    Progress(size_t totalEntries)
        : currentPos(0), prevPrintedId(0), totalEntries(totalEntries) {}

    Progress() : currentPos(0), prevPrintedId(0), totalEntries(SIZE_MAX) {}

    void reset(size_t totalEntries) {}

    void updateProgress() {}
  };

 private:
  const int level;
  std::string buffer;
  bool interactive;
};

#endif
