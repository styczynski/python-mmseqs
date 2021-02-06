#ifndef DEBUG_H
#define DEBUG_H


#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>

#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstddef>
#include <iostream>
#include <mmseqs/commons/mathUtil.h>
#include <mmseqs/commons/timer.h>

struct FatalException : public std::exception {
   std::string message_body;
   FatalException(std::string message) : message_body(message) {}
   ~FatalException() throw () {}
   const char* what() const throw() { return message_body.c_str(); }
};

class Log {
 public:
  static const int NOTHING = 0;
  static const int ERROR = 1;
  static const int WARNING = 2;
  static const int INFO = 3;

  static int debugLevel;
  std::shared_ptr<spdlog::logger> logger_instance;

  explicit Log(std::string file_log_output) {

     auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
     console_sink->set_level(spdlog::level::trace);
     console_sink->set_pattern("[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v");

     if (file_log_output.size() > 0) {
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("mmseqs_log.log", true);
        file_sink->set_level(spdlog::level::trace);
        spdlog::sinks_init_list sinks = {console_sink, file_sink};
        logger_instance = std::make_shared<spdlog::logger>("mmseqs", sinks);
     } else {
        logger_instance = std::make_shared<spdlog::logger>("mmseqs", console_sink);
     }

     level = Log::INFO;
     setLogLevel(level);
  };

  ~Log() {}

  template<typename FormatString, typename... Args>
  void error(const FormatString &fmt, Args&&...args) {
      logger_instance->error(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void info(const FormatString &fmt, Args&&...args) {
      logger_instance->info(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void warn(const FormatString &fmt, Args&&...args) {
      logger_instance->warn(fmt, std::forward<Args>(args)...);
  }

  template<typename FormatString, typename... Args>
  void debug(const FormatString &fmt, Args&&...args) {
      logger_instance->debug(fmt, std::forward<Args>(args)...);
  }

  // log to str and return it
  template<typename FormatString, typename... Args>
  static std::string log_to_str(const FormatString &msg, Args&&... args) {
    std::ostringstream oss;
    auto oss_sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(oss);
    spdlog::logger oss_logger("pattern_tester", oss_sink);
    oss_logger.set_level(spdlog::level::info);
    oss_logger.info(msg, std::forward<Args>(args)...);
    return oss.str();
  }

  template<typename FormatString, typename... Args>
  void failure(const FormatString &fmt, Args&&...args) {
      logger_instance->error(fmt, std::forward<Args>(args)...);
      // TODO: Implement MPI cleanup
//      int __status = (exitCode);                \
//    if (MMseqsMPI::active && __status == 0) { \
//      MPI_Finalize();                         \
//      MMseqsMPI::active = false;              \
//    }                                         \
//    std::cerr.flush();                        \
//    std::cout.flush();                        \
//    exit(__status);
      throw FatalException(Log::log_to_str(fmt, std::forward<Args>(args)...));
  }

  void setLogLevel(int i) {
    debugLevel = i;
    if (i == Log::NOTHING) {
        logger_instance->set_level(spdlog::level::critical);
    } else if (i == Log::ERROR) {
        logger_instance->set_level(spdlog::level::err);
    } else if (i == Log::WARNING) {
        logger_instance->set_level(spdlog::level::warn);
    } else if (i == Log::INFO) {
        logger_instance->set_level(spdlog::level::debug);
    } else {
        // Level not supported
        logger_instance->set_level(spdlog::level::critical);
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
  int level;
  std::string buffer;
  bool interactive;
};

#endif
