//
// Created by jeremy on 8/11/23.
//

#include "util/Logger.h"

#include <algorithm>
#include <cstdio>
#include "unistd.h"

namespace cafmaker
{

  // -----------------------------------------------------------------------
  Logger::Logger(std::string preamble, THRESHOLD thresh, std::ostream& stream)
    : fPreamble(std::move(preamble)), fThresh(thresh), fStream(stream)
  {
    // will only use colors if stream is a terminal.
    // see https://stackoverflow.com/questions/18081392/discrimination-between-file-and-console-streams
    fIsTerm = ( (&stream == &std::cout && isatty( fileno(stdout) ))
                || (&stream == &std::cerr && isatty( fileno(stderr)  ))
                || (&stream == &std::clog && isatty( fileno(stderr)  )) );
  }

  // -----------------------------------------------------------------------
  // the methods for each THRESHOLD are duplicates.
  // use a macro to eliminate some of the copy-paste
  #define CAFMAKER_INTERNAL_LOGGER_METHOD(threshold, COLOR) \
  const Logger & Logger::threshold() const \
  {                                                         \
      if ( !(fIsMuted = THRESHOLD::threshold < fThresh) )          \
      {                                             \
        const std::string & preamble = !fPendingPreamble.empty() ? fPendingPreamble : fPreamble; \
        (*this) << (fIsTerm ? "\033[94m" : "") << preamble << (!preamble.empty() ? " " : "") << (fIsTerm ? "\033[" #COLOR "m" : "") << #threshold << (fIsTerm ? "\033[00m" : "") << ": ";             \
      }                                             \
      fPendingPreamble.clear();                     \
      return *this;                                 \
  }

  CAFMAKER_INTERNAL_LOGGER_METHOD(VERBOSE, 95)   // magenta
  CAFMAKER_INTERNAL_LOGGER_METHOD(DEBUG, 95)     // magenta
  CAFMAKER_INTERNAL_LOGGER_METHOD(INFO, 32)      // magenta
  CAFMAKER_INTERNAL_LOGGER_METHOD(WARNING, 33)   // yellow
  CAFMAKER_INTERNAL_LOGGER_METHOD(ERROR, 33)     // yellow
  CAFMAKER_INTERNAL_LOGGER_METHOD(FATAL, 91)     // red


  // -----------------------------------------------------------------------
  Logger::THRESHOLD Logger::parseStringThresh(std::string threshStr)
  {
    std::for_each(threshStr.begin(), threshStr.end(), [](char &ch){ch = ::toupper(ch); });

    if (threshStr == "VERBOSE")
      return THRESHOLD::VERBOSE;
    else if (threshStr == "DEBUG")
      return THRESHOLD::DEBUG;
    else if (threshStr == "INFO")
      return THRESHOLD::INFO;
    else if (threshStr == "WARNING")
      return THRESHOLD::WARNING;
    else if (threshStr == "ERROR")
      return THRESHOLD::ERROR;
    else if (threshStr == "FATAL")
      return THRESHOLD::FATAL;

    throw std::runtime_error("Unrecognized log threshold: " + threshStr);
  }

  // -----------------------------------------------------------------------
  Logger & LOG_S()
  {
    static Logger logger = Logger("GLOBAL");
    return logger;
  }

  // -----------------------------------------------------------------------
  const Logger & LOG_S(const std::string& preamble)
  {
    return LOG_S() << Logger::Preamble(preamble);
  }

}
