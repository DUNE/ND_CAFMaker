/// \file Logger.h
///
/// Very rudimentary logging facility.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    Aug. 2023


#ifndef ND_CAFMAKER_LOGGER_H
#define ND_CAFMAKER_LOGGER_H

#include <iostream>

namespace cafmaker
{
  /// Rudimentary logger facility.  Currently only supports stdout...
  /// (would need some re-engineering to directly write to file
  ///  so that different instances didn't clobber one another)
  class Logger
  {
    public:
      enum class THRESHOLD { VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL };
      static THRESHOLD parseStringThresh(std::string threshStr);

      /// Insert a Preamble object into a Logger() stream
      /// to change the preamble used in the next invocation
      struct Preamble
      {
        explicit Preamble(std::string t)
          : text(std::move(t))
        {}

        std::string text;
      };

      explicit Logger(std::string preamble,
                      THRESHOLD thresh = THRESHOLD::INFO,
                      std::ostream& stream = std::cout);

      void SetPreamble(std::string preamble)   { fPreamble = std::move(preamble); }

      THRESHOLD GetThreshold() const           { return fThresh; }
      void      SetThreshold(THRESHOLD thresh) { fThresh = thresh; }

      /// Force a write to the output.
      template <typename T>
      const Logger & operator<<(const T & obj) const
      {
        // if passed a Preamble, we're just adjusting the preamble rather than writing to the stream
        if constexpr (std::is_same_v<T, Logger::Preamble>)
          fPendingPreamble = obj.text;
        else
        {
          if (!fIsMuted)
            fStream << obj;
        }
        return *this;
      }

      const Logger & VERBOSE() const;
      const Logger & DEBUG() const;
      const Logger & INFO() const;
      const Logger & WARNING() const;
      const Logger & ERROR() const;
      const Logger & FATAL() const;

    private:
      std::string fPreamble;        ///< write this in front of every message

      mutable std::string fPendingPreamble; ///< temporary preamble used only for the next message
      mutable bool fIsMuted = false;  ///<  is this stream currently muted?

      std::ostream & fStream;  ///<  where output will be written

      THRESHOLD fThresh;      ///<  current log threshold
      bool fIsTerm;           ///<  is this output stream a terminal?
  };

  /// Retrieve the global logger for stream use (preamble setting)
  const Logger & LOG_S(const std::string& preamble);

  /// Retrieve the global logger object for general use (including setting the log threshold)
  Logger & LOG_S();

} // cafmaker

#endif //ND_CAFMAKER_LOGGER_H
