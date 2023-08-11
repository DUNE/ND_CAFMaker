/// \file Logger.h
///
/// Base class for any algorithm that wants to inherit Logger functionality
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    Aug. 2023

#ifndef ND_CAFMAKER_LOGGABLE_H
#define ND_CAFMAKER_LOGGABLE_H

#include <string>

#include "util/Logger.h"

namespace cafmaker
{

  /// Base class for other classes that want their own Logger instances
  class Loggable
  {
    public:
      Loggable(std::string name="", cafmaker::Logger::THRESHOLD defaultThresh=cafmaker::Logger::THRESHOLD::INFO)
        : LOG(std::move(name), defaultThresh)
      {}

      void SetLogThrehsold(cafmaker::Logger::THRESHOLD thresh) { LOG.SetThreshold(thresh); }

      virtual ~Loggable() = default;

    protected:
      cafmaker::Logger LOG;
  };

} // cafmaker

#endif //ND_CAFMAKER_LOGGABLE_H
