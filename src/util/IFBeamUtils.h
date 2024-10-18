/// \file IFBeamUtils.cxx
///
/// Utility functions for IFBeam database
///
/// \author  S. Kumaran <s.kumaran@uci.edu>
/// \date    Oct. 2024

#ifndef ND_CAFMAKER_IFBeamUtils_H
#define ND_CAFMAKER_IFBeamUtils_H

#include "reco/IRecoBranchFiller.h"
#include <string>

namespace cafmaker
{
  namespace util
  {

  std::string toISO8601(double time_sec);
  double getTriggerTime(const Trigger& trigger);
  size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* userp);
  }

}

#endif //ND_CAFMAKER_IFBeamUtils_H
