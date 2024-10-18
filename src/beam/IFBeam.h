/// \file IFBeam.h
///
/// IFBeam class to query beam information
///
/// \author  S. Kumaran <s.kumaran@uci.edu>
/// \date    Oct. 2024

#ifndef IFBEAM_H
#define IFBEAM_H

#include <Params.h>
#include "reco/IRecoBranchFiller.h"
#include <string>
#include <map>
#include <vector>


namespace cafmaker
{
  class IFBeam {
  public:
      using BeamSpills = std::map<double, double>;
      using TriggerGroup = std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>;
  
      IFBeam(const std::vector<TriggerGroup>& groupedTriggers, bool is_data);   
  
      double getPOT(const cafmaker::Params& par, const TriggerGroup & groupedTrigger, int ii);
   
  
  private:
      const std::string potDevice = "E:TRTGTD";
      BeamSpills beamSpills;
  
      void loadBeamSpills(const std::vector<TriggerGroup>& groupedTriggers);
      std::string createUrl(const std::string& min_time_iso, const std::string& max_time_iso);
      double unitToFactor(const std::string& unit);
  };
   
}
#endif // IFBEAM_H
