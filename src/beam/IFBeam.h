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
      using BeamInfo = std::map<double, std::vector<double>>;
      using DeviceMap = std::map<const std::string, BeamInfo>;
      using TriggerGroup = std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>;
  
      IFBeam(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers, bool is_data);   

      std::vector<double> getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const BeamInfo& data);
      std::vector<double> getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii,const std::string& deviceName);
  
  private:
      DeviceMap deviceMap = {
      //_________POT _______________________
        {"E:TRTGTD", {}},
        {"E:TOR101", {}},
        {"E:TR101D", {}},
      //_________Horn Current_______________
        {"E:NSLINA", {}},
        {"E:NSLINB", {}},
        {"E:NSLINC", {}},
        {"E:NSLIND", {}},
        {"E:NSLIN", {}},
      //________Horn polarity_______________
        {"E:HRNDIR", {}},
      //__________Beam Position_____________
        {"E:HPTGT[]", {}},
        {"E:HITGT[]", {}},
        {"E:HP121[]", {}},
        {"E:VPTGT[]", {}},
        {"E:VITGT[]", {}},
        {"E:VP121[]", {}},
      //__________Beam Width________________
        {"E:MTGTDS[]", {}}
      };

      BeamInfo retrieveInfoFromDataBase(const std::string url);
      void loadBeamSpills(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers);
      std::string createUrl(const std::string potDevice, const std::string& min_time_iso, const std::string& max_time_iso);
      double unitToFactor(const std::string& unit);
  };
   
}
#endif // IFBEAM_H
