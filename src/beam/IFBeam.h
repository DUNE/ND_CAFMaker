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

      double getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const BeamSpills& data);
      double getPOT(const cafmaker::Params& par, const TriggerGroup & groupedTrigger, int ii);
      double getHornI(const cafmaker::Params& par, const TriggerGroup & groupedTrigger, int ii);
   
  
  private:
      // _________POT _______________________
      // pot value from toroids : E:TOR101, TORTGT
      const std::string potDevice = "E:TRTGTD";
      const std::string potDevice2 = "E:TOR101";
      // _________Horn Current_______________
      // The Horn current is the sum of the normalized and shifted values for the horn
      // current from the four strip-lines : E:NSLINA, E:NSLINB, E:NSLINC, E:NSLIND
      // see: https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/IFDBSpillInfo_module.cc#L685
      // the horn current is calculated as I_horn = (E:NSLINA-(+0.01))/0.9951 + (E:NSLINB-(-0.14))/0.9957 + (E:NSLINC-(-0.05))/0.9965 + (E:NSLIND-(-0.07))/0.9945
      const std::string hornCurrentDeviceA = "E:NSLINA";
      const std::string hornCurrentDeviceB = "E:NSLINB";
      const std::string hornCurrentDeviceC = "E:NSLINC";
      const std::string hornCurrentDeviceD = "E:NSLIND";
      // ________Horn polarity_______________
      // horn polarity E:HRNDIR
      BeamSpills beamSpills;
      BeamSpills hornCurrent;
  
      void loadBeamSpills(const std::vector<TriggerGroup>& groupedTriggers);
      std::string createUrl(const std::string potDevice, const std::string& min_time_iso, const std::string& max_time_iso);
      double unitToFactor(const std::string& unit);
  };
   
}
#endif // IFBEAM_H
