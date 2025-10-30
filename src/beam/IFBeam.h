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
      using BeamInfo = std::map<double, std::vector<double>>;
      using TriggerGroup = std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>;
  
      IFBeam(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers, bool is_data);   

      double getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const BeamSpills& data);
      std::vector<double> getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const BeamInfo& data);
      double getPOT(const cafmaker::Params& par, const TriggerGroup & groupedTrigger, int ii);
      double getHornI(const cafmaker::Params& par, const TriggerGroup & groupedTrigger, int ii);
      std::vector<double> getHorizontalPosTGT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii);
      std::vector<double> getHorizontalIntTGT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii);
      std::vector<double> getHorizontalPos121(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii);
      std::vector<double> getVerticalPosTGT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii);
      std::vector<double> getVerticalIntTGT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii);
      std::vector<double> getVerticalPos121(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii);
      std::vector<double> getMultiWireInfo(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii);
   
  
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
      //__________Beam Position_____________
      const std::string horizontalPosTGTDevice = "E:HPTGT[]";
      const std::string horizontalIntTGTDevice = "E:HITGT[]";
      const std::string horizontalPos121Device = "E:HP121[]";
      const std::string verticalPosTGTDevice = "E:VPTGT[]";
      const std::string verticalIntTGTDevice = "E:VITGT[]";
      const std::string verticalPos121Device = "E:VP121[]";
      //__________Beam Width________________
      const std::string multiwireDevice = "E:MTGTDS[]";

      BeamSpills beamSpills;
      BeamSpills hornCurrentA;
      BeamSpills hornCurrentB;
      BeamSpills hornCurrentC;
      BeamSpills hornCurrentD;
      BeamSpills hornI; // linear combination
      BeamInfo horizontalPosTGT;
      BeamInfo horizontalIntTGT;
      BeamInfo horizontalPos121;
      BeamInfo verticalPosTGT;
      BeamInfo verticalIntTGT;
      BeamInfo verticalPos121;
      BeamInfo multiwireInfo;
  
      void retrieveInfoFromDataBase(const std::string url, BeamSpills& data);
      void retrieveInfoFromDataBase(const std::string url, BeamInfo& data);
      void loadBeamSpills(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers);
      std::string createUrl(const std::string potDevice, const std::string& min_time_iso, const std::string& max_time_iso);
      double unitToFactor(const std::string& unit);
  };
   
}
#endif // IFBEAM_H
