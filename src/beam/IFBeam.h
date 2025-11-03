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
      using TriggerGroup = std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>;
  
      IFBeam(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers, bool is_data);   

      std::vector<double> getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const BeamInfo& data);
      std::vector<double> getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii,const std::string& deviceName);
  
  private:
      // _________POT _______________________
      const std::string potTRTGTDDevice = "E:TRTGTD";
      const std::string potTOR101Device = "E:TOR101";
      const std::string potTR101DDevice = "E:TR101D";
      // _________Horn Current_______________
      // The Horn current is the sum of the normalized and shifted values for the horn current from the four strip-lines : E:NSLINA, E:NSLINB, E:NSLINC, E:NSLIND
      // see: https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/IFDBSpillInfo_module.cc#L685
      // the horn current is calculated as I_horn = (E:NSLINA-(+0.01))/0.9951 + (E:NSLINB-(-0.14))/0.9957 + (E:NSLINC-(-0.05))/0.9965 + (E:NSLIND-(-0.07))/0.9945
      const std::string hornCurrentDeviceA = "E:NSLINA";
      const std::string hornCurrentDeviceB = "E:NSLINB";
      const std::string hornCurrentDeviceC = "E:NSLINC";
      const std::string hornCurrentDeviceD = "E:NSLIND";
      const std::string hornCurrentDevice = "E:NSLIN";
      // ________Horn polarity_______________
      const std::string hornDirDevice = "E:HRNDIR";
      //__________Beam Position_____________
      const std::string horizontalPosTGTDevice = "E:HPTGT[]";
      const std::string horizontalIntTGTDevice = "E:HITGT[]";
      const std::string horizontalPos121Device = "E:HP121[]";
      const std::string verticalPosTGTDevice = "E:VPTGT[]";
      const std::string verticalIntTGTDevice = "E:VITGT[]";
      const std::string verticalPos121Device = "E:VP121[]";
      //__________Beam Width________________
      const std::string multiwireDevice = "E:MTGTDS[]";

      BeamInfo potTRTGTD;
      BeamInfo potTOR101;
      BeamInfo potTR101D;
      BeamInfo hornCurrentA;
      BeamInfo hornCurrentB;
      BeamInfo hornCurrentC;
      BeamInfo hornCurrentD;
      BeamInfo hornI; // linear combination of A,B,C,D
      BeamInfo hornDir; 
      BeamInfo horizontalPosTGT;
      BeamInfo horizontalIntTGT;
      BeamInfo horizontalPos121;
      BeamInfo verticalPosTGT;
      BeamInfo verticalIntTGT;
      BeamInfo verticalPos121;
      BeamInfo multiwireInfo;
  
      void retrieveInfoFromDataBase(const std::string url, BeamInfo& data);
      void loadBeamSpills(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers);
      std::string createUrl(const std::string potDevice, const std::string& min_time_iso, const std::string& max_time_iso);
      double unitToFactor(const std::string& unit);
  };
   
}
#endif // IFBEAM_H
