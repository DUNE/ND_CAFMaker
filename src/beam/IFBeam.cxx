/// \file IFBeam.cxx
///
/// IFBeam class to query beam information
///
/// \author  S. Kumaran <s.kumaran@uci.edu>
/// \date    Oct. 2024

#include "IFBeam.h"
#include "util/Logger.h"
#include "util/IFBeamUtils.h"
#include <curl/curl.h>
#include <iostream>
#include <regex>
#include <sstream>
#include <limits>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace cafmaker
{

  
  IFBeam::IFBeam(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers, bool is_data) {
      if (is_data)
  	loadBeamSpills(par, groupedTriggers); // Load all beam spills upon instantiation if data
  }
  
  std::string IFBeam::createUrl(const std::string potDevice, const std::string& min_time_iso, const std::string& max_time_iso) {
      std::ostringstream url_stream;
      url_stream << "https://dbdata3vm.fnal.gov:9443/ifbeam/data/data?v=" << potDevice
                 << "&e=" << "e,a9"
                 << "&t0=" << min_time_iso
                 << "&t1=" << max_time_iso
                 << "&f=json";
      return url_stream.str();
  }
  
  double IFBeam::unitToFactor(const std::string& unit) { //because the exponent part in IFBeam is stored as a string under "value" :/
      std::regex re("E(\\d+)");
      std::smatch match;
  
      if (std::regex_search(unit, match, re) && match.size() > 1) {
          int exponent = std::stoi(match.str(1));
          return std::pow(10, exponent);
      } else {
          cafmaker::LOG_S("Fetch beam information").WARNING() << "Unknown unit in beam database: " << unit << "\n";
          return 1.0;
      }
  }

  void IFBeam::retrieveInfoFromDataBase(const std::string url, BeamInfo& data) {

      std::cout << __func__ << " line " << __LINE__ <<": data has " << data.size() << " entries\n";
      double ms_to_s = 1e-3;
 
      CURL* curl;
      CURLcode res;
      std::string readBuffer; 
      std::cout << "Fetching beam information from: " << url << "\n";
      curl = curl_easy_init();
      
      if (curl) {
          curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
          curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, util::WriteCallback);
          curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
          res = curl_easy_perform(curl);
          if (res != CURLE_OK) {
              cafmaker::LOG_S("Fetch beam information").ERROR() << "curl_easy_perform() failed: " << curl_easy_strerror(res) << "\n";
              std::abort();
          }
          curl_easy_cleanup(curl);
          try {
              auto json_data = json::parse(readBuffer);
              // std::cout << "readBuffer " << readBuffer << "\n";
              for (const auto& spill : json_data["rows"]) {
                  double time = spill["clock"].get<double>() * ms_to_s;
                  std::string unit = spill["units"];
                  auto& c_array = spill["c"];

                  std::vector<double> values;
                  for (const auto& c_entry : c_array) 
                  {
                    if (c_entry["v"].is_number())
                    {
                      double this_v = c_entry["v"].get<double>() * unitToFactor(unit);
                      values.push_back(this_v);
                    }
                  }
                  data[time] = values;
                  std::cout << " data at time : " << time << ", has " << values.size() << " measurements \n";
              }
          } catch (const json::exception& e) {
              cafmaker::LOG_S("Fetch beam information").ERROR() << "Failed to parse JSON data: " << e.what() << "\n";
  	    std::abort();
          }
      }
  }

  void IFBeam::retrieveInfoFromDataBase(const std::string url, BeamSpills& data) {

      double ms_to_s = 1e-3;
 
      CURL* curl;
      CURLcode res;
      std::string readBuffer; 
      std::cout << "Fetching beam information from: " << url << "\n";
      curl = curl_easy_init();
      
      if (curl) {
          curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
          curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, util::WriteCallback);
          curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
          res = curl_easy_perform(curl);
          if (res != CURLE_OK) {
              cafmaker::LOG_S("Fetch beam information").ERROR() << "curl_easy_perform() failed: " << curl_easy_strerror(res) << "\n";
              std::abort();
          }
          curl_easy_cleanup(curl);
          try {
              auto json_data = json::parse(readBuffer);
              // std::cout << "readBuffer " << readBuffer << "\n";
              for (const auto& spill : json_data["rows"]) {
                  double time = spill["clock"].get<double>() * ms_to_s;
                  std::string unit = spill["units"];
                  double value = spill["value"].get<double>() * unitToFactor(unit);
                  data[time] = value;
              }
          } catch (const json::exception& e) {
              cafmaker::LOG_S("Fetch beam information").ERROR() << "Failed to parse JSON data: " << e.what() << "\n";
  	    std::abort();
          }
      }
  }
  
  void IFBeam::loadBeamSpills(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers) { //todo: this should return other information as well like horn current, position, etc., querying all devices and storing in a map
      beamSpills.clear();
      hornCurrentA.clear();
      hornCurrentB.clear();
      hornCurrentC.clear();
      hornCurrentD.clear();
      hornI.clear();
      horizontalPosTGT.clear();
      horizontalIntTGT.clear();
      horizontalPos121.clear();
      verticalPosTGT.clear();
      verticalIntTGT.clear();
      verticalPos121.clear();
      multiwireInfo.clear();

      double dt = 5.0; //time window to query before and after first and last spill, respectively
      double min_time = std::numeric_limits<double>::max();
      double max_time = std::numeric_limits<double>::lowest();
      for (const auto& group : groupedTriggers) {
          for (const auto& trig : group) {
              if (!trig.first->IsBeamTrigger(trig.second.triggerType)) continue;
              double trigger_time = util::getTriggerTime(trig.second);

              min_time = std::min(min_time, trigger_time - dt);
              max_time = std::max(max_time, trigger_time + dt);
          }
      }
      if (min_time == std::numeric_limits<double>::max() || max_time == std::numeric_limits<double>::lowest())
      {
        std::cerr<<"No beam trigger loaded, are you sure you want to use IFBeam database? \n";
        return;
      }
      std::string min_time_iso = util::toISO8601(min_time);
      std::string max_time_iso = util::toISO8601(max_time);
      std::string url_pot = createUrl(potDevice, min_time_iso, max_time_iso);
      std::string url_hornI_A = createUrl(hornCurrentDeviceA, min_time_iso, max_time_iso);
      std::string url_hornI_B = createUrl(hornCurrentDeviceB, min_time_iso, max_time_iso);
      std::string url_hornI_C = createUrl(hornCurrentDeviceC, min_time_iso, max_time_iso);
      std::string url_hornI_D = createUrl(hornCurrentDeviceD, min_time_iso, max_time_iso);
      std::string url_horizontalPosTGT = createUrl(horizontalPosTGTDevice, min_time_iso, max_time_iso);
      std::string url_horizontalIntTGT = createUrl(horizontalIntTGTDevice, min_time_iso, max_time_iso);
      std::string url_horizontalPos121 = createUrl(horizontalPos121Device, min_time_iso, max_time_iso);
      std::string url_verticalPosTGT = createUrl(verticalPosTGTDevice, min_time_iso, max_time_iso);
      std::string url_verticalIntTGT = createUrl(verticalIntTGTDevice, min_time_iso, max_time_iso);
      std::string url_verticalPos121 = createUrl(verticalPos121Device, min_time_iso, max_time_iso);
      std::string url_multiwire = createUrl(multiwireDevice, min_time_iso, max_time_iso);
  
      retrieveInfoFromDataBase(url_pot, beamSpills);
      retrieveInfoFromDataBase(url_hornI_A, hornCurrentA);
      retrieveInfoFromDataBase(url_hornI_B, hornCurrentB);
      retrieveInfoFromDataBase(url_hornI_C, hornCurrentC);
      retrieveInfoFromDataBase(url_hornI_D, hornCurrentD);
      retrieveInfoFromDataBase(url_horizontalPosTGT, horizontalPosTGT);
      retrieveInfoFromDataBase(url_horizontalIntTGT, horizontalIntTGT);
      retrieveInfoFromDataBase(url_horizontalPos121, horizontalPos121);
      retrieveInfoFromDataBase(url_verticalPosTGT, verticalPosTGT);
      retrieveInfoFromDataBase(url_verticalIntTGT, verticalIntTGT);
      retrieveInfoFromDataBase(url_verticalPos121, verticalPos121);
      retrieveInfoFromDataBase(url_multiwire, multiwireInfo);

      // see: https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/IFDBSpillInfo_module.cc#L685
      // the horn current is calculated as I_horn = (E:NSLINA-(+0.01))/0.9951 + (E:NSLINB-(-0.14))/0.9957 + (E:NSLINC-(-0.05))/0.9965 + (E:NSLIND-(-0.07))/0.9945
      
      double currentA=0.0, currentB=0.0, currentC=0.0, currentD=0.0;

      for(const auto& pairA : hornCurrentA) {
        auto timeA = pairA.first;
        auto currentA = pairA.second; 
        
        if(currentA == 0.) continue;

        for(const auto& pairB : hornCurrentB) {
          auto timeB = pairB.first;
          if (abs(timeA - timeB) <= par().cafmaker().beamMatchDT()) {
            currentB = pairB.second;
            break;
          }
        }

        if(currentB == 0.) continue;

        for(const auto& pairC : hornCurrentC) {
          auto timeC = pairC.first;
          if (abs(timeA - timeC) <= par().cafmaker().beamMatchDT()) 
          {
            currentC = pairC.second;
            break;
          }
        }

        if(currentC == 0.) continue;

        for(const auto& pairD : hornCurrentD) {
          auto timeD = pairD.first;
          if (abs(timeA - timeD) <= par().cafmaker().beamMatchDT()){
            currentD = pairD.second;
            break;
          }
        }

        if(currentD == 0.) continue;
        
        hornI[timeA] = (currentA-0.01)/0.9951 + (currentB+0.14)/0.9957 + (currentC+0.05)/0.9965 + (currentD+0.07)/0.9945;
      }

  }
  

  double IFBeam::getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const BeamSpills& data) {
      double pot = 0.0;
      if (!(groupedTrigger.front().first->IsBeamTrigger(groupedTrigger.front().second.triggerType))) return 0.0;
      
      auto it = std::find_if(data.begin(), data.end(),
          [par, &groupedTrigger](const auto& spill) {
              return std::all_of(groupedTrigger.cbegin(), groupedTrigger.cend(),
                  [par, &spill](const auto& groupedTrigger) {
                      return std::abs(util::getTriggerTime(groupedTrigger.second) - spill.first) < par().cafmaker().beamMatchDT();
                  });
          });
  

      if (it != data.end()) {
          pot = it->second;
      } 
      else {
          bool any_matched = false;
          std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> matched_triggers;
          std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> unmatched_triggers;
          for (auto trig : groupedTrigger) {
              //Only fetch pot for beam trigger. 
              if (!trig.first->IsBeamTrigger(trig.second.triggerType)) return 0;
              bool matched = std::any_of(data.cbegin(), data.cend(),
                  [par, &trig](const auto& spill) {
                      return std::abs(util::getTriggerTime(trig.second) - spill.first) < par().cafmaker().beamMatchDT();
                  });
              if (matched) {
                  any_matched = true;
                  matched_triggers.push_back(trig);
              } else {
                  unmatched_triggers.push_back(trig);
              }
          }
  
          std::stringstream log_message;
  
          if (any_matched) {
              log_message << "Only some triggers match beam spill for trigger group " << ii << ":\n"
                          << "Matched triggers: \n";
              for (auto trig : matched_triggers)
                  log_message << std::fixed << trig.first->GetName() << " " << util::getTriggerTime(trig.second) <<" "<<trig.second.triggerType<< "\n";
  
              log_message << "Unmatched triggers: \n";
              for (auto trig : unmatched_triggers)
                  log_message << std::fixed << trig.first->GetName() << " " << util::getTriggerTime(trig.second)<<" "<<trig.second.triggerType << "\n";
              cafmaker::LOG_S("Fill beam POT information").ERROR() << log_message.str() << "\n";
              std::abort();

          } 
          else {
              log_message << "No matching spill found for trigger group " << ii << " with triggers: \n";
              for (auto trig : groupedTrigger)
                  log_message << std::fixed << trig.first->GetName() << " " << util::getTriggerTime(trig.second) << "\n";
              cafmaker::LOG_S("Fill beam POT information").WARNING() << log_message.str() << "\n";
          }
      }
  
     return pot;
  }

  std::vector<double> IFBeam::getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const BeamInfo& data) {
      std::cout << "function " << __func__ << ", line " << __LINE__ << "\n";

      std::vector<double> values;
      if (!(groupedTrigger.front().first->IsBeamTrigger(groupedTrigger.front().second.triggerType))) return {0.0};
      
      std::cout << "function " << __func__ << ", line " << __LINE__ << "\n";

      auto it = std::find_if(data.begin(), data.end(),
          [par, &groupedTrigger](const auto& spill) {
              return std::all_of(groupedTrigger.cbegin(), groupedTrigger.cend(),
                  [par, &spill](const auto& groupedTrigger) {
                      return std::abs(util::getTriggerTime(groupedTrigger.second) - spill.first) < par().cafmaker().beamMatchDT();
                  });
          });
  
      std::cout << "function " << __func__ << ", line " << __LINE__ << "\n";

      if (it != data.end()) {
          values = it->second;
      } 
      else {
          std::cout << "function " << __func__ << ", line " << __LINE__ << "\n";
          bool any_matched = false;
          std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> matched_triggers;
          std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> unmatched_triggers;
          for (auto trig : groupedTrigger) {
              //Only fetch values for beam trigger. 
              if (!trig.first->IsBeamTrigger(trig.second.triggerType)) return {0.0};
              bool matched = std::any_of(data.cbegin(), data.cend(),
                  [par, &trig](const auto& spill) {
                      return std::abs(util::getTriggerTime(trig.second) - spill.first) < par().cafmaker().beamMatchDT();
                  });
              if (matched) {
                  any_matched = true;
                  matched_triggers.push_back(trig);
              } else {
                  unmatched_triggers.push_back(trig);
              }
          }
  
          std::cout << "function " << __func__ << ", line " << __LINE__ << "\n";
          std::stringstream log_message;
  
          if (any_matched) {
              log_message << "Only some triggers match beam spill for trigger group " << ii << ":\n"
                          << "Matched triggers: \n";
              for (auto trig : matched_triggers)
                  log_message << std::fixed << trig.first->GetName() << " " << util::getTriggerTime(trig.second) <<" "<<trig.second.triggerType<< "\n";
  
              log_message << "Unmatched triggers: \n";
              for (auto trig : unmatched_triggers)
                  log_message << std::fixed << trig.first->GetName() << " " << util::getTriggerTime(trig.second)<<" "<<trig.second.triggerType << "\n";
              cafmaker::LOG_S("Fill beam POT information").ERROR() << log_message.str() << "\n";
              std::abort();

          } 
          else {
              log_message << "No matching spill found for trigger group " << ii << " with triggers: \n";
              for (auto trig : groupedTrigger)
                  log_message << std::fixed << trig.first->GetName() << " " << util::getTriggerTime(trig.second) << "\n";
              cafmaker::LOG_S("Fill beam POT information").WARNING() << log_message.str() << "\n";
          }
      }
  
     return values;
  }

  double IFBeam::getPOT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, beamSpills);
  }
 
  double IFBeam::getHornI(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, hornI);
  }

   std::vector<double> IFBeam::getHorizontalPosTGT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, horizontalPosTGT);
  } 

 std::vector<double> IFBeam::getHorizontalIntTGT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, horizontalIntTGT);
  } 

  std::vector<double> IFBeam::getHorizontalPos121(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, horizontalPos121);
  } 

  std::vector<double> IFBeam::getVerticalIntTGT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, verticalIntTGT);
  } 

  std::vector<double> IFBeam::getVerticalPosTGT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, verticalPosTGT);
  } 

  std::vector<double> IFBeam::getVerticalPos121(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, verticalPos121);
  } 

  std::vector<double> IFBeam::getMultiWireInfo(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
    return getData(par, groupedTrigger, ii, multiwireInfo);
  } 
}
