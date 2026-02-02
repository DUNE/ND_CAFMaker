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
          return 1.0;
      }
  }

  IFBeam::BeamInfo IFBeam::retrieveInfoFromDataBase(const std::string url) {

      cafmaker::LOG_S("Fetchig beam info from url : ").VERBOSE() << url << "\n";

      double ms_to_s = 1e-3;

      BeamInfo data;
 
      CURL* curl;
      CURLcode res;
      std::string readBuffer; 
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
              }
          } catch (const json::exception& e) {
              cafmaker::LOG_S("Fetch beam information").ERROR() << "Failed to parse JSON data: " << e.what() << "\n";
  	    std::abort();
          }
      }
      return data;
  }

  void IFBeam::loadBeamSpills(const cafmaker::Params& par, const std::vector<TriggerGroup>& groupedTriggers) { 

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

      for (auto& device : deviceMap)
      {
        device.second.clear();
        std::string url_device = createUrl(device.first, min_time_iso, max_time_iso);
        device.second = retrieveInfoFromDataBase(url_device);
      }

      // see: https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/IFDBSpillInfo_module.cc#L685
      // the horn current is calculated as I_horn = (E:NSLINA-(+0.01))/0.9951 + (E:NSLINB-(-0.14))/0.9957 + (E:NSLINC-(-0.05))/0.9965 + (E:NSLIND-(-0.07))/0.9945
      
      double currentA=0.0, currentB=0.0, currentC=0.0, currentD=0.0;

      for(const auto& pairA : deviceMap["E:NSLINA"]) {
        auto timeA = pairA.first;
        auto currentA = pairA.second.at(0); 
        
        if(currentA == 0.) continue;

        for(const auto& pairB : deviceMap["E:NSLINB"]) {
          auto timeB = pairB.first;
          if (abs(timeA - timeB) <= par().cafmaker().beamMatchDT()) {
            currentB = pairB.second.at(0);
            break;
          }
        }

        if(currentB == 0.) continue;

        for(const auto& pairC : deviceMap["E:NSLINC"]) {
          auto timeC = pairC.first;
          if (abs(timeA - timeC) <= par().cafmaker().beamMatchDT()) 
          {
            currentC = pairC.second.at(0);
            break;
          }
        }

        if(currentC == 0.) continue;

        for(const auto& pairD : deviceMap["E:NSLIND"]) {
          auto timeD = pairD.first;
          if (abs(timeA - timeD) <= par().cafmaker().beamMatchDT()){
            currentD = pairD.second.at(0);
            break;
          }
        }

        if(currentD == 0.) continue;
        
        deviceMap["E:NSLIN"][timeA] = {(currentA-0.01)/0.9951 + (currentB+0.14)/0.9957 + (currentC+0.05)/0.9965 + (currentD+0.07)/0.9945};
      }

  }

  std::vector<double> IFBeam::getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const BeamInfo& data) {

      std::vector<double> values;
      if (!(groupedTrigger.front().first->IsBeamTrigger(groupedTrigger.front().second.triggerType))) return {0.0};
      
      auto it = std::find_if(data.begin(), data.end(),
          [par, &groupedTrigger](const auto& spill) {
              return std::all_of(groupedTrigger.cbegin(), groupedTrigger.cend(),
                  [par, &spill](const auto& groupedTrigger) {
                      return std::abs(util::getTriggerTime(groupedTrigger.second) - spill.first) < par().cafmaker().beamMatchDT();
                  });
          });

      if (it != data.end()) {
          values = it->second;
      } 
      else {
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

  std::vector<double> IFBeam::getData(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii, const std::string& deviceName) {

    if(deviceMap.count(deviceName)>0)
      return getData(par, groupedTrigger, ii, deviceMap[deviceName]);

    cafmaker::LOG_S("Trying retrieving info from unkwnown device ").WARNING() << deviceName << ", returning empty vector\n";

    return {};
  }
}
