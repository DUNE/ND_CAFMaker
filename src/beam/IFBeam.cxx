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
#include <regex>
#include <sstream>
#include <limits>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace cafmaker
{

  
  IFBeam::IFBeam(const std::vector<TriggerGroup>& groupedTriggers, bool is_data) {
      if (is_data)
  	loadBeamSpills(groupedTriggers); // Load all beam spills upon instantiation if data
  }
  
  std::string IFBeam::createUrl(const std::string& min_time_iso, const std::string& max_time_iso) {
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
  
  void IFBeam::loadBeamSpills(const std::vector<TriggerGroup>& groupedTriggers) { //todo: this should return other information as well like horn current, position, etc., querying all devices and storing in a map
      beamSpills.clear();
      double dt = 5.0; //time window to query before and after first and last spill, respectively
      double ms_to_s = 1e-3;
      double min_time = std::numeric_limits<double>::max();
      double max_time = std::numeric_limits<double>::lowest();
      for (const auto& group : groupedTriggers) {
          for (const auto& trig : group) {
              if (trig.second.triggerType !=1) continue;
              double trigger_time = util::getTriggerTime(trig.second);

              // HACK
              if (trigger_time > 1e12)
                continue;
              min_time = std::min(min_time, trigger_time - dt);
              max_time = std::max(max_time, trigger_time + dt);
          }
      }
      std::string min_time_iso = util::toISO8601(min_time);
      std::string max_time_iso = util::toISO8601(max_time);
      std::string url = createUrl(min_time_iso, max_time_iso);
  
      
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
              for (const auto& spill : json_data["rows"]) {
                  double time = spill["clock"].get<double>() * ms_to_s;
                  std::string unit = spill["units"];
                  double pot = spill["value"].get<double>() * unitToFactor(unit);
                  beamSpills[time] = pot;
              }
          } catch (const json::exception& e) {
              cafmaker::LOG_S("Fetch beam information").ERROR() << "Failed to parse JSON data: " << e.what() << "\n";
  	    std::abort();
          }
      }
  }
  
  double IFBeam::getPOT(const cafmaker::Params& par, const TriggerGroup& groupedTrigger, int ii) {
      double pot = 0.0;
      if (groupedTrigger.front().second.triggerType != 1) return 0.0;
   //   for (auto trig : groupedTrigger) {std::cout<<trig.second.triggerType<<" ";}
  //  std::cout<<std::endl;
      auto it = std::find_if(beamSpills.begin(), beamSpills.end(),
          [par, &groupedTrigger](const auto& spill) {
              return std::all_of(groupedTrigger.cbegin(), groupedTrigger.cend(),
                  [par, &spill](const auto& groupedTrigger) {
                      return std::abs(util::getTriggerTime(groupedTrigger.second) - spill.first) < par().cafmaker().beamMatchDT();
                  });
          });
  
      if (it != beamSpills.end()) {
          pot = it->second;
      } 
      else {
          bool any_matched = false;
          std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> matched_triggers;
          std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> unmatched_triggers;
          for (auto trig : groupedTrigger) {
              //check if the trigger is a type 1 trigger (beam trigger)
              if (trig.second.triggerType != 1) return 0;
              bool matched = std::any_of(beamSpills.cbegin(), beamSpills.cend(),
                  [par, &trig](const auto& spill) {
                      return std::abs(util::getTriggerTime(trig.second) - spill.first) < par().cafmaker().beamMatchDT();
                  });
              std::cout<<std::endl<<std::endl; 
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
             /* log_message << "No matching spill found for trigger group " << ii << " with triggers: \n";
              for (auto trig : groupedTrigger)
                  log_message << std::fixed << trig.first->GetName() << " " << util::getTriggerTime(trig.second) << "\n";
              cafmaker::LOG_S("Fill beam POT information").WARNING() << log_message.str() << "\n";*/
          }
      }
  
     return pot;
  }
}
