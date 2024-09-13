#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <regex>

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/parse.h"


#include "CAF.h"
#include "Params.h"
#include "reco/IRecoBranchFiller.h"
#include "util/Logger.h"
#include "util/Progress.h"
#include "IFBeam.h"


std::string pot_prim_device = "E:TRTGTD";


size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* userp){  //write https response data into a string
    userp->append((char*)contents, size * nmemb);
    return size * nmemb;
}

double unitToFactor(const std::string& unit) { // convert the "units" (E14) field in IFBeam database to exponential
    std::regex re("E(\\d+)");
    std::smatch match;

    if (std::regex_search(unit, match, re) && match.size() > 1) {
        int exponent = std::stoi(match.str(1));
        return std::pow(10, exponent);
    } else {
        std::cerr << "Unknown or malformed unit in beam database: " << unit << std::endl;
        return 1.0; 
    }

}

double GetTriggerTime(const cafmaker::Trigger& trigger) {
    return trigger.triggerTime_s + 1e-9 * trigger.triggerTime_ns;
}

std::string CreateUrl(std::string min_time_iso, std::string max_time_iso, std::string device_name){

    std::ostringstream url_stream;
    url_stream << "https://dbdata3vm.fnal.gov:9443/ifbeam/data/data?v=" << device_name
               << "&e=" << "e,a9"
               << "&t0=" << min_time_iso
               << "&t1=" << max_time_iso
               << "&f=json";

   std::string url = url_stream.str();

   return url;
}

std::string toISO8601(double time_sec) { // the IFBeam request only takes isoformat, so the unix time should be converted to iso format
    auto tp = std::chrono::system_clock::from_time_t(static_cast<time_t>(time_sec));
    auto in_time_t = std::chrono::system_clock::to_time_t(tp);
    
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%dT%H:%M:%S");

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()) % 1000;
    if (ms.count() > 0) {
        ss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    }
    ss << (std::localtime(&in_time_t)->tm_gmtoff >= 0 ? '+' : '-') 
       << std::setfill('0') << std::setw(2) << std::abs(std::localtime(&in_time_t)->tm_gmtoff) / 3600
       << ':' << std::setfill('0') << std::setw(2) << (std::abs(std::localtime(&in_time_t)->tm_gmtoff) % 3600) / 60;

    return ss.str();
}

BeamSpills loadBeamSpills(const std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>>& groupedTriggers) { //load beam spills using https request    
    BeamSpills beam_spills;
    double dt = 5.0; // +/- 5 seconds
    double min_time = std::numeric_limits<double>::max();
    double max_time = std::numeric_limits<double>::lowest();

    for (const auto& group : groupedTriggers) {
        for (const auto& trig : group) {
            double trigger_time = GetTriggerTime(trig.second);
            min_time = std::min(min_time, trigger_time - dt);
            max_time = std::max(max_time, trigger_time + dt);
        }
    }

    std::string min_time_iso = toISO8601(min_time);
    std::string max_time_iso = toISO8601(max_time);

    std::string url = CreateUrl(min_time_iso, max_time_iso, pot_prim_device);
     
    std::cout << "Loading beam spills from " << url << std::endl; 
    CURL* curl;
    CURLcode res;
    std::string readBuffer;

    curl = curl_easy_init();
    if (curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        res = curl_easy_perform(curl);
        if (res != CURLE_OK) {
            std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
        }
        curl_easy_cleanup(curl);
        try { 
            auto json_data = json::parse(readBuffer);
            for (const auto& spill : json_data["rows"]) {
                double time = spill["clock"].get<double>() / 1000.0; 
		std::string unit = spill["units"];
                double pot = spill["value"].get<double>()* unitToFactor(unit);
                beam_spills[time] = pot;
            }
        } catch (const json::exception& e) {
            std::cerr << "Failed to parse JSON data: " << e.what() << std::endl;
        }
    }

    return beam_spills;
}


    
double getPOT(const cafmaker::Params& par, std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>& groupedTrigger, const BeamSpills& beam_spills, int ii) {
    double pot = 0.0;

    auto it = std::find_if(beam_spills.begin(), beam_spills.end(),
        [par, &groupedTrigger](const auto& spill) {
            return std::all_of(groupedTrigger.cbegin(), groupedTrigger.cend(),
                [par, &spill](const auto& groupedTrigger) {
                    return std::abs(GetTriggerTime(groupedTrigger.second) - spill.first) < par().cafmaker().beamMatchDT();
                });
        });

    if (it != beam_spills.end()) {
        pot = it->second;
    } else {
        bool any_matched = false;
        std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> matched_triggers;
        std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> unmatched_triggers;

        for (auto trig : groupedTrigger) {
            bool matched = std::any_of(beam_spills.cbegin(), beam_spills.cend(),
                [par, &trig](const auto& spill) {
                    return std::abs(GetTriggerTime(trig.second) - spill.first) < par().cafmaker().beamMatchDT();
                });

            if (matched) {
                any_matched = true;
                matched_triggers.push_back(trig);
            } else {
                unmatched_triggers.push_back(trig);
            }
        }

        auto LOG = [&]() -> const cafmaker::Logger & { return cafmaker::LOG_S("Beam spill matching"); };
        std::stringstream log_message;

        if (any_matched) {
            log_message << "Only some triggers match beam spill for trigger group " << ii << ":\n"
                        << "Matched triggers: \n";
            for (auto trig : matched_triggers)
                log_message << std::fixed << trig.first->GetName() << " " << GetTriggerTime(trig.second) << "\n";

            log_message << "Unmatched triggers: \n";
            for (auto trig : unmatched_triggers)
                log_message << std::fixed << trig.first->GetName() << " " << GetTriggerTime(trig.second) << "\n";

            LOG().ERROR() << log_message.str() << "\n";
            std::abort();
        } else {
            log_message << "No matching spill found for trigger group " << ii << " with triggers: \n";
            for (auto trig : groupedTrigger)
                log_message << std::fixed << trig.first->GetName() << " " << GetTriggerTime(trig.second) << "\n";
            LOG().WARNING() << log_message.str() << "\n";
        }
    }

    return pot;
}

