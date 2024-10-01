#include "util/Logger.h"
#include "util/Progress.h"
#include "IFBeam.h"
#include <curl/curl.h>
#include <chrono>
#include <regex>
#include <sstream>
#include <iomanip>
#include <limits>

size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* userp) { //write https response data into a string
    userp->append((char*)contents, size * nmemb);
    return size * nmemb;
}


double getTriggerTime(const cafmaker::Trigger& trigger) {
    return trigger.triggerTime_s + 1e-9 * trigger.triggerTime_ns;
}

std::string toISO8601(double time_sec) { //IFBeam query requires ISO format
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

IFBeam::IFBeam(const std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>>& groupedTriggers, bool is_data) {
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
        std::cerr << "Unknown unit in beam database: " << unit << std::endl;
        return 1.0;
    }
}

void IFBeam::loadBeamSpills(const std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>>& groupedTriggers) { //todo: this should return other information as well like horn current, position, etc., querying all devices and storing in a map
    beamSpills.clear();
    double dt = 5.0; //time window to query before and after first and last spill, respectively
    double ms_to_s = 1e-3;
    double min_time = std::numeric_limits<double>::max();
    double max_time = std::numeric_limits<double>::lowest();
    for (const auto& group : groupedTriggers) {
        for (const auto& trig : group) {
            double trigger_time = getTriggerTime(trig.second);
            min_time = std::min(min_time, trigger_time - dt);
            max_time = std::max(max_time, trigger_time + dt);
        }
    }
    std::string min_time_iso = toISO8601(min_time);
    std::string max_time_iso = toISO8601(max_time);
    std::string url = createUrl(min_time_iso, max_time_iso);

    
    CURL* curl;
    CURLcode res;
    std::string readBuffer;

    std::cout << "Fetching beam information from: " << url << "\n";
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
                double time = spill["clock"].get<double>() * ms_to_s;
                std::string unit = spill["units"];
                double pot = spill["value"].get<double>() * unitToFactor(unit);
                beamSpills[time] = pot;
            }
        } catch (const json::exception& e) {
            std::cerr << "Failed to parse JSON data: " << e.what() << std::endl;
        }
    }
}

double IFBeam::getPOT(const cafmaker::Params& par, std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>& groupedTrigger, int ii) {
    double pot = 0.0;

    auto it = std::find_if(beamSpills.begin(), beamSpills.end(),
        [par, &groupedTrigger](const auto& spill) {
            return std::all_of(groupedTrigger.cbegin(), groupedTrigger.cend(),
                [par, &spill](const auto& groupedTrigger) {
                    return std::abs(getTriggerTime(groupedTrigger.second) - spill.first) < par().cafmaker().beamMatchDT();
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
            bool matched = std::any_of(beamSpills.cbegin(), beamSpills.cend(),
                [par, &trig](const auto& spill) {
                    return std::abs(getTriggerTime(trig.second) - spill.first) < par().cafmaker().beamMatchDT();
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
                log_message << std::fixed << trig.first->GetName() << " " << getTriggerTime(trig.second) << "\n";

            log_message << "Unmatched triggers: \n";
            for (auto trig : unmatched_triggers)
                log_message << std::fixed << trig.first->GetName() << " " << getTriggerTime(trig.second) << "\n";

            LOG().ERROR() << log_message.str() << "\n";
            std::abort();
        } 
        else {
            log_message << "No matching spill found for trigger group " << ii << " with triggers: \n";
            for (auto trig : groupedTrigger)
                log_message << std::fixed << trig.first->GetName() << " " << getTriggerTime(trig.second) << "\n";
            LOG().WARNING() << log_message.str() << "\n";
        }
    }

   return pot;
}
