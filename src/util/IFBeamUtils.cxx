/// \file IFBeamUtils.h
///
/// Utility functions for IFBeam database
///
/// \author  S. Kumaran <s.kumaran@uci.edu>
/// \date    Oct. 2024

#include <util/IFBeamUtils.h>
#include <chrono>
#include <sstream>
#include <iomanip>

namespace cafmaker
{
  namespace util
  {
    size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* userp) { //write https response data into a string
        userp->append((char*)contents, size * nmemb);
        return size * nmemb;
    }
    
    
    double getTriggerTime(const Trigger& trigger) {
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
  }
}
