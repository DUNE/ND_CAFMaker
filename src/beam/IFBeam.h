#ifndef IFBEAM_H
#define IFBEAM_H

#include <string>
#include <map>
#include <vector>
#include <nlohmann/json.hpp>

#include "Params.h"
#include "CAF.h"
#include "reco/IRecoBranchFiller.h"


class IFBeam {
public:
    using BeamSpills = std::map<double, double>;
    IFBeam(const std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>>& groupedTriggers, bool is_data);
    
    double getPOT(const cafmaker::Params& par, std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>& groupedTrigger, int ii);

private:
    const std::string potDevice = "E:TRTGTD";
    BeamSpills beamSpills;

    void loadBeamSpills(const std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>>& groupedTriggers);
    std::string createUrl(const std::string& min_time_iso, const std::string& max_time_iso);
    double unitToFactor(const std::string& unit);
};
 
std::string toISO8601(double time_sec);
double getTriggerTime(const cafmaker::Trigger& trigger);

#endif // IFBEAM_H
