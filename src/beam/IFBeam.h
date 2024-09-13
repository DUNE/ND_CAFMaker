#ifndef IFBEAM_H
#define IFBEAM_H

#include <string>
#include <map>
#include <vector>
#include <nlohmann/json.hpp>

#include "Params.h"
#include "CAF.h"

using BeamSpills = std::map<double, double>;
using json = nlohmann::json;

BeamSpills loadBeamSpills(const std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>>& groupedTriggers);   

double getPOT(const cafmaker::Params& par, 
              std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>& groupedTrigger,
              const BeamSpills& beam_spills, 
              int ii);

#endif // IFBEAM_H
