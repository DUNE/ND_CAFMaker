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

bool loadBeamSpills(const std::string& filename, BeamSpills& beam_spills);

double getPOT(const cafmaker::Params& par, 
              std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>& groupedTrigger, 
              int ii);

#endif // IFBEAM_H
