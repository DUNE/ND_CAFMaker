/// \file FillTruthForSAND.h
///
/// Fill truth branches for SAND data.
///
/// \author  L.Di Noto based on code by J. Wolcott <jwolcott@fnal.gov> and C. Marshall <chris.marshall@rochester.edu>
/// \date    Apr. 2022

#ifndef ND_CAFMAKER_FILLTRUTHFORSAND_H
#define ND_CAFMAKER_FILLTRUTHFORSAND_H

#include "fwd.h"
#include "struct.h"
// fixme: this will need to be put back to the actual response_helper type when DIRT-II finishes model recommendations
#include <string>
namespace nusyst                      //from NDlar version
{
  using response_helper = std::string;
}

void fillTruthForSAND(caf::StandardRecord& sr,
               TTree *intree,
               TTree * gtree,
               event * evt,
               const genie::NtpMCEventRecord * mcrec,
               const cafmaker::Params &par,
               nusyst::response_helper& rh);

#endif //ND_CAFMAKER_FILLTRUTHSAND_H
