/// \file FillTruth.h
///
/// Fill truth branches.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>, based on code by C. Marshall <chris.marshall@rochester.edu>
/// \date    Jan. 2022

#ifndef ND_CAFMAKER_FILLTRUTH_H
#define ND_CAFMAKER_FILLTRUTH_H

#include "fwd.h"

// fixme: this will need to be put back to the actual response_helper type when DIRT-II finishes model recommendations
#include <string>
namespace nusyst
{
  using response_helper = std::string;
}

void fillTruth(int ii,
	       caf::StandardRecord& sr,
               TTree * gtree,
               const genie::NtpMCEventRecord * mcrec,
               const cafmaker::Params &par,
               nusyst::response_helper& rh);

#endif //ND_CAFMAKER_FILLTRUTH_H
