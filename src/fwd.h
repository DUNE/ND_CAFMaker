/// \file fwd.h
///
/// Forward declarations that are annoying to write every time
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    Jan. 2022

#ifndef ND_CAFMAKER_FWD_H
#define ND_CAFMAKER_FWD_H

class TTree;

// forward-declaring the fhiclcpp stuff
// takes some gyrations since they're templated
namespace fhicl
{
  template <typename T, typename KeysToIgnore>
  class Table;
}

namespace caf
{
  class StandardRecord;
}

namespace cafmaker
{
  class dumpTree;

  class FhiclConfig;
  using Params = fhicl::Table<cafmaker::FhiclConfig, void>;
}

namespace genie
{
  class NtpMCEventRecord;
}

#endif //ND_CAFMAKER_FWD_H
