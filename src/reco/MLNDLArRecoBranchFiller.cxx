#include "MLNDLArRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "Params.h"

namespace cafmaker
{

  // ------------------------------------------------------------------------------
  // todo: possibly build some mechanism for customizing the dataset names in the file here
  MLNDLArRecoBranchFiller::MLNDLArRecoBranchFiller(const std::string &h5filename)
    : fTrackFiller(h5filename, "tracks"), fShowerFiller(h5filename, "showers")
  {
    // if we got this far, nothing bad happened trying to open the file or dataset
    SetConfigured(true);
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::_FillRecoBranches(caf::StandardRecord &sr, const cafmaker::dumpTree &dt, const cafmaker::params &par) const
  {
    assert(par.n >= 0);
    for (std::size_t evtIdx = 0; evtIdx < static_cast<std::size_t>(par.n); evtIdx++)
    {
      fTrackFiller.FillSR(sr, evtIdx);
      fShowerFiller.FillSR(sr, evtIdx);
     } // for (evtIdx)

  } // MLNDLArRecoBranchFiller::_FillRecoBranches()

} // namespace cafmaker