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
  void
  MLNDLArRecoBranchFiller::_FillRecoBranches(std::size_t evtIdx,
                                             caf::StandardRecord &sr,
                                             const cafmaker::dumpTree &dt,
                                             const cafmaker::params &par) const
  {
    fTrackFiller.FillSR(sr, evtIdx);
    fShowerFiller.FillSR(sr, evtIdx);
  } // MLNDLArRecoBranchFiller::_FillRecoBranches()

} // namespace cafmaker