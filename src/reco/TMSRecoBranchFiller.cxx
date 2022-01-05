#include "TMSRecoBranchFiller.h"

namespace cafmaker
{
  TMSRecoBranchFiller::TMSRecoBranchFiller(const std::string &tmsRecoFilename)
    : fTMSRecoFile(tmsRecoFilename.c_str())
  {
    if (!fTMSRecoFile.IsZombie())
      SetConfigured(true);
  }

  // ---------------------------------------------------------------------------

  void TMSRecoBranchFiller::_FillRecoBranches(std::size_t evtIdx,
                                              caf::StandardRecord &sr,
                                              const cafmaker::dumpTree &dt,
                                              const cafmaker::Params &par) const
  {
    // here we copy all the TMS reco into the SRTMS branch of the StandardRecord object.

  }

}