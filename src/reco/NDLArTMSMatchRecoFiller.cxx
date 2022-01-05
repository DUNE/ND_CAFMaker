#include "NDLArTMSMatchRecoFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRNDTrackAssn.h"
#include "MLNDLArRecoBranchFiller.h"
#include "TMSRecoBranchFiller.h"


namespace cafmaker
{
  NDLArTMSMatchRecoFiller::NDLArTMSMatchRecoFiller()
  {
    // nothing to do
    SetConfigured(true);
  }

  void NDLArTMSMatchRecoFiller::_FillRecoBranches(std::size_t evtIdx,
                                                  caf::StandardRecord &sr,
                                                  const cafmaker::dumpTree &dt,
                                                  const cafmaker::Params &par) const
  {
    // match tracks using the info that should have been filled by the ND-LAr and TMS reco fillers
    MatchTracks(sr);
  }

// -----------------------------------------------------

  void NDLArTMSMatchRecoFiller::MatchTracks(caf::StandardRecord &sr) const
  {
    // actually do track matching here:
    //   compare impact parameters, angles, whatever


    // then, fill in the ndmatch object with the matches.
    // for example (obviously you'd use calculated vaoues, not hard-coded ones...):
    caf::SRNDTrackAssn match;
    match.larid = 4;
    match.tmsid = 3;
    match.transdispl = 1.865; // whatever
    match.angdispl = 8.4;  // whatever

    sr.nd.trkmatch.emplace_back(match);
    sr.nd.ntrkmatch++;
  }

}