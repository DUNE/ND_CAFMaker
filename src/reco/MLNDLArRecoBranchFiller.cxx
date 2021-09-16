//
// Created by jeremy on 9/10/21.
//

#include "MLNDLArRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "Params.h"

namespace cafmaker
{

  // ------------------------------------------------------------------------------
  MLNDLArRecoBranchFiller::MLNDLArRecoBranchFiller(const std::string &h5filename,
                                                   const std::string &h5dataset)
    : fH5file(h5filename, h5dataset)
  {
    // if we got this far, nothing bad happened trying to open the file or dataset
    SetConfigured(true);
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::_FillRecoBranches(caf::StandardRecord &sr, const cafmaker::dumpTree &dt, const cafmaker::params &par) const
  {

    for (std::size_t evtIdx = 0; evtIdx < par.n; evtIdx++)
    {
      std::vector<cafmaker::Track> recoTracks = fH5file.EventTracks(evtIdx);
      sr.ndlar.ntracks = recoTracks.size();
      for (const auto & track : recoTracks)
      {
        caf::SRTrack srTrack;
        srTrack.start = track.start;
        srTrack.end = track.end;
        srTrack.end_dir = track.end_dir;
        sr.ndlar.tracks.emplace_back(std::move(srTrack));
      } // for (track)
    } // for (evtIdx)

  } // MLNDLArRecoBranchFiller::_FillRecoBranches()

} // namespace cafmaker