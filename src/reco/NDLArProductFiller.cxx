#include "NDLArProductFiller.h"

#include <iostream>

#include "duneanaobj/StandardRecord/StandardRecord.h"

namespace cafmaker
{

  // -------------------------------------------------------------

  template <>
  void NDLArProductFiller<caf::SRTrack>::FillSR(caf::StandardRecord &sr, std::size_t evtIdx) const
  {
    //std::cout << "Filling tracks for event " << evtIdx << std::endl;
    std::vector<caf::SRTrack> recoTracks = EventProducts(evtIdx);
    sr.nd.lar.ntracks = recoTracks.size();
    sr.nd.lar.tracks = std::move(recoTracks);
  }

  // -------------------------------------------------------------

  template <>
  std::vector<caf::SRTrack>
  NDLArProductFiller<caf::SRTrack>::ParseVals(const float * buffer, std::size_t nRows, std::size_t nCols)
  {
    std::vector<caf::SRTrack> tracksOut;
    tracksOut.reserve(nRows);

    // avoid looking these up every time since they won't change
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_start_x)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_start_y)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_start_z)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_end_x)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_end_y)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_end_z)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_start_dir_x)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_start_dir_y)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_start_dir_z)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_end_dir_x)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_end_dir_y)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_end_dir_z)
    BUFFER_LOOKUP_VAR(caf::SRTrack, trk_visE)

    for (std::size_t rowIdx = 0; rowIdx < nRows; rowIdx++)
    {
      caf::SRTrack tr;
      std::size_t rowOffset = rowIdx * nCols;

      // just build the offset in ...
      const float * offBuf = &buffer[rowOffset];

      tr.start = {offBuf[trk_start_x], offBuf[trk_start_y], offBuf[trk_start_z]};
      tr.end = {offBuf[trk_end_x], offBuf[trk_end_y], offBuf[trk_end_z]};
      tr.dir = {offBuf[trk_start_dir_x], offBuf[trk_start_dir_y], offBuf[trk_start_dir_z]};
      tr.enddir = {offBuf[trk_end_dir_x], offBuf[trk_end_dir_y], offBuf[trk_end_dir_z]};
      tr.Evis = offBuf[trk_visE];

      //std::cout << "  filling track: " << tr << std::endl;

      tracksOut.emplace_back(std::move(tr));
    }

    return tracksOut;
  }

  // -------------------------------------------------------------
  // -------------------------------------------------------------

  template <>
  void NDLArProductFiller<caf::SRShower>::FillSR(caf::StandardRecord &sr, std::size_t evtIdx) const
  {
    //std::cout << "Filling showers for event " << evtIdx << std::endl;
    std::vector<caf::SRShower> recoShw = EventProducts(evtIdx);
    sr.nd.lar.nshowers = recoShw.size();
    sr.nd.lar.showers = std::move(recoShw);
  }

  // -------------------------------------------------------------

  template <>
  std::vector<caf::SRShower>
  NDLArProductFiller<caf::SRShower>::ParseVals(const float * buffer, std::size_t nRows, std::size_t nCols)
  {
    std::vector<caf::SRShower> shwOut;
    shwOut.reserve(nRows);

    BUFFER_LOOKUP_VAR(caf::SRShower, shw_start_x)
    BUFFER_LOOKUP_VAR(caf::SRShower, shw_start_y)
    BUFFER_LOOKUP_VAR(caf::SRShower, shw_start_z)
    BUFFER_LOOKUP_VAR(caf::SRShower, shw_dir_x)
    BUFFER_LOOKUP_VAR(caf::SRShower, shw_dir_y)
    BUFFER_LOOKUP_VAR(caf::SRShower, shw_dir_z)
    BUFFER_LOOKUP_VAR(caf::SRShower, shw_visE)

    for (std::size_t rowIdx = 0; rowIdx < nRows; rowIdx++)
    {
      caf::SRShower shw;
      std::size_t rowOffset = rowIdx * nCols;

      // just build the offset in ...
      const float * offBuf = &buffer[rowOffset];

      shw.start = {offBuf[shw_start_x], offBuf[shw_start_y], offBuf[shw_start_z]};
      shw.direction = {offBuf[shw_dir_x], offBuf[shw_dir_y], offBuf[shw_dir_z]};
      shw.Evis = offBuf[shw_visE];

      shwOut.emplace_back(std::move(shw));
    }

    return shwOut;
  }

}
