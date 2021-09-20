#include "NDLArProductFiller.h"

#include <iostream>

#include "duneanaobj/StandardRecord/StandardRecord.h"

namespace cafmaker
{
  template <>
  const std::vector<std::string> EXPECTED_COLUMN_NAMES<caf::SRTrack>
  {
      "trk_start_x",
      "trk_start_y",
      "trk_start_z",
      "trk_end_x",
      "trk_end_y",
      "trk_end_z",
      "trk_end_dir_x",
      "trk_end_dir_y",
      "trk_end_dir_z",
      "trk_visE"
  };

  // -------------------------------------------------------------
  template <>
  const std::vector<std::string> EXPECTED_COLUMN_NAMES<caf::SRShower>
  {
      "shw_start_x",
      "shw_start_y",
      "shw_start_z",
      "shw_dir_x",
      "shw_dir_y",
      "shw_dir_z",
      "shw_visE"
  };

  // -------------------------------------------------------------

  template <>
  void NDLArProductFiller<caf::SRTrack>::FillSR(caf::StandardRecord &sr, std::size_t evtIdx) const
  {
    std::vector<caf::SRTrack> recoTracks = EventProducts(evtIdx);
    sr.ndlar.ntracks = recoTracks.size();
    sr.ndlar.tracks = std::move(recoTracks);
  }

  // -------------------------------------------------------------

  template <>
  std::vector<caf::SRTrack>
  NDLArProductFiller<caf::SRTrack>::ParseVals(const float * buffer, std::size_t nRows, std::size_t nCols)
  {
    std::vector<caf::SRTrack> tracksOut;
    tracksOut.reserve(nRows);

    for (std::size_t rowIdx = 0; rowIdx < nRows; rowIdx++)
    {
      caf::SRTrack tr;
      std::size_t rowOffset = rowIdx * nCols;

      // n.b.: if you change anything here,
      //       make SURE the column assumptions matchs the checks in ValidateColumns()
      tr.start = {buffer[rowOffset + 0], buffer[rowOffset + 1], buffer[rowOffset + 2]};
      tr.end = {buffer[rowOffset + 3], buffer[rowOffset + 4], buffer[rowOffset + 5]};
      tr.end_dir = {buffer[rowOffset + 6], buffer[rowOffset + 7], buffer[rowOffset + 8]};
      tr.Evis = buffer[rowOffset+9];

      tracksOut.emplace_back(std::move(tr));
    }

    return tracksOut;
  }

  // -------------------------------------------------------------
  // -------------------------------------------------------------

  template <>
  void NDLArProductFiller<caf::SRShower>::FillSR(caf::StandardRecord &sr, std::size_t evtIdx) const
  {
    std::vector<caf::SRShower> recoShw = EventProducts(evtIdx);
    sr.ndlar.nshowers = recoShw.size();
    sr.ndlar.showers = std::move(recoShw);
  }

  // -------------------------------------------------------------

  template <>
  std::vector<caf::SRShower>
  NDLArProductFiller<caf::SRShower>::ParseVals(const float * buffer, std::size_t nRows, std::size_t nCols)
  {
    std::vector<caf::SRShower> shwOut;
    shwOut.reserve(nRows);

    for (std::size_t rowIdx = 0; rowIdx < nRows; rowIdx++)
    {
      caf::SRShower shw;
      std::size_t rowOffset = rowIdx * nCols;

      // todo: probably should come up with a scheme for validating
      //       that these offsets really match up with the data we think they do
      shw.start = {buffer[rowOffset + 0], buffer[rowOffset + 1], buffer[rowOffset + 2]};
      shw.direction = {buffer[rowOffset + 3], buffer[rowOffset + 4], buffer[rowOffset + 5]};
      shw.Evis = buffer[rowOffset + 6];

      shwOut.emplace_back(std::move(shw));
    }

    return shwOut;
  }

}