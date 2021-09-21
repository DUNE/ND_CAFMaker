#include <cassert>
#include <iostream>

#include "reco/NDLArProductFiller.h"
#include "reco/NDLArSummaryH5DatasetReader.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"


int main( int argc, char const *argv[] )
{
  // this test macro only accepts two arguments:
  //   * the name of the hdf5 file to read in
  //   * the name of the HDF dataset inside it
  assert(argc == 3);

  cafmaker::NDLArProductFiller<caf::SRTrack> trackFiller(argv[1], argv[2], "Event", "column_names");

  std::set<std::size_t> evs = trackFiller.DatasetReader().Events();
  std::cout << "Dataset '" << argv[2] << "' within HDF5 file '" << argv[1] << "' contains " << evs.size() << " events"
            << " ranging from " << *evs.begin() << " to " << *evs.rbegin() << std::endl;
  //for (const auto & ev : h5reader.Events())
  //  std::cout << "  " << ev;
  //std::cout << std::endl;

  for (std::size_t evtIdx = 0; evtIdx < std::min(50ul, *std::max_element(evs.begin(), evs.end())); evtIdx++)
  {
    caf::StandardRecord sr;

    trackFiller.FillSR(sr, evtIdx);
    if (sr.ndlar.ntracks < 1)
      continue;

    std::cout << "Tracks for event " << evtIdx << ":" << std::endl;
    for (const auto & trk : sr.ndlar.tracks)
    {
      std::cout << "  start = " << trk.start
                << ",  end = " << trk.end
                << ",  end dir = " << trk.end_dir
                << std::endl;
    }
  }
}