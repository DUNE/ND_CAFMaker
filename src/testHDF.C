#include <cassert>
#include <iostream>

#include "reco/NDLArSummaryH5DatasetReader.h"


int main( int argc, char const *argv[] )
{
  // this test macro only accepts two arguments:
  //   * the name of the hdf5 file to read in
  //   * the name of the HDF dataset inside it
  assert(argc == 3);

  cafmaker::NDLArSummaryH5DatasetReader h5reader(argv[1], argv[2], "Event");

  std::set<std::size_t> evs = h5reader.Events();
  std::cout << "Dataset '" << argv[2] << "' within HDF5 file '" << argv[1] << "' contains " << evs.size() << " events"
            << " ranging from " << *evs.begin() << " to " << *evs.rbegin() << std::endl;
  //for (const auto & ev : h5reader.Events())
  //  std::cout << "  " << ev;
  //std::cout << std::endl;

  for (std::size_t evtIdx = 1; evtIdx < 50; evtIdx++)
  {
    std::vector<cafmaker::Track> tracks = h5reader.EventTracks(evtIdx);
    if (tracks.empty())
      continue;

    std::cout << "Tracks for event " << evtIdx << ":" << std::endl;
    for (const auto & trk : tracks)
    {
      std::cout << "  start = " << trk.start
                << ",  end = " << trk.end
                << ",  end dir = " << trk.end_dir
                << std::endl;
    }
  }
}