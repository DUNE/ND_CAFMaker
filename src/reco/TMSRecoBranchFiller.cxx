#include "TMSRecoBranchFiller.h"

namespace cafmaker
{
  TMSRecoBranchFiller::TMSRecoBranchFiller(const std::string &tmsRecoFilename)
  {
    fTMSRecoFile = new TFile(tmsRecoFilename.c_str(), "READ");
    name = std::string("TMS");
    if (!fTMSRecoFile->IsZombie()) {
      SetConfigured(true);
      // Save pointer to input tree
      TMSRecoTree = dynamic_cast<TTree*>(fTMSRecoFile->Get("Line_Candidates"));
      TMSRecoTree->SetBranchAddress("nLines", &_nLines);
      TMSRecoTree->SetBranchAddress("nHitsInTrack", _nHitsInTrack);
      TMSRecoTree->SetBranchAddress("TrackLength", _TrackLength);
      TMSRecoTree->SetBranchAddress("TotalTrackEnergy", _TotalTrackEnergy);
      TMSRecoTree->SetBranchAddress("Occupancy", _Occupancy);

      TMSRecoTree->SetBranchAddress("DirectionX_Upstream", _DirectionX_Upstream);
      TMSRecoTree->SetBranchAddress("DirectionZ_Upstream", _DirectionZ_Upstream);

      TMSRecoTree->SetBranchAddress("DirectionX_Downstream", _DirectionX_Downstream);
      TMSRecoTree->SetBranchAddress("DirectionZ_Downstream", _DirectionZ_Downstream);

      TMSRecoTree->SetBranchAddress("TrackHitPos", _TrackHitPos);
    } else {
      fTMSRecoFile = NULL;
      TMSRecoTree = NULL;
    }
  }

  // ---------------------------------------------------------------------------

  // here we copy all the TMS reco into the SRTMS branch of the StandardRecord object.
  void TMSRecoBranchFiller::_FillRecoBranches(std::size_t evtIdx,
                                              caf::StandardRecord &sr,
                                              const cafmaker::Params &par) const
  {
    // Get nth entry from tree
    TMSRecoTree->GetEntry(evtIdx);

    // First set number of tracks
    sr.nd.tms.ntracks = _nLines;
    sr.nd.tms.tracks.resize(_nLines);

    // Fill in the track info 
    for (int i = 0; i < _nLines; ++i) {
      double prevz = -9E10;
      // Hit info
      // Loop might be unnecessary if you really want optimisation...
      for (int j = 0; j < _nHitsInTrack[i]; ++j) {
        double z = _TrackHitPos[i][j][0];
        //double x = _TrackHitPos[i][j][1];
        //std::cout << "line " << i << ", hit " << j << ", x, z = " << x << ", " << z << std::endl;

        // Should all be ordered in z, check this
        if (z < prevz) {
          std::cerr << "hits in z not ordered" << std::endl;
          throw;
        }
        prevz = z;
      }

      // Save first and last hit in track
      sr.nd.tms.tracks[i].start = caf::SRVector3D(_TrackHitPos[i][0][1], -999, _TrackHitPos[i][0][0]);
      sr.nd.tms.tracks[i].end   = caf::SRVector3D(_TrackHitPos[i][_nHitsInTrack[i]-1][1], -999, _TrackHitPos[i][_nHitsInTrack[i]-1][0]);

      // Track info
      sr.nd.tms.tracks[i].len_gcm3  = _TrackLength[i];
      sr.nd.tms.tracks[i].qual      = _Occupancy[i];
      sr.nd.tms.tracks[i].E         = _TotalTrackEnergy[i];

      // Get the directions
      sr.nd.tms.tracks[i].dir     = caf::SRVector3D(_DirectionX_Upstream[i], -999, _DirectionZ_Upstream[i]);
      sr.nd.tms.tracks[i].enddir  = caf::SRVector3D(_DirectionX_Downstream[i], -999, _DirectionZ_Downstream[i]);
    }

  }
} // end namespace
