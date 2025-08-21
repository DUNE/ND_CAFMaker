#include "NDLArTMSMatchRecoFiller.h"

namespace cafmaker
{
  NDLArTMSMatchRecoFiller::NDLArTMSMatchRecoFiller()
    : IRecoBranchFiller("LArTMSMatcher")
  {
    // nothing to do
    SetConfigured(true);
  }

  void
  NDLArTMSMatchRecoFiller::_FillRecoBranches(const Trigger &trigger,
                                             caf::StandardRecord &sr,
                                             const cafmaker::Params &par,
                                             const TruthMatcher *truthMatcher) const
  {
    // match tracks using the info that should have been filled by the ND-LAr and TMS reco fillers
    unsigned int n_TMS_tracks = sr.nd.tms.ntracks;  //# tracks in TMS in evtIdx
    unsigned int n_LAr_tracks = sr.nd.lar.ntracks;  //# tracks in LAr in evtIdx

    sr.nd.ntrkmatch = 0;
    
    //todo: pick better thresholds
    //thresholds in x and y that define a match
    xthreshold = 7; //cm
    ythreshold = 20; //cm

    for(unsigned int ilar = 0; ilar < n_LAr_tracks; ++ilar ) {
      double z2_lar = sr.nd.lar.tracks[ilar].end.z;     //cm

      // Check LAr track exits through the end; if not move to next track
      if (z2_lar < 800) continue;

      // Now get the other variables
      double x2_lar = sr.nd.lar.tracks[ilar].end.x;     //cm
      double y2_lar = sr.nd.lar.tracks[ilar].end.y;     //cm
      double xdir2_lar = sr.nd.lar.tracks[ilar]. enddir.x; //dir cosines
      double ydir2_lar = sr.nd.lar.tracks[ilar].enddir.y; //dir cosines
      double zdir2_lar= sr.nd.lar.tracks[ilar].enddir.z; //dir cosines

      // Loop over the TMS tracks
      for(unsigned int itms = 0; itms < n_TMS_tracks; ++itms ) {
        double z1_tms = sr.nd.tms.tracks[itms].start.z;   //cm

        // Check that this TMS track starts in first few planes; if not move to next track
        if (z1_tms > 1200) continue;

        double x1_tms = sr.nd.tms.tracks[itms].start.x; //cm
        double y1_tms = sr.nd.tracks[itms].start.y; //cm


        // Propagate the LAr track into the TMS, by drawing a straight line into the TMS using the 
        // end direction in the LAr. Propagate the track to the same z as the TMS is, and 
        // compare the predicted vs actual position in x and y
        double x_residual = x2_lar + (xdir2_lar / zdir2_lar) * (z1_tms - z2_lar) - x1_tms;
        double y_residual = y2_lar + (ydir2_lar / zdir2_lar) * (z1_tms - z2_lar) - y1_tms;

        // Old from when LAr directions were not working
//        double len = sqrt((x2_lar-x1_lar)*(x2_lar-x1_lar) + (z2_lar-z1_lar)*(z2_lar-z1_lar));
//        double slope_x_lar = (x2_lar-x1_lar)/len;
//        double slope_z_lar = (z2_lar-z1_lar)/len;
//        double costheta = slope_x_lar * slope_x_tms + slope_z_lar * slope_z_tms;

        // Check matching
        if (abs(x_residual) < xthreshold && abs(y_residual) > ythreshold) {
          // Make the match
          caf::SRNDTrackAssn match;
          match.larid = ilar;
          match.tmsid = itms;
//          match.transdispl = residual;
//          match.angdispl = costheta;
//          sr.nd.trkmatch.emplace_back(match);
          sr.nd.ntrkmatch += 1;
//          // Can also have multiple track matches, so don't break at the end of the match
//          // Should pick the best track match in selection code
        }
      }
    }
  }

  // todo: this is a placeholder
  std::deque<Trigger> NDLArTMSMatchRecoFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    return std::deque<Trigger>();
  }

}
