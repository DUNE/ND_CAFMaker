#include "NDLArTMSMatchRecoFiller.h"

namespace cafmaker
{
  NDLArTMSMatchRecoFiller::NDLArTMSMatchRecoFiller()
  {
    // nothing to do
    SetConfigured(true);
    name = "LArTMSMatcher";
  }

  void
  NDLArTMSMatchRecoFiller::_FillRecoBranches(std::size_t evtIdx,
                                             caf::StandardRecord &sr,
                                             const cafmaker::Params &par,
                                             const TruthMatcher *truthMatcher) const
  {
//    // match tracks using the info that should have been filled by the ND-LAr and TMS reco fillers
//    unsigned int n_TMS_tracks = sr.nd.tms.ntracks;  //# tracks in TMS in evtIdx
//    unsigned int n_LAr_tracks = sr.nd.lar.ntracks;  //# tracks in LAr in evtIdx
//
//    sr.nd.ntrkmatch = 0;
//
//    for(unsigned int ilar = 0; ilar < n_LAr_tracks; ++ilar ) {
//      double z2_lar = sr.nd.lar.tracks[ilar].end.z;     //cm
//
//      // Check LAr track exits through the end; if not move to next track
//      if (z2_lar < 800) continue;
//
//      // Now get the other variables
//      double x1_lar = sr.nd.lar.tracks[ilar].start.x;   //cm
//      double z1_lar = sr.nd.lar.tracks[ilar].start.z;   //cm
//      double x2_lar = sr.nd.lar.tracks[ilar].end.x;     //cm
//
//      // The slope of the track in LAr
//      // Unfortunately the saved unit vectors in the LAr reco are corrupted
//      // This should be replaced with the actual saved unit vectors
//      // Hence we calculate the slope from the start and end position of the track in LAr instead
//      double slope_lar = (x2_lar-x1_lar)/(z2_lar-z1_lar);
//
//      // Loop over the TMS tracks
//      for(unsigned int itms = 0; itms < n_TMS_tracks; ++itms ) {
//        double z1_tms = sr.nd.tms.tracks[itms].start.z;   //cm
//
//        // Check that this TMS track starts in first few planes; if not move to next track
//        if (z1_tms > 1200) continue;
//
//        double x1_tms = sr.nd.tms.tracks[itms].start.x;   //cm
//
//        // Use the most upstream slope of the TMS track 
//        double slope_x_tms = sr.nd.tms.tracks[itms].dir.x;
//        double slope_z_tms = sr.nd.tms.tracks[itms].dir.z;
//
//        // Propagate the LAr track into the TMS, by drawing a straight line into the TMS using the 
//        // slope in the LAr. Propagate the track to the same z as the TMS is, and 
//        // compare the predicted vs actual position in x
//        double residual = x1_lar + slope_lar * (z1_tms - z1_lar) - x1_tms;
//
//        // Calculate the unit vectors
//        double len = sqrt((x2_lar-x1_lar)*(x2_lar-x1_lar) + (z2_lar-z1_lar)*(z2_lar-z1_lar));
//        double slope_x_lar = (x2_lar-x1_lar)/len;
//        double slope_z_lar = (z2_lar-z1_lar)/len;
//        double costheta = slope_x_lar * slope_x_tms + slope_z_lar * slope_z_tms;
//
//        // Check matching; 40cm and 0.95 costheta between the unit vectors in LAr and TMS
//        if (abs(residual) < 40. && fabs(costheta) > 0.95) {
//          // Make the match
//          caf::SRNDTrackAssn match;
//          match.larid = ilar;
//          match.tmsid = itms;
//          match.transdispl = residual;
//          match.angdispl = costheta;
//          sr.nd.trkmatch.emplace_back(match);
//          sr.nd.ntrkmatch += 1;
//          // Can also have multiple track matches, so don't break at the end of the match
//          // Should pick the best track match in selection code
//        }
//      }
//    }
  }
}
