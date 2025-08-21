#include "NDLArTMSMatchRecoFiller.h"

namespace cafmaker
{
  NDLArTMSMatchRecoFiller::NDLArTMSMatchRecoFiller()
    : IRecoBranchFiller("LArTMSMatcher")
  {
    // nothing to do
    SetConfigured(true);
  }

  bool NDLArTMSMatchRecoFiller::Spatial_Match(caf::SRTrack track_lar, caf::SRTrack track_tms, double &xthreshold, double &ythreshold) const
  {
        double z2_lar = track_lar.end.z;     //cm

        // Check LAr track exits through the end; if not move to next track
        if (z2_lar < 800){
          return false;
        }

        // Now get the other variables
        double x2_lar = track_lar.end.x;     //cm
        double y2_lar = track_lar.end.y;     //cm
        double xdir2_lar = track_lar.enddir.x; //dir cosines
        double ydir2_lar = track_lar.enddir.y; //dir cosines
        double zdir2_lar= track_lar.enddir.z; //dir cosines

        double z1_tms = track_tms.start.z;   //cm

        // Check that this TMS track starts in first few planes; if not move to next track
        if (z1_tms > 1200){
          return false;
        }

        double x1_tms = track_tms.start.x; //cm
        double y1_tms = track_tms.start.y; //cm


        // Propagate the LAr track into the TMS, by drawing a straight line into the TMS using the 
        // end direction in the LAr. Propagate the track to the same z as the TMS is, and 
        // compare the predicted vs actual position in x and y
        double x_residual = x2_lar + (xdir2_lar / zdir2_lar) * (z1_tms - z2_lar) - x1_tms;
        double y_residual = y2_lar + (ydir2_lar / zdir2_lar) * (z1_tms - z2_lar) - y1_tms;

        return(abs(x_residual) < xthreshold && abs(y_residual) > ythreshold);
  }

  void NDLArTMSMatchRecoFiller::_FillRecoBranches(const Trigger &trigger,
                                             caf::StandardRecord &sr,
                                             const cafmaker::Params &par,
                                             const TruthMatcher *truthMatcher) const
  {
    // match tracks using the info that should have been filled by the ND-LAr and TMS reco fillers
    //sr.nd.ntrkmatch = 0; //outdated?
    
    //todo: pick better thresholds
    //thresholds in x and y that define a match
    double xthreshold = 7; //cm
    double ythreshold = 20; //cm

    for(unsigned int idlp = 0; idlp < sr.nd.lar.ndlp; ++idlp ) {
      caf::SRNDLArInt dlp = sr.nd.lar.dlp[idlp];

      for(unsigned int ilar = 0; ilar < dlp.ntracks; ++ilar){
        dlp.tracks[ilar];
        
        // Loop over the TMS tracks
        for(unsigned int ixn_tms = 0; ixn_tms < sr.nd.tms.nixn; ++ixn_tms ) {
          
          caf::SRTMSInt tms_int = sr.nd.tms.ixn[ixn_tms];
          unsigned int n_tms_tracks = tms_int.ntracks;

          for(unsigned int itms = 0; itms < n_tms_tracks; itms++){

            // Check matching
            if (Spatial_Match(dlp.tracks[ilar], tms_int.tracks[itms], xthreshold, ythreshold)) {
              // Make the match
              caf::SRNDLArID larid;
              larid.ixn = idlp;
              larid.idx = ilar;
              caf::SRTMSID tmsid;
              tmsid.ixn = ixn_tms;
              tmsid.idx = itms;

              caf::SRNDTrackAssn match;
              match.larid = larid;
              match.tmsid = tmsid;
    //          match.transdispl = residual;
    //          match.angdispl = costheta;
              sr.nd.trkmatch.extrap.emplace_back(match);
              //sr.nd.ntrkmatch += 1; //outdated?
    //          // Can also have multiple track matches, so don't break at the end of the match
    //          // Should pick the best track match in selection code
            }
          }
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
