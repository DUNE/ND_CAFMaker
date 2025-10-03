#include "NDLArTMSMatchRecoFiller.h"

namespace cafmaker
{
  NDLArTMSMatchRecoFiller::NDLArTMSMatchRecoFiller()
    : IRecoBranchFiller("LArTMSMatcher")
  {
    // nothing to do
    SetConfigured(true);
  }

  // We need a comparator to go through the std::map and there's no default comparator for SRNDLar.
  struct LarIdComparator
  {
    bool operator()(const caf::SRNDLArID &lhs, const caf::SRNDLArID &rhs) const
    {
      // Define your custom comparison logic here
      // In this example, compare based on 'id' only
      if (lhs.ixn < rhs.ixn)
        return true;
      else
      {
        if (lhs.ixn == rhs.ixn)
          return lhs.idx < rhs.idx;
      }
      return false;
    }
  };

  double NDLArTMSMatchRecoFiller::Angular_Match(double xdir1_tms, double ydir1_tms, double zdir1_tms, double xdir2_lar, double ydir2_lar, double zdir2_lar) const
  {
      double dir_dot = (xdir1_tms * xdir2_lar) + (ydir1_tms * ydir2_lar) + (zdir1_tms * zdir2_lar);
      double abs_tms = sqrt((xdir1_tms * xdir1_tms) + (ydir1_tms * ydir1_tms) + (zdir1_tms * zdir1_tms));
      double abs_lar = sqrt((xdir2_lar * xdir2_lar) + (ydir2_lar * ydir2_lar) + (zdir2_lar * zdir2_lar));

      double delta_theta = acos(dir_dot / (abs_tms * abs_lar)) * (180/3.14); 

      return delta_theta;
  }

  bool NDLArTMSMatchRecoFiller::Spatial_Match(caf::SRTrack track_lar, caf::SRTrack track_tms, double &xthreshold, double &ythreshold, double &theta_xthreshold, double &theta_ythreshold, double &residual, double &ang_residual) const
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
    double xdir1_tms = track_tms.dir.x; //directional cosines
    double ydir1_tms = track_tms.dir.y;
    double zdir1_tms = track_tms.dir.z;

    //calculate the angle between track start and end directions in both the xz- and yz-planes
    double theta_x = Angular_Match(xdir1_tms, 0, zdir1_tms, xdir2_lar, 0, zdir2_lar);
    double theta_y = Angular_Match(0,ydir1_tms, zdir1_tms, 0, ydir2_lar, zdir2_lar);

    ang_residual = Angular_Match(xdir1_tms, ydir1_tms, zdir1_tms, xdir2_lar, ydir2_lar, zdir2_lar);

    // Propagate the LAr track into the TMS, by drawing a straight line into the TMS using the 
    // end direction in the LAr. Propagate the track to the same z as the TMS is, and 
    // compare the predicted vs actual position in x and y
    double x_residual = x2_lar + (xdir2_lar / zdir2_lar) * (z1_tms - z2_lar) - x1_tms;
    double y_residual = y2_lar + (ydir2_lar / zdir2_lar) * (z1_tms - z2_lar) - y1_tms;

    residual = sqrt(pow(x_residual, 2) + pow(y_residual, 2));

    if(abs(x_residual) < xthreshold && abs(y_residual) > ythreshold && abs(theta_x) < theta_xthreshold && abs(theta_y) < theta_ythreshold){
      return true;
    }
    else{
      return false;
    }
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
    double theta_xthreshold = 5; //degrees
    double theta_ythreshold = 10; //degrees

    std::map<caf::SRNDLArID, caf::SRNDTrackAssn, LarIdComparator> mult_map_spine, mult_map_pandora; // Lar to handle track matching multiplicity

    for(unsigned int ixn_tms = 0; ixn_tms < sr.nd.tms.nixn; ixn_tms++){
      caf::SRTMSInt tms_int = sr.nd.tms.ixn[ixn_tms];
      unsigned int n_tms_tracks = tms_int.ntracks;

      for(unsigned int itms = 0; itms < n_tms_tracks; ++itms){

        for(unsigned int ixn_dlp = 0; ixn_dlp < sr.nd.lar.ndlp; ixn_dlp++){
          caf::SRNDLArInt dlp = sr.nd.lar.dlp[ixn_dlp];
          for (unsigned int ilar = 0; ilar < dlp.ntracks; ++ilar){
            dlp.tracks[ilar];
            tms_int.tracks[itms];
            double residual = 0;
            double ang_residual = 0;

            if(Spatial_Match(dlp.tracks[ilar], tms_int.tracks[itms], xthreshold, ythreshold, theta_xthreshold, theta_ythreshold, residual, ang_residual))
            {
              caf::SRNDLArID larid;
              larid.ixn = ixn_dlp;
              larid.idx = ilar;
              larid.reco = caf::kDeepLearnPhys;
              caf::SRTMSID tmsid;
              tmsid.ixn = ixn_tms;
              tmsid.idx = itms;

              caf::SRNDTrackAssn match;
              match.larid = larid;
              match.tmsid = tmsid;
              match.transdispl = residual;
              match.angdispl = ang_residual;

              if (isnan(mult_map_spine[larid].transdispl))
                mult_map_spine[larid] = match;
              else
              {
                if (mult_map_spine[larid].transdispl > match.transdispl)
                  mult_map_spine[larid] = match;
              }
            }
          }
        }

        for(unsigned int ixn_pandora = 0; ixn_pandora < sr.nd.lar.npandora; ixn_pandora++){
          caf::SRNDLArInt pandora = sr.nd.lar.pandora[ixn_pandora];
          for (unsigned int ilar = 0; ilar < pandora.ntracks; ++ilar){
            pandora.tracks[ilar];
            tms_int.tracks[itms];
            double residual = 0;
            double ang_residual = 0;

            if(Spatial_Match(pandora.tracks[ilar], tms_int.tracks[itms], xthreshold, ythreshold, theta_xthreshold, theta_ythreshold, residual, ang_residual))
            {
              caf::SRNDLArID larid;
              larid.ixn = ixn_pandora;
              larid.idx = ilar;
              larid.reco = caf::kPandoraNDLAr;
              caf::SRTMSID tmsid;
              tmsid.ixn = ixn_tms;
              tmsid.idx = itms;

              caf::SRNDTrackAssn match;
              match.larid = larid;
              match.tmsid = tmsid;
              match.transdispl = residual;
              match.angdispl = ang_residual;

              if (isnan(mult_map_pandora[larid].transdispl))
                mult_map_pandora[larid] = match;
              else
              {
                if (mult_map_pandora[larid].transdispl > match.transdispl)
                  mult_map_pandora[larid] = match;
              }
            }
          }
        }
      }
    }

    for (auto m : mult_map_spine)
    {
      sr.nd.trkmatch.extrap.emplace_back(m.second);
      sr.nd.trkmatch.nextrap += 1;
    }

    for (auto m : mult_map_pandora)
    {
      sr.nd.trkmatch.extrap.emplace_back(m.second);
      sr.nd.trkmatch.nextrap += 1;
    }
  }

  std::deque<Trigger> NDLArTMSMatchRecoFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    return std::deque<Trigger>();
  }

}
