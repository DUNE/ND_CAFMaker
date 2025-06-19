#include "NDLArMINERvAMatchRecoFiller.h"

namespace cafmaker
{
  NDLArMINERvAMatchRecoFiller::NDLArMINERvAMatchRecoFiller(double _z_extr, double _d_x, double _d_y, double _d_thetax, double d_theta_y)
      : IRecoBranchFiller("LArMINERvAMatcher")
  {
    // setup matching criteria
    z_extr = _z_extr;
    d_x = _d_x;
    d_y = _d_y;
    d_thetax = _d_thetax;
    d_thetay = d_theta_y;
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

  bool NDLArMINERvAMatchRecoFiller::Passes_cut(caf::SRTrack track_minerva, caf::SRTrack track_Lar, double &costheta, double &residual) const
  {
    double x1_minerva = track_minerva.start.x;
    double x2_minerva = track_minerva.end.x;
    double y1_minerva = track_minerva.start.y;
    double y2_minerva = track_minerva.end.y;
    double z1_minerva = track_minerva.start.z;
    double z2_minerva = track_minerva.end.z;

    double x1_lar = track_Lar.start.x;
    double x2_lar = track_Lar.end.x;
    double y1_lar = track_Lar.start.y;
    double y2_lar = track_Lar.end.y;
    double z1_lar = track_Lar.start.z;
    double z2_lar = track_Lar.end.z;

    /*
    The experimental setup: Liquid Argon Detector is placed betbeen two MINERvA planes.
    To define matching criteria it is needed to find angles between LAr and MINERvA tracks.
    For LAr detector resolution is diffeent in X and Y direction, therefore it is needed to find angles between tracks
    as finction of the angle in X direction and as the function of an angle in Y direction. Distances between tracks
    will be calculated as distancec between extrapolated points - points of intersection of LAr and MINERvA tracls with
    the plane (parallel to plane XY) of LAr detector.
    */

    double tg_theta_mn_x = (x2_minerva - x1_minerva) / (z2_minerva - z1_minerva); // tangent of an angle between minerva track and X-axis
    double tg_theta_mn_y = (y2_minerva - y1_minerva) / (z2_minerva - z1_minerva); // tangent of an angle between minerva track and Y-axis
    double theta_mn_x = atan(tg_theta_mn_x);                                      // angle between minerva track and X-axis
    double theta_mn_y = atan(tg_theta_mn_y);                                      // angle between minerva track and Y-axis

    double tg_theta_nd_x = (x2_lar - x1_lar) / (z2_lar - z1_lar); // tangent of the angle between LAr track and X-axis
    double tg_theta_nd_y = (y2_lar - y1_lar) / (z2_lar - z1_lar); // tangent of the angle between LAr track and Y-axis
    double theta_nd_x = atan(tg_theta_nd_x);                      // angle between LAr track and X-axis
    double theta_nd_y = atan(tg_theta_nd_y);                      // angle between LAr track and Y-axis

    double delta_theta_x = theta_mn_y - theta_nd_y;
    double delta_theta_y = theta_mn_x - theta_nd_x;

    // Extrapolating Both tracks to the same point z = zextr (here it's the front of Lar)
    double t_mn = (z_extr - z1_minerva) / (z2_minerva - z1_minerva);
    double x_mn = t_mn * (x2_minerva - x1_minerva) + x1_minerva; // X-coordinate of extrapolated point of LAr track
    double y_mn = t_mn * (y2_minerva - y1_minerva) + y1_minerva; // Y-coordinate of extrapolated point of LAr track

    double t_nd = (z_extr - z1_lar) / (z2_lar - z1_lar); // parametr of the equation of the line (LAr track)
    double x_nd = t_nd * (x2_lar - x1_lar) + x1_lar;     // X-coordinate of extrapolated point of LAr track
    double y_nd = t_nd * (y2_lar - y1_lar) + y1_lar;     // Y-coordinate of extrapolated point of LAr track

    double dist_x = (x_mn - x_nd); // distance between X-coordinates of extrapolated points of minerva and LAr tracks
    double dist_y = (y_mn - y_nd); // distance between Y-coordinates of extrapolated points of minerva and LAr tracks

    residual = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
    costheta = ((x2_minerva - x1_minerva) * (x2_lar - x1_lar) +
                (y2_minerva - y1_minerva) * (y2_lar - y1_lar) +
                (z2_minerva - z1_minerva) * (z2_lar - z1_lar)) /
               (track_Lar.len_cm * track_minerva.len_cm); // angle between minerva and Lar tracks

    return (abs(delta_theta_x) < d_thetax && abs(delta_theta_y) < d_thetay && abs(dist_x) < d_y && abs(dist_y) < d_x);
  }

  void NDLArMINERvAMatchRecoFiller::_FillRecoBranches(const Trigger &trigger,
                                                      caf::StandardRecord &sr,
                                                      const cafmaker::Params &par,
                                                      const TruthMatcher *truthMatcher) const
  {
    // match tracks using the info that should have been filled by the ND-LAr and MINERvA reco filled

    std::map<caf::SRNDLArID, caf::SRNDTrackAssn, LarIdComparator> mult_map_spine, mult_map_pandora; // Lar to handle track matching multiplicity

    for (unsigned int ixn_minerva = 0; ixn_minerva < sr.nd.minerva.nixn; ixn_minerva++)
    {

      caf::SRMINERvAInt Mnv_int = sr.nd.minerva.ixn[ixn_minerva];
      unsigned int n_minerva_tracks = Mnv_int.ntracks; // # tracks in minerva in evtIdx

      for (unsigned int iminerva = 0; iminerva < n_minerva_tracks; ++iminerva)
      {

        for (unsigned int ixn_dlp = 0; ixn_dlp < sr.nd.lar.ndlp; ixn_dlp++) // SPINE Tracks
        {
          caf::SRNDLArInt dlp = sr.nd.lar.dlp[ixn_dlp];
          for (unsigned int ilar = 0; ilar < dlp.ntracks; ++ilar)
          {
            dlp.tracks[ilar];
            Mnv_int.tracks[iminerva];
            double residual = 0;
            double costheta = 0;
            if (Passes_cut(Mnv_int.tracks[iminerva], dlp.tracks[ilar], costheta, residual))
            {
              caf::SRMINERvAID mnvid;
              mnvid.ixn = ixn_minerva;
              mnvid.idx = iminerva;
              caf::SRNDLArID larid;
              larid.ixn = ixn_dlp;
              larid.idx = ilar;
              larid.reco = caf::kDeepLearnPhys;

              // Make the match

              caf::SRNDTrackAssn match;
              match.larid = larid;
              match.minervaid = mnvid;
              match.transdispl = residual;
              match.angdispl = costheta;

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

        for (unsigned int ixn_pandora = 0; ixn_pandora < sr.nd.lar.npandora; ixn_pandora++) // Pandora Tracks
        {
          caf::SRNDLArInt pandora = sr.nd.lar.pandora[ixn_pandora];
          for (unsigned int ilar = 0; ilar < pandora.ntracks; ++ilar)
          {
            pandora.tracks[ilar];
            Mnv_int.tracks[iminerva];
            double residual = 0;
            double costheta = 0;
            if (Passes_cut(Mnv_int.tracks[iminerva], pandora.tracks[ilar], costheta, residual))
            {
              caf::SRMINERvAID mnvid;
              mnvid.ixn = ixn_minerva;
              mnvid.idx = iminerva;
              caf::SRNDLArID larid;
              larid.ixn = ixn_pandora;
              larid.idx = ilar;
              larid.reco = caf::kPandoraNDLAr;

              // Make the match

              caf::SRNDTrackAssn match;
              match.larid = larid;
              match.minervaid = mnvid;
              match.transdispl = residual;
              match.angdispl = costheta;

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

  std::deque<Trigger> NDLArMINERvAMatchRecoFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    return std::deque<Trigger>();
  }
}
