#include "NDLArMINERvAMatchRecoFiller.h"

namespace cafmaker
{
  NDLArMINERvAMatchRecoFiller::NDLArMINERvAMatchRecoFiller()
      : IRecoBranchFiller("LArMINERvAMatcher")
  {
    // nothing to do
    SetConfigured(true);
  }

  void NDLArMINERvAMatchRecoFiller::_FillRecoBranches(const Trigger &trigger,
                                                      caf::StandardRecord &sr,
                                                      const cafmaker::Params &par,
                                                      const TruthMatcher *truthMatcher) const

  {
    // match tracks using the info that should have been filled by the ND-LAr and MINERvA reco filled

    double z_extr = -70; // Extrapolated position compariton. Here it's the front of the Lar modules.

    double d_x = 17 ; // Maximum residual in x coordinate [cm];
    double d_y = 19 ; // Maximum residual in y coordinate [cm];
    double d_thetax = .08 ; // Maximum Angle difference wrt to x axis [rad]; 
    double d_thetay = .09 ; // Maximum Angle difference wrt to y axis [rad];
    
    std::map<caf::SRNDLArID, caf::SRNDTrackAssn> mult_map; //Lar to handle track matching multiplicity

    for (unsigned int ixn_dlp = 0; ixn_dlp < sr.nd.lar.ndlp; ixn_dlp++) //SPINE Tracks
    {
      caf::SRNDLArInt dlp = sr.nd.lar.dlp[ixn_dlp];
      for (unsigned int ilar = 0; ilar < dlp.ntracks; ++ilar)
      {
        // starting and ending positions of LAr tracks
        double x1_lar = dlp.tracks[ilar].start.x;
        double x2_lar = dlp.tracks[ilar].end.x;
        double y1_lar = dlp.tracks[ilar].start.y;
        double y2_lar = dlp.tracks[ilar].end.y;
        double z1_lar = dlp.tracks[ilar].start.z;
        double z2_lar = dlp.tracks[ilar].end.z;

        // Loop over the minerva tracks

        for (unsigned int ixn_minerva = 0; ixn_minerva < sr.nd.minerva.nixn; ixn_minerva++)
        {

          caf::SRMINERvAInt Mnv_int = sr.nd.minerva.ixn[ixn_minerva];
          unsigned int n_minerva_tracks = Mnv_int.ntracks; // # tracks in minerva in evtIdx

          for (unsigned int iminerva = 0; iminerva < n_minerva_tracks; ++iminerva)
          {

            // starting and ending positions of minerva tracks
            double x1_minerva = Mnv_int.tracks[iminerva].start.x;
            double x2_minerva = Mnv_int.tracks[iminerva].end.x;
            double y1_minerva = Mnv_int.tracks[iminerva].start.y;
            double y2_minerva = Mnv_int.tracks[iminerva].end.y;
            double z1_minerva = Mnv_int.tracks[iminerva].start.z;
            double z2_minerva = Mnv_int.tracks[iminerva].end.z;

            /*
            The experimental setup: Liquid Argon Detector is placed betbeen two MINERvA planes.
            To define matching criteria it is needed to find angles between LAr and MINERvA tracks.
            For LAr detector resolution is diffeent in X and Y direction, therefore it is needed to find angles between tracks
            as finction of the angle in X direction and as the function of an angle in Y direction. Distances between tracks
            will be calculated as distancec between extrapolated points - points of intersection of LAr and MINERvA tracls with
            the plane (parallel to plane XY) of LAr detector.
            */

            double tg_theta_mn_x; // tangent of an angle between minerva track and X-axis
            tg_theta_mn_x = (x2_minerva - x1_minerva) / (z2_minerva - z1_minerva);
            double tg_theta_mn_y; // tangent of an angle between minerva track and Y-axis
            tg_theta_mn_y = (y2_minerva - y1_minerva) / (z2_minerva - z1_minerva);
            double theta_mn_x; // angle between minerva track and X-axis
            double theta_mn_y; // angle between minerva track and Y-axis
            theta_mn_x = atan(tg_theta_mn_x);
            theta_mn_y = atan(tg_theta_mn_y);

            double tg_theta_nd_x; // tangent of the angle between LAr track and X-axis
            double tg_theta_nd_y; // tangent of the angle between LAr track and Y-axis
            tg_theta_nd_x = (x2_lar - x1_lar) / (z2_lar - z1_lar);
            tg_theta_nd_y = (y2_lar - y1_lar) / (z2_lar - z1_lar);
            double theta_nd_x; // angle between LAr track and X-axis
            double theta_nd_y; // angle between LAr track and Y-axis
            theta_nd_x = atan(tg_theta_nd_x);
            theta_nd_y = atan(tg_theta_nd_y);

            double delta_theta_x;
            double delta_theta_y;
            delta_theta_y = theta_mn_y - theta_nd_y;
            delta_theta_x = theta_mn_x - theta_nd_x;
            double t_mn; // parametre of the equation of the line (minerva track)

            //Extrapolating Both tracks to the same point z = zextr (here it's the front of Lar)
            t_mn = (z_extr - z1_minerva) / (z2_minerva - z1_minerva);

            double x_mn; // X-coordinate of extrapolated point of LAr track
            double y_mn; // Y-coordinate of extrapolated point of LAr track
            // double z_mn;
            x_mn = t_mn * (x2_minerva - x1_minerva) + x1_minerva;
            y_mn = t_mn * (y2_minerva - y1_minerva) + y1_minerva;
            // z_mn = t_mn * (z2_minerva - z1_minerva) + z1_minerva;

            double t_nd; // parametr of the equation of the line (LAr track)
            t_nd = (z_extr - z1_lar) / (z2_lar - z1_lar);

            double x_nd; // X-coordinate of extrapolated point of LAr track
            double y_nd; // Y-coordinate of extrapolated point of LAr track
            // double z_nd;
            x_nd = t_nd * (x2_lar - x1_lar) + x1_lar;
            y_nd = t_nd * (y2_lar - y1_lar) + y1_lar;
            // z_nd = t_nd * (z2_lar - z1_lar) + z1_lar;

            double dist_x; // distance between X-coordinates of extrapolated points of minerva and LAr tracks
            double dist_y; // distance between Y-coordinates of extrapolated points of minerva and LAr tracks
            dist_x = (x_mn - x_nd);
            dist_y = (y_mn - y_nd);

            double len_mn = sqrt((x2_minerva - x1_minerva) * (x2_minerva - x1_minerva) + (z2_minerva - z1_minerva) * (z2_minerva - z1_minerva) + (y2_minerva - y1_minerva) * (y2_minerva - y1_minerva)); // magnitude of minerva track
            double len_nd = sqrt((x2_lar - x1_lar) * (x2_lar - x1_lar) + (z2_lar - z1_lar) * (z2_lar - z1_lar) + (y2_lar - y1_lar) * (y2_lar - y1_lar));                                                 // magnitude of LAr track

            double residual = 0;
            residual = sqrt(pow(dist_x, 2) + pow(dist_y, 2));

            double costheta = 0;
            costheta = ((x2_minerva - x1_minerva) * (x2_lar - x1_lar) + (y2_minerva - y1_minerva) * (y2_lar - y1_lar) + (z2_minerva - z1_minerva) * (z2_lar - z1_lar)) / (len_nd * len_mn); // angle between minerva and Lar tracks

            if (abs(delta_theta_x) < d_thetax && abs(delta_theta_y) < d_thetay && abs(dist_x) < d_y && abs(dist_y) < d_x )
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

              if (isnan(mult_map[larid].transdispl)) mult_map[larid] = match;
              else 
              {
                if (mult_map[larid].transdispl >  match.transdispl) mult_map[larid] = match;
              }
              mult_map[larid].push_back(match); // Add to the track ixn_dlp, ilar 1 Minerva match to handle multiplicity.
            }
          }
        }
      }
    }
    for (auto m : mult)
    {
      sr.nd.trkmatch.extrap.emplace_back(m.second);
      sr.nd.trkmatch.nextrap += 1;
    }
  }
  std::deque<Trigger> NDLArMINERvAMatchRecoFiller::GetTriggers(int triggerType) const
  {
    return std::deque<Trigger>();
  }

}
