#include "NDLArTMSUniqueMatchRecoFiller.h"
#include <cmath>
#include "TRandom3.h"
#include <set>

#include "Math/Functor.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Rtypes.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDSymfwd.h"
#include "TMatrixDfwd.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TPolyLine3D.h"
#include "TSpline.h"
#include "TVectorDfwd.h"
#include "TVectorT.h"

namespace cafmaker
{
  bool Track_match_sorter(const caf::SRNDTrackAssn trackMatch1, const caf::SRNDTrackAssn trackMatch2) {
    double fScore1 = trackMatch1.matchScore;
    double fScore2 = trackMatch2.matchScore;

    if (fScore1 < fScore2) {
      return true;
    }
    else {
      return false;
    }
  }

  NDLArTMSUniqueMatchRecoFiller::NDLArTMSUniqueMatchRecoFiller(const double sigmaX, const double sigmaY, const double sigmaThX, const double sigmaThY, const bool useTime, const double meanT, const double sigmaT, const double fCut, const bool useSmearTime)
    : IRecoBranchFiller("LArTMSMatcher")
  {
    sigma_x = sigmaX;
    sigma_y = sigmaY;
    sigma_angle_x = sigmaThX;
    sigma_angle_y = sigmaThY;
    use_time = useTime;
    mean_t = meanT;
    sigma_t = sigmaT;
    f_cut = fCut;
    use_smear_time = useSmearTime;
    // nothing to do
    SetConfigured(true);
  }

  std::vector<float> NDLArTMSUniqueMatchRecoFiller::Project_track(const caf::SRTrack track, const bool forward) const
  {
    float x, y, z;

    float dir_x, dir_y, dir_z;
    
    float proj_z;
    float proj_x;
    float proj_y;

    if (forward) { // projects a LAr track forward to TMS
       x = track.end.x;
       y = track.end.y;
       z = track.end.z;

       dir_x = track.enddir.x;
       dir_y = track.enddir.y;
       dir_z = track.enddir.z;

       proj_z = tms_z_lim1 - z;
       proj_x = dir_x*proj_z/dir_z + x;
       proj_y = dir_y*proj_z/dir_z + y;

       std::cout << "    Calculated Track End Point Z " << z << std::endl;
       std::cout << "    Calculated TMS Start Z " << tms_z_lim1 << std::endl;
       std::cout << "    Calculated Track Z Distance " << tms_z_lim1 - z << std::endl;
       std::cout << "    Calculated Track End Dir Z " << dir_z << std::endl;
       std::cout << "    Calculated Z Distance / Dir Z " << (tms_z_lim1 - z)/dir_z << std::endl;
       std::cout << "    Calculated Track End Point X " << x << std::endl;
       std::cout << "    Calculated Track End Dir X " << dir_x << std::endl;
       std::cout << "    Calculated X Displacement " << dir_x*(tms_z_lim1 - z)/dir_z << std::endl;
       std::cout << "    Calculated X Projection " << x + dir_x*(tms_z_lim1 - z)/dir_z << std::endl;
       std::cout << "    Calculated Track End Point X " << x << std::endl;
       std::cout << "    Calculated Track End Dir X " << dir_x << std::endl;
       std::cout << "    Calculated X Displacement " << dir_x*(tms_z_lim1 - z)/dir_z << std::endl;
       std::cout << "    Calculated X Projection " << x + dir_x*(tms_z_lim1 - z)/dir_z << std::endl;
    }
    else { // projects a TMS track backward to LAr
       x = track.start.x;
       y = track.start.y;
       z = track.start.z;

       dir_x = track.dir.x;
       dir_y = track.dir.y;
       dir_z = track.dir.z;

       proj_z = z - lar_z_lim2;
       proj_x = -dir_x*proj_z/dir_z + x;
       proj_y = -dir_y*proj_z/dir_z + y;
    }
    std::vector<float> proj_point = {proj_x, proj_y, proj_z};
    return proj_point;
  }

  std::vector<float> NDLArTMSUniqueMatchRecoFiller::Angle_between_tracks(const caf::SRTrack tms_track, const caf::SRTrack lar_track) const
  {
      float tms_dir_x = tms_track.dir.x;
      float tms_dir_y = tms_track.dir.y;
      float tms_dir_z = tms_track.dir.z;

      float lar_dir_x = lar_track.enddir.x;
      float lar_dir_y = lar_track.enddir.y;
      float lar_dir_z = lar_track.enddir.z;

      float xz_dot_prod = tms_dir_x*lar_dir_x + tms_dir_z*lar_dir_z;
      if (xz_dot_prod != 0) {
        xz_dot_prod = xz_dot_prod/(sqrt(pow(tms_dir_x,2)+pow(tms_dir_z,2))*sqrt(pow(lar_dir_x,2)+pow(lar_dir_z,2)));
      }
      float yz_dot_prod = tms_dir_y*lar_dir_y + tms_dir_z*lar_dir_z;
      if (yz_dot_prod != 0) {
        yz_dot_prod = yz_dot_prod/(sqrt(pow(tms_dir_y,2)+pow(tms_dir_z,2))*sqrt(pow(lar_dir_y,2)+pow(lar_dir_z,2)));
      }
      float dot_prod = tms_dir_x*lar_dir_x + tms_dir_y*lar_dir_y + tms_dir_z*lar_dir_z;
      float angle_x = 180.0/TMath::Pi() * acos(xz_dot_prod);
      float angle_y = 180.0/TMath::Pi() * acos(yz_dot_prod);
      float angle_overall = 180.0/TMath::Pi() * acos(dot_prod);
      std::vector<float> angles = {angle_x,angle_y,angle_overall};
      return angles;
  }

  bool NDLArTMSUniqueMatchRecoFiller::Consider_TMS_track(const caf::SRTrack tms_track, const double tms_z_cutoff) const
  {
    float x_start = tms_track.start.x;
    float y_start = tms_track.start.y;
    float z_start = tms_track.start.z;

    if ((x_start > tms_x_lim1)&&(x_start < tms_x_lim2) &&
        (y_start > tms_y_lim1)&&(y_start < tms_y_lim2) &&
        (z_start > tms_z_lim1 -20)&&(z_start < tms_z_lim1 + tms_z_cutoff) && // checks track begins within fiducial volume and close enough to front
      
        (Project_track(tms_track,false)[0] > lar_x_lim1)&&(Project_track(tms_track,false)[0] < lar_x_lim2) &&
        (Project_track(tms_track,false)[1] > lar_y_lim1)&&(Project_track(tms_track,false)[1] < lar_y_lim2)) // checks that direction would have allowed it to originate from LAr
          {
            return true;
          } 
    else {
            if (x_start < tms_x_lim1){
              std::cout << "    X_START < TMS_X_LIM1" << std::endl;
              std::cout << "     x_start " << x_start << std::endl;
              std::cout << "     tms_x_lim1 " << tms_x_lim1 << std::endl;
            }
            if (x_start > tms_x_lim2){
              std::cout << "    X_START > TMS_X_LIM2" << std::endl;
              std::cout << "     x_start " << x_start << std::endl;
              std::cout << "     tms_x_lim2 " << tms_x_lim2 << std::endl;
            }
            if (y_start < tms_y_lim1){
              std::cout << "    Y_START < TMS_Y_LIM1" << std::endl;
              std::cout << "     y_start " << y_start << std::endl;
              std::cout << "     tms_y_lim1 " << tms_y_lim1 << std::endl;
            }
            if (y_start > tms_y_lim2){
              std::cout << "    Y_START > TMS_Y_LIM2" << std::endl;
              std::cout << "     y_start " << y_start << std::endl;
              std::cout << "     tms_y_lim2 " << tms_y_lim2 << std::endl;
            }
            if (z_start < tms_z_lim1-20){
              std::cout << "    Z_START < TMS_Z_LIM1-20" << std::endl;
              std::cout << "     z_start " << z_start << std::endl;
              std::cout << "     tms_z_lim1 " << tms_z_lim1 << std::endl;
            }
            if (z_start > tms_z_lim1 + tms_z_cutoff){
              std::cout << "    Z_START > TMS_Z_LIM1 + TMS_Z_CUTOFF" << std::endl;
              std::cout << "     z_start " << z_start << std::endl;
              std::cout << "     tms_z_lim1 " << tms_z_lim1 << std::endl;
              std::cout << "     tms_z_cutoff " << tms_z_cutoff << std::endl;
            }



            if (Project_track(tms_track,false)[0] < lar_x_lim1){
              std::cout << "    PROJECT_TRACK(TMS_TRACK,FALSE)[0] < LAR_X_LIM1" << std::endl;
              std::cout << "     Project_track(tms_track,false)[0] " << Project_track(tms_track,false)[0] << std::endl;
              std::cout << "     lar_x_lim1 " << lar_x_lim1 << std::endl;
            }
            if (Project_track(tms_track,false)[0] > lar_x_lim2){
              std::cout << "    PROJECT_TRACK(TMS_TRACK,FALSE)[0] > LAR_X_LIM2" << std::endl;
              std::cout << "     Project_track(tms_track,false)[0] " << Project_track(tms_track,false)[0] << std::endl;
              std::cout << "     lar_x_lim2 " << lar_x_lim2 << std::endl;
            } 
            if (Project_track(tms_track,false)[1] < lar_y_lim1){
              std::cout << "    PROJECT_TRACK(TMS_TRACK,FALSE)[1] < LAR_Y_LIM1" << std::endl;
              std::cout << "     Project_track(tms_track,false)[1] " << Project_track(tms_track,false)[1] << std::endl;
              std::cout << "     lar_y_lim1 " << lar_y_lim1 << std::endl;
            } 
            if (Project_track(tms_track,false)[1] > lar_y_lim2){
              std::cout << "     PROJECT_TRACK(TMS_TRACK,FALSE)[1] > LAR_Y_LIM2" << std::endl;
              std::cout << "      Project_track(tms_track,false)[1] " << Project_track(tms_track,false)[1] << std::endl;
              std::cout << "      lar_y_lim2 " << lar_y_lim2 << std::endl;
            } 

            return false;
    }
  }

  bool NDLArTMSUniqueMatchRecoFiller::Consider_LAr_track(const caf::SRTrack lar_track, const double lar_z_cutoff) const
  {
    float x_start = lar_track.start.x;
    float y_start = lar_track.start.y;
    float z_start = lar_track.start.z;

    float x_end = lar_track.end.x;
    float y_end = lar_track.end.y;
    float z_end = lar_track.end.z;

    if ((x_start > lar_x_lim1)&&(x_start < lar_x_lim2) &&
        (y_start > lar_y_lim1)&&(y_start < lar_y_lim2) &&
        (z_start > lar_z_lim1)&&(z_start < lar_z_lim2) && // checks track begins within fiducial volume

        (x_end > lar_x_lim1)&&(x_end < lar_x_lim2) &&
        (y_end > lar_y_lim1)&&(y_end < lar_y_lim2) &&
        (z_end > lar_z_lim2 - lar_z_cutoff)&&(z_end < lar_z_lim2) && // checks track ends close enough to back of LAr
      
        (Project_track(lar_track,true)[0] > tms_x_lim1)&&(Project_track(lar_track,true)[0] < tms_x_lim2) &&
        (Project_track(lar_track,true)[1] > tms_y_lim1)&&(Project_track(lar_track,true)[1] < tms_y_lim2)) // checks that direction would allow it to hit TMS
        {
          return true;
        } 
    else {
            if (x_start < lar_x_lim1){
              std::cout << "    X_START < LAr_X_LIM1" << std::endl;
              std::cout << "     x_start " << x_start << std::endl;
              std::cout << "     lar_x_lim1 " << lar_x_lim1 << std::endl;
            }
            if (x_start > lar_x_lim2){
              std::cout << "    X_START > LAr_X_LIM2" << std::endl;
              std::cout << "     x_start " << x_start << std::endl;
              std::cout << "     lar_x_lim2 " << lar_x_lim2 << std::endl;
            }
            if (y_start < lar_y_lim1){
              std::cout << "    Y_START < LAr_Y_LIM1" << std::endl;
              std::cout << "     y_start " << y_start << std::endl;
              std::cout << "     lar_y_lim1 " << lar_y_lim1 << std::endl;
            }
            if (y_start > lar_y_lim2){
              std::cout << "    Y_START > LAr_Y_LIM2" << std::endl;
              std::cout << "     y_start " << y_start << std::endl;
              std::cout << "     lar_y_lim2 " << lar_y_lim2 << std::endl;
            }
            if (z_start < lar_z_lim1){
              std::cout << "    Z_START < LAr_Z_LIM1" << std::endl;
              std::cout << "     z_start " << z_start << std::endl;
              std::cout << "     lar_z_lim1 " << lar_z_lim1 << std::endl;
            }
            if (z_start > lar_z_lim2){
              std::cout << "    Z_START > LAr_Z_LIM2" << std::endl;
              std::cout << "     z_start " << z_start << std::endl;
              std::cout << "     lar_z_lim2 " << lar_z_lim2 << std::endl;
            }

            if (x_end < lar_x_lim1){
              std::cout << "    X_END < LAr_X_LIM1" << std::endl;
              std::cout << "     x_end " << x_end << std::endl;
              std::cout << "     lar_x_lim1 " << lar_x_lim1 << std::endl;
            }
            if (x_end > lar_x_lim2){
              std::cout << "    X_END > LAr_X_LIM2" << std::endl;
              std::cout << "     x_end " << x_end << std::endl;
              std::cout << "     lar_x_lim2 " << lar_x_lim2 << std::endl;
            }
            if (y_end < lar_y_lim1){
              std::cout << "    Y_END < LAr_Y_LIM1" << std::endl;
              std::cout << "     y_end " << y_end << std::endl;
              std::cout << "     lar_y_lim1 " << lar_y_lim1 << std::endl;
            }
            if (y_end > lar_y_lim2){
              std::cout << "    Y_END > LAr_Y_LIM2" << std::endl;
              std::cout << "     y_end " << y_end << std::endl;
              std::cout << "     lar_y_lim2 " << lar_y_lim2 << std::endl;
            }
            if (z_end < lar_z_lim2 - lar_z_cutoff){
              std::cout << "    Z_END < LAr_Z_LIM2 - LAr_Z_CUTOFF" << std::endl;
              std::cout << "     z_end " << z_end << std::endl;
              std::cout << "     lar_z_lim2 " << lar_z_lim2 << std::endl;
              std::cout << "     lar_z_cutoff " << lar_z_cutoff << std::endl;
            }
            if (z_end > lar_z_lim2){
              std::cout << "    Z_END > LAr_Z_LIM2" << std::endl;
              std::cout << "     z_end " << z_end << std::endl;
              std::cout << "     lar_z_lim2 " << lar_z_lim2 << std::endl;
            }

            if (Project_track(lar_track,true)[0] < tms_x_lim1){
              std::cout << "    PROJECT_TRACK(LAr_TRACK,TRUE)[0] < TMS_X_LIM1" << std::endl;
              std::cout << "     Project_track(lar_track,true)[0] " << Project_track(lar_track,true)[0] << std::endl;
              std::cout << "     tms_x_lim1 " << tms_x_lim1 << std::endl;
            }
            if (Project_track(lar_track,true)[0] > tms_x_lim2){
              std::cout << "    PROJECT_TRACK(LAr_TRACK,TRUE)[0] > TMS_X_LIM2" << std::endl;
              std::cout << "     Project_track(lar_track,true)[0] " << Project_track(lar_track,true)[0] << std::endl;
              std::cout << "     tms_x_lim2 " << tms_x_lim2 << std::endl;
            } 
            if (Project_track(lar_track,true)[1] < tms_y_lim1){
              std::cout << "    PROJECT_TRACK(LAr_TRACK,TRUE)[1] < TMS_Y_LIM1" << std::endl;
              std::cout << "     Project_track(lar_track,true)[1] " << Project_track(lar_track,true)[1] << std::endl;
              std::cout << "     tms_y_lim1 " << tms_y_lim1 << std::endl;
            } 
            if (Project_track(lar_track,true)[1] > tms_y_lim2){
              std::cout << "    PROJECT_TRACK(LAr_TRACK,TRUE)[1] > TMS_Y_LIM2" << std::endl;
              std::cout << "     Project_track(lar_track,true)[1] " << Project_track(lar_track,true)[1] << std::endl;
              std::cout << "     tms_y_lim2 " << tms_y_lim2 << std::endl;
            } 

          return false;
    }
  }

	constexpr auto range_gramper_cm()
        {
	std::array<double, 73> Range_grampercm{
                {0.9833,   1.36,      1.786,    2.507,   3.321,   4.859,   6.598,   8.512,   10.58,   12.78,
                15.1,     17.52,     20.04,    25.31,   30.84,   36.59,   42.5,    54.73,   67.32,   86.66,
                106.3,    139.4,     172.5,    205.6,   238.5,   271.1,   303.5,   335.7,   367.7,   431.0,
                493.4,    555.2,     616.3,    736.8,   855.2,   1030.0,  1202.0,  1482.0,  1758.0,  2029.0,
                2297.0,   2562.0,    2825.0,   3085.0,  3343.0,  3854.0,  4359.0,  4859.0,  5354.0,  6333.0,
                7298.0,   8726.0,    10130.0,  12430.0, 14690.0, 16920.0, 19100.0, 21260.0, 23380.0, 25480.0,
                27550.0,  31610.0,   35580.0,  39460.0, 43260.0, 50620.0, 57680.0, 67780.0, 77340.0, 92220.0,
                1.06e+05, 1.188e+05, 1.307e+05}};
        for (double& value : Range_grampercm) {
                value /= NDLArTMSUniqueMatchRecoFiller::LArDen; // convert to cm
        }
    return Range_grampercm;
    }

        constexpr auto Range_grampercm = range_gramper_cm();
        constexpr std::array<double, 73> KE_MeV{
        {10.0,    12.0,    14.0,    17.0,    20.0,    25.0,    30.0,    35.0,    40.0,    45.0,
        50.0,    55.0,    60.0,    70.0,    80.0,    90.0,    100.0,   120.0,   140.0,   170.0,
        200.0,   250.0,   300.0,   350.0,   400.0,   450.0,   500.0,   550.0,   600.0,   700.0,
        800.0,   900.0,   1000.0,  1200.0,  1400.0,  1700.0,  2000.0,  2500.0,  3000.0,  3500.0,
        4000.0,  4500.0,  5000.0,  5500.0,  6000.0,  7000.0,  8000.0,  9000.0,  10000.0, 12000.0,
        14000.0, 17000.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0, 45000.0, 50000.0, 55000.0,
        60000.0, 70000.0, 80000.0, 90000.0, 1e+05,   1.2e+05, 1.4e+05, 1.7e+05, 2e+05,   2.5e+05,
        3e+05,   3.5e+05, 4e+05}};
        TGraph const KEvsR{73, Range_grampercm.data(), KE_MeV.data()};
        TSpline3 const KEvsR_spline3{"KEvsRS", &KEvsR};


  double NDLArTMSUniqueMatchRecoFiller::Muon_LAr_KE_Reco(const float trk_length) const
  { // Function for reconstructing the kinetic energy of a muon from the distance it would travel in LAr
	// From larreco/RecoAlg/TrackMomentumCalculator.cxx

	return KEvsR_spline3.Eval(trk_length);
  }

  void NDLArTMSUniqueMatchRecoFiller::Create_matches(std::vector<caf::SRNDTrackAssn> possibleMatches, bool Pandora, caf::StandardRecord &sr) const
  {
    std::sort(possibleMatches.begin(),possibleMatches.end(),Track_match_sorter);

    std::vector<caf::SRNDLArID> matched_lar; 
    std::vector<caf::SRTMSID> matched_tms; // stores LAr and TMS indices that have already been matched

    for (unsigned int match_idx = 0; match_idx < possibleMatches.size(); match_idx++) {
      std::cout << "Checking Match " << match_idx << "/" << possibleMatches.size() << std::endl;
      caf::SRNDTrackAssn track_match = possibleMatches[match_idx];
      std::cout << " LAr Interaction " << track_match.larid.ixn << " Index " << track_match.larid.idx << std::endl;
      std::cout << " TMS Interaction " << track_match.tmsid.ixn << " Index " << track_match.tmsid.idx << std::endl;
      double score = track_match.matchScore;
      std::cout << " Match Score " << score << std::endl;
      std::cout << " fCut " << f_cut << std::endl;
      if (score > f_cut) {
        std::cout << " SCORE > FCUT" << std::endl;
        continue; // this was break before, I think that was wrong
      }
      caf::SRNDLArID larid = track_match.larid;
      bool seen_lar = false; // checks if this LAr track has been matched already
      for (auto const seen_larid : matched_lar) { // I think this might be wrong. May need to fix
        if (seen_larid.ixn == larid.ixn && seen_larid.idx == larid.idx) {
          seen_lar = true;
          break;
        }
      }
      if (seen_lar) {
        continue;
      }
      caf::SRTMSID tmsid = track_match.tmsid;
      bool seen_tms = false; // checks if this TMS track has been matched already
      for (auto const seen_tmsid : matched_tms) {
        if (seen_tmsid.ixn == tmsid.ixn && seen_tmsid.idx == tmsid.idx) {
          seen_tms = true;
          break;
        }
      }
      if (seen_tms) {
        continue;
      }

      matched_tms.push_back(tmsid);
      matched_lar.push_back(larid);
      sr.nd.trkmatch.extrap.push_back(track_match); // adds successfully matched pair to StandardRecord of track matches
      sr.nd.trkmatch.nextrap += 1;
      caf::SRTrack joint_track = track_match.trk;
      caf::SRTrack lar_track;
      if (Pandora) {
         lar_track = sr.nd.lar.pandora[larid.ixn].tracks[larid.idx];
      } // if the Pandora flag is true, then the LAr tracks are listed within sr.nd.lar.pandora
      else {
         lar_track = sr.nd.lar.dlp[larid.ixn].tracks[larid.idx];
      }     // otherwise, they're listed within sr.nd.lar.dlp
      caf::SRTrack tms_track = sr.nd.tms.ixn[tmsid.ixn].tracks[tmsid.idx];     // this is the TMS track
      joint_track.start = lar_track.start;      // starting point of joint track is starting point of LAr track (Pandora or SPINE)
      joint_track.end = tms_track.end;          // ending point of joint track is ending point of TMS track
      joint_track.dir = lar_track.dir;          // starting direction of joint track is starting direction of LAr track (Pandora or SPINE)
      joint_track.enddir = tms_track.enddir;    // end direction of joint track is end direction of TMS track
      joint_track.time = lar_track.time;        // TODO: once we have reco LAr time working properly for both Pandora and SPINE this should be switched to lar_track.time from tms_track.time
      joint_track.Evis = lar_track.Evis + tms_track.Evis;
      joint_track.charge = tms_track.charge;

      // TODO: The following values pertain to the dead region between LAr and TMS. They may not remain accurate as geometry changes. Long-term solution is to directly query the geometry file instead
      float gap_dist = sqrt(pow(tms_track.start.x - lar_track.end.x,2)+pow(tms_track.start.y - lar_track.end.y,2)+pow(tms_track.start.z - lar_track.end.z,2));
      double LArTMSGap = TMSStart - ActiveLArEnd;
      double DeadLArDist = DeadLArEnd - ActiveLArEnd;
      double DeadLArFrac = DeadLArDist / LArTMSGap; // What fraction of the distance between LAr and TMS is dead LAr?
      joint_track.len_cm = lar_track.len_cm + tms_track.len_gcm2/LArDen + DeadLArFrac*gap_dist; // Divide the TMS areal density by the LAr density and add the dead LAr length

      // TODO: Split this straight line distance (gap_dist) into segments as it passes through each subsequent material - only the dead LAr has been implemented so far
      joint_track.len_gcm2 = (lar_track.len_cm + DeadLArFrac*gap_dist)*LArDen + tms_track.len_gcm2;
      // TODO: add the rest of the joint_track attributes (qual, truth, truthOverlap)
      double KE_mu = Muon_LAr_KE_Reco(joint_track.len_cm); // Gives muon kinetic energy in MeV
      float M_mu = 105.658; // Muon mass in MeV
      joint_track.E = KE_mu + M_mu;
    }
  }

  std::vector<caf::SRNDTrackAssn> NDLArTMSUniqueMatchRecoFiller::Compute_match_scores(const caf::SRNDLArInt ixn, const unsigned int ixn_lar, const unsigned int n_tracks, const unsigned int ixn_tms, const unsigned int itms, const double lar_z_cutoff, const caf::SRTrack tms_trk, caf::StandardRecord &sr, const cafmaker::Trigger &trigger, const float time_smear, std::set<int> matchIDs) const
  { // given a TMS track and a LAr interaction, computes the match scores between that TMS track and all LAr tracks in the interaction
    std::vector<caf::SRNDTrackAssn> potentialMatchList;

    for (unsigned int itrk = 0; itrk < n_tracks; itrk++) {
      std::cout << "   LAr track " << itrk << "/" << n_tracks << std::endl;
      caf::SRTrack trk = ixn.tracks[itrk];
      bool trueMatch = false;

      if (!Consider_LAr_track(trk,lar_z_cutoff)) {
        std::cout << "    LAr TRACK FAILED" << std::endl;
        continue; //skips the lar track if it isn't suitable according to the function
      }

      std::vector<float> proj_vec = Project_track(trk,true);
      std::cout << "    Track End Point Z " << trk.end.z << std::endl;
      std::cout << "    TMS Start Z " << tms_z_lim1 << std::endl;
      std::cout << "    Track Z Distance " << tms_z_lim1 - trk.end.z << std::endl;
      std::cout << "    Track End Dir Z " << trk.enddir.z << std::endl;
      std::cout << "    Z Distance / Dir Z " << (tms_z_lim1 - trk.end.z)/trk.enddir.z << std::endl;
      std::cout << "    Track End Point X " << trk.end.x << std::endl;
      std::cout << "    Track End Dir X " << trk.enddir.x << std::endl;
      std::cout << "    X Displacement " << trk.enddir.x*(tms_z_lim1 - trk.end.z)/trk.enddir.z << std::endl;
      std::cout << "    X Projection " << trk.end.x + trk.enddir.x*(tms_z_lim1 - trk.end.z)/trk.enddir.z << std::endl;
      std::cout << "    Calculated X Projection " << proj_vec[0] << std::endl;
      std::cout << "    TMS Start X " << tms_trk.start.x << std::endl;
      float delta_x = tms_trk.start.x - proj_vec[0];
      std::cout << "    Delta X " << delta_x << std::endl;
      std::cout << "    Sigma X " << sigma_x << std::endl;
      std::cout << "    X Term " << pow(delta_x/sigma_x,2) << std::endl;

      std::cout << "    Z Distance / Dir Z " << (tms_z_lim1 - trk.end.z)/trk.enddir.z << std::endl;
      std::cout << "    Track End Point Y " << trk.end.y << std::endl;
      std::cout << "    Track End Dir Y " << trk.enddir.y << std::endl;
      std::cout << "    Y Displacement " << trk.enddir.y*(tms_z_lim1 - trk.end.z)/trk.enddir.z << std::endl;
      std::cout << "    Y Projection " << trk.end.y + trk.enddir.y*(tms_z_lim1 - trk.end.z)/trk.enddir.z << std::endl;
      std::cout << "    Calculated Y Projection " << proj_vec[1] << std::endl;
      std::cout << "    TMS Start Y " << tms_trk.start.y << std::endl;

      float delta_y = tms_trk.start.y - proj_vec[1];
      std::cout << "    Delta Y " << delta_y << std::endl;
      std::cout << "    Sigma Y " << sigma_y << std::endl;
      std::cout << "    Y Term " << pow(delta_y/sigma_y,2) << std::endl;

      std::vector<float> angles = Angle_between_tracks(tms_trk,trk);
      angles[0] = std::copysign(angles[0],delta_x);
      std::cout << "    Delta Theta X " << angles[0] << std::endl;
      std::cout << "    Sigma Theta X " << sigma_angle_x << std::endl;
      std::cout << "    Theta X Term " << pow(angles[0]/sigma_angle_x,2) << std::endl;
      angles[1] = std::copysign(angles[1],delta_y); // if the TMS track has a larger X or Y value than the projected LAr track, then the X or Y angle is positive. Otherwise, it's negative
      std::cout << "    Delta Theta Y " << angles[1] << std::endl;
      std::cout << "    Sigma Theta Y " << sigma_angle_y << std::endl;
      std::cout << "    Theta Y Term " << pow(angles[1]/sigma_angle_y,2) << std::endl;

      double matchScore = std::numeric_limits<double>::max(); // initialize match score to max value

      double lar_time = 0;
      caf::SRVector3D start_pos;
      double delta_t = 0; // initialize time of LAr track and time difference between it and TMS track to 0

      float angle_x = angles[0];
      float angle_y = angles[1]; // x and y components of LAr and TMS tracks
      matchScore = pow(delta_x/sigma_x,2) + pow(delta_y/sigma_y,2) + pow(angle_x/sigma_angle_x,2)+ pow(angle_y/sigma_angle_y,2);
      bool truthFail = false;
      std::vector<float> tOv = trk.truthOverlap;
      std::vector<caf::TrueParticleID> truIDs = trk.truth;
      if (tOv.empty()) {
        std::cout << "    tOv EMPTY" << std::endl;
        truthFail = true;
      }
      if (truIDs.empty()) {
        std::cout << "    truIDs EMPTY" << std::endl;
        truthFail = true;
      }
      if (truIDs.size() != tOv.size()) {
        std::cout << "    truIDs.size() != tOv.size()" << std::endl;
        truthFail = true;
      }
      if (!truthFail) {
        int idx_max = std::distance(tOv.begin(),std::max_element(tOv.begin(),tOv.end()));
        caf::TrueParticleID partID = truIDs[idx_max]; // ID of true particle that makes up the majority of the track
        const auto& matchedPart = FindParticle(sr.mc,partID); // gets the particle object corresponding to the ID
        if (matchedPart == nullptr) {
          std::cout << "    matchedPart was a null pointer" << std::endl;
          truthFail = true;
        }
        //else { // commented this out because I need the reco particle ID not the true particle ID
        //  trk.part = matchedPart->G4ID; // save the particle ID to the track
        //}
      }

      if (use_time) {
        // this handles time-based matching
        std::cout << "    USING TIME" << std::endl;
        if (!use_smear_time) {// this occurs when the LAr time has filled in the file
          std::cout << "    NO SMEAR TIME" << std::endl;
          lar_time = trk.time;
        }
        if (use_smear_time) {// this triggers if we're using a file where the LAr time hasn't filled
          std::cout << "    USING SMEAR TIME" << std::endl;
          if (!truthFail) {
            int idx_max = std::distance(tOv.begin(),std::max_element(tOv.begin(),tOv.end()));
            caf::TrueParticleID partID = truIDs[idx_max]; // ID of true particle that makes up the majority of the track
            const auto& matchedPart = FindParticle(sr.mc,partID); // gets the particle object corresponding to the ID
            lar_time = matchedPart->time - 1e9*trigger.triggerTime_s - trigger.triggerTime_ns + time_smear; // adds gaussian smear to the true time with std 10 ns
            std::cout << "     LAr Time " << lar_time << std::endl;
            double tms_time = tms_trk.time;
            std::cout << "     TMS Time " << tms_time << std::endl;
            delta_t = tms_time - lar_time;
            std::cout << "     Delta T " << delta_t << std::endl;
            std::cout << "     Mean T " << mean_t << std::endl;
            std::cout << "     Sigma T " << sigma_t << std::endl;
            std::cout << "     T Term " << pow((delta_t-mean_t)/sigma_t,2) << std::endl;
            matchScore += pow((delta_t-mean_t)/sigma_t,2); // adds the time difference term to the matchScore

            // Following code is for checking if the particle IDs for matching tracks themselves match. This allows you to identify true matches
            std::vector<float> tOvTMS = tms_trk.truthOverlap;
            std::vector<caf::TrueParticleID> truIDsTMS = tms_trk.truth;
            if (!tOvTMS.empty() && !truIDsTMS.empty() && tOvTMS.size() == truIDsTMS.size()) {
              int idx_max_TMS = std::distance(tOvTMS.begin(),std::max_element(tOvTMS.begin(),tOvTMS.end()));
              caf::TrueParticleID partIDTMS = truIDsTMS[idx_max_TMS];
              const auto& TMSPart = FindParticle(sr.mc,partIDTMS);
              if (TMSPart != nullptr) {
                if (matchedPart->G4ID==TMSPart->G4ID) {
                  trueMatch = true;
                  matchIDs.insert(matchedPart->G4ID); // adds the ID to the set of matchIDs we're keeping track of. We already know the LAr and TMS track have the same ID due to the check above
                  std::cout << "     TRUE MATCH" << std::endl;
                }
              }
              else {
                std::cout << "     TMSPart was a null pointer" << std::endl;
              }
            }
            else {
              if (tOvTMS.empty()) {
                std::cout << "     tOvTMS was empty" << std::endl;
              }
              if (truIDsTMS.empty()) {
                std::cout << "     truIDsTMS was empty" << std::endl;
              }
              if (tOvTMS.size() != truIDsTMS.size()) {
                std::cout << "     tOvTMS and truIDsTMS not the same size" << std::endl;
                std::cout << "     tOvTMS size " << tOvTMS.size() << std::endl;
                std::cout << "     truIDsTMS size " << truIDsTMS.size() << std::endl;
              }
            }
          }
          else {
            std::cout << "     Problem finding the LAr particle" << std::endl;
          }
        }
      }
      else {
        std::cout << "     NOT USING TIME" << std::endl;
      }
      std::cout << "     Match Score " << matchScore << std::endl;
      std::cout << "     fCut " << f_cut << std::endl;
      caf::SRTMSID tmsid;
      tmsid.ixn = ixn_tms;
      tmsid.idx = itms;
      caf::SRNDLArID larid;
      larid.reco = caf::kPandoraNDLAr;
      larid.ixn = ixn_lar;
      larid.idx = itrk;

      caf::SRNDTrackAssn potential_match;
      if (use_time) {
        potential_match.matchType = caf::NDRecoMatchType::kUniqueWithTime;
      }
      else {
        potential_match.matchType = caf::NDRecoMatchType::kUniqueNoTime;
      }
      potential_match.tmsid = tmsid;
      potential_match.larid = larid;
      potential_match.matchScore = matchScore;
      potential_match.transdispl = sqrt(pow(delta_x,2)+pow(delta_y,2));
      potential_match.cosangdispl = cos(TMath::Pi()/180.0 * angles[2]);
      potential_match.trueMatch = trueMatch;
      potential_match.deltaX = delta_x;
      potential_match.deltaY = delta_y;
      potential_match.deltaThetaX = angles[0];
      potential_match.deltaThetaY = angles[1];
      potential_match.deltaT = delta_t;
      potentialMatchList.push_back(potential_match);
    }

    return potentialMatchList;
  }

  void NDLArTMSUniqueMatchRecoFiller::_FillRecoBranches(const cafmaker::Trigger &trigger,
                                             caf::StandardRecord &sr,
                                             const cafmaker::Params &/*par*/,
                                             const TruthMatcher */*truthMatcher*/) const
  {
    TRandom3 rng( static_cast<unsigned int>(trigger.triggerTime_ns ));

    std::vector<caf::SRNDTrackAssn> possiblePandoraMatches; // vector will store all possible matched tracks between Pandora and TMS
    std::vector<caf::SRNDTrackAssn> possibleSPINEMatches; // vector will store all possible matched tracks between SPINE and TMS

    double tms_z_cutoff = 20;
    double lar_z_cutoff = 20; // tracks must overlap last/first 20 cm of the detectors

    std::set<int> matchIDs; // will store the IDs of true matches found

    for (unsigned int ixn_tms = 0; ixn_tms < sr.nd.tms.nixn; ixn_tms++) {
      std::cout << "TMS Interaction " << ixn_tms << "/" << sr.nd.tms.nixn << std::endl;
      caf::SRTMSInt tms_int = sr.nd.tms.ixn[ixn_tms];
      unsigned int n_tms_tracks = tms_int.ntracks;
      std::cout << "NUMBER OF TMS TRACKS " << n_tms_tracks << std::endl;

      for (unsigned int itms = 0; itms < n_tms_tracks; itms++){
        std::cout << " TMS Track " << itms << "/" << n_tms_tracks << std::endl;
        caf::SRTrack tms_trk = tms_int.tracks[itms];

        if (!Consider_TMS_track(tms_trk,tms_z_cutoff)) {
          std::cout << " TMS TRACK FAILED" << std::endl;
          continue; // skips the TMS track if it isn't suitable according to the function
        }

        for (unsigned int ixn_pan = 0; ixn_pan < sr.nd.lar.npandora; ixn_pan++){
          std::cout << "  Pandora Interaction " << ixn_pan << "/" << sr.nd.lar.npandora << std::endl;
          caf::SRNDLArInt pan_int = sr.nd.lar.pandora[ixn_pan];
          unsigned int n_pan_tracks = pan_int.ntracks;

	  float smearTime = rng.Gaus(0.,10.); // Used for the cheated LAr time
          std::vector<caf::SRNDTrackAssn> panTrkAssns = Compute_match_scores(pan_int, ixn_pan, n_pan_tracks, ixn_tms, itms, lar_z_cutoff, tms_trk, sr, trigger, smearTime, matchIDs);

          copy(panTrkAssns.begin(), panTrkAssns.end(), back_inserter(possiblePandoraMatches));
        }

        for (unsigned int ixn_dlp = 0; ixn_dlp < sr.nd.lar.ndlp; ixn_dlp++){
          std::cout << "  SPINE Interaction " << ixn_dlp << "/" << sr.nd.lar.ndlp << std::endl;
          caf::SRNDLArInt dlp_int = sr.nd.lar.dlp[ixn_dlp];
          unsigned int n_dlp_tracks = dlp_int.ntracks;

	  float smearTime = rng.Gaus(0.,10.); // Used for the cheated LAr time 
          std::vector<caf::SRNDTrackAssn> dlpTrkAssns = Compute_match_scores(dlp_int, ixn_dlp, n_dlp_tracks, ixn_tms, itms, lar_z_cutoff, tms_trk, sr, trigger, smearTime, matchIDs);

          copy(dlpTrkAssns.begin(), dlpTrkAssns.end(), back_inserter(possibleSPINEMatches));
        }
      }
    }

    if (possiblePandoraMatches.size() > 0) {
      Create_matches(possiblePandoraMatches,true,sr); // tells the matcher that it's working with Pandora LAr tracks
      }

    if (possibleSPINEMatches.size() > 0) {
      Create_matches(possibleSPINEMatches,false,sr); // tells the matcher that it's not working with Pandora LAr tracks (therefore, SPINE tracks)
      }

    // Check that the script works as intended
    unsigned int num_matches = 0;
    unsigned int num_true_matches = 0;
    for (unsigned int match_no = 0; match_no < sr.nd.trkmatch.nextrap; match_no++){
      caf::SRNDTrackAssn match_pair = sr.nd.trkmatch.extrap[match_no];
      if (match_pair.matchType == caf::NDRecoMatchType::kUniqueWithTime || match_pair.matchType == caf::NDRecoMatchType::kUniqueNoTime){
        // One of the match pairs this algorithm found, so proceed with the counting
        num_matches += 1;
        if (match_pair.trueMatch){
          num_true_matches += 1;
          caf::SRNDLArID match_trk_id = match_pair.larid;
          caf::SRTrack matched_lar_trk = sr.nd.lar.Reco<caf::SRTrack>(match_trk_id);
          std::cout << "TOTAL RECO ENERGY " << match_pair.trk.E << std::endl;
          std::cout << "LAr RECO ENERGY " << matched_lar_trk.E << std::endl;
          std::vector<float> match_tOv = matched_lar_trk.truthOverlap;
          std::vector<caf::TrueParticleID> matched_truIDs = matched_lar_trk.truth;
          int match_idx_max = std::distance(match_tOv.begin(),std::max_element(match_tOv.begin(),match_tOv.end()));
          caf::TrueParticleID match_partID = matched_truIDs[match_idx_max]; // ID of true particle that makes up the majority of the track
          const auto& matchedParticle = FindParticle(sr.mc,match_partID); // gets the particle object corresponding to the ID
          float match_true_E = matchedParticle->p.E; // energy of the true particle
          std::cout << "TOTAL TRUE ENERGY " << match_true_E << std::endl;
          }
        }
      }
    std::cout << "TOTAL MATCHES FOUND THIS SPILL " << num_matches << std::endl;
    std::cout << "TRUE MATCHES FOUND THIS SPILL " << num_true_matches << std::endl;
    if (num_matches > 0){
      float pur = num_true_matches/num_matches;
      std::cout << "PURITY THIS SPILL " << pur << std::endl;
      }
    }

  // todo: this is a placeholder
  std::deque<cafmaker::Trigger> NDLArTMSUniqueMatchRecoFiller::GetTriggers(int /*triggerType*/, bool /*beamOnly*/) const
  {
    return std::deque<cafmaker::Trigger>();
  }

}
