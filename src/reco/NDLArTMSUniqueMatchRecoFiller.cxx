#include "NDLArTMSUniqueMatchRecoFiller.h"
#include <cmath>
#include <random>
#include <chrono>

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

  NDLArTMSUniqueMatchRecoFiller::NDLArTMSUniqueMatchRecoFiller(const double sigmaX, const double sigmaY, const bool singleAngle, const double sigmaTh, const double sigmaThX, const double sigmaThY, const bool useTime, const double meanT, const double sigmaT, const double fCut)
    : IRecoBranchFiller("LArTMSMatcher")
  {
    sigma_x = sigmaX;
    sigma_y = sigmaY;
    single_angle = singleAngle;
    sigma_angle = sigmaTh;
    sigma_angle_x = sigmaThX;
    sigma_angle_y = sigmaThY;
    use_time = useTime;
    mean_t = meanT;
    sigma_t = sigmaT;
    f_cut = fCut;
    // nothing to do
    SetConfigured(true);
  }

  std::vector<double> NDLArTMSUniqueMatchRecoFiller::Project_track(const caf::SRTrack track, const bool forward) const
  {
    double x, y, z;

    double dir_x, dir_y, dir_z;
    
    double proj_z;
    double proj_x;
    double proj_y;

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
    std::vector<double> proj_point = {proj_x, proj_y, proj_z};
    return proj_point;
  }

  std::vector<double> NDLArTMSUniqueMatchRecoFiller::Angle_between_tracks(const caf::SRTrack tms_track, const caf::SRTrack lar_track) const
  {
      double tms_dir_x = tms_track.dir.x;
      double tms_dir_y = tms_track.dir.y;
      double tms_dir_z = tms_track.dir.z;

      double lar_dir_x = lar_track.enddir.x;
      double lar_dir_y = lar_track.enddir.y;
      double lar_dir_z = lar_track.enddir.z;

      double xz_dot_prod = tms_dir_x*lar_dir_x + tms_dir_z*lar_dir_z;
      if (xz_dot_prod != 0) {
        double xz_dot_prod = xz_dot_prod/(sqrt(pow(tms_dir_x,2)+pow(tms_dir_z,2))*sqrt(pow(lar_dir_x,2)+pow(lar_dir_z,2)));
      }
      double yz_dot_prod = tms_dir_y*lar_dir_y + tms_dir_z*lar_dir_z;
      if (yz_dot_prod != 0) {
        double yz_dot_prod = yz_dot_prod/(sqrt(pow(tms_dir_y,2)+pow(tms_dir_z,2))*sqrt(pow(lar_dir_y,2)+pow(lar_dir_z,2)));
      }
      double dot_prod = tms_dir_x*lar_dir_x + tms_dir_y*lar_dir_y + tms_dir_z*lar_dir_z;
      double angle_x = 180.0/TMath::Pi() * acos(xz_dot_prod);
      double angle_y = 180.0/TMath::Pi() * acos(yz_dot_prod);
      double angle_overall = 180.0/TMath::Pi() * acos(dot_prod);
      std::vector<double> angles = {angle_x,angle_y,angle_overall};
      return angles;
  }

  bool NDLArTMSUniqueMatchRecoFiller::Consider_TMS_track(const caf::SRTrack tms_track, const double tms_z_cutoff) const
  {
    double x_start = tms_track.start.x;
    double y_start = tms_track.start.y;
    double z_start = tms_track.start.z;

    if ((x_start > tms_x_lim1)&&(x_start < tms_x_lim2) &&
        (y_start > tms_y_lim1)&&(y_start < tms_y_lim2) &&
        (z_start > tms_z_lim1)&&(z_start < tms_z_lim1 + tms_z_cutoff) && // checks track begins within fiducial volume and close enough to front
      
        (Project_track(tms_track,false)[0] > lar_x_lim1)&&(Project_track(tms_track,false)[0] < lar_x_lim2) &&
        (Project_track(tms_track,false)[1] > lar_y_lim1)&&(Project_track(tms_track,false)[1] < lar_y_lim2)) // checks that direction would have allowed it to originate from LAr
          {
            return true;
          } 
    else {
      return false;
    }
  }

  bool NDLArTMSUniqueMatchRecoFiller::Consider_LAr_track(const caf::SRTrack lar_track, const double lar_z_cutoff) const
  {
    double x_start = lar_track.start.x;
    double y_start = lar_track.start.y;
    double z_start = lar_track.start.z;

    double x_end = lar_track.end.x;
    double y_end = lar_track.end.y;
    double z_end = lar_track.end.z;

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
      return false;
    }
  }

  void NDLArTMSUniqueMatchRecoFiller::Create_matches(std::vector<caf::SRNDTrackAssn> possibleMatches, caf::StandardRecord &sr) const
  {
    std::sort(possibleMatches.begin(),possibleMatches.end(),Track_match_sorter);

    std::vector<caf::SRNDLArID> matched_lar; 
    std::vector<caf::SRTMSID> matched_tms; // stores LAr and TMS indices that have already been matched

    for (unsigned int match_idx = 0; match_idx < possibleMatches.size(); match_idx++) {
      caf::SRNDTrackAssn track_match = possibleMatches[match_idx];
      double score = track_match.matchScore;
      if (score > f_cut) {
        break;}
      caf::SRNDLArID larid = track_match.larid;
      bool seen_lar = false; // checks if this LAr track has been matched already
      for (auto const seen_larid : matched_lar) {
        if (seen_larid.ixn == larid.ixn && seen_larid.idx == larid.idx) {
          seen_lar = true;
          break;
        }
      }
      if (seen_lar) {
        continue;}
      caf::SRTMSID tmsid = track_match.tmsid;
      bool seen_tms = false; // checks if this TMS track has been matched already
      for (auto const seen_tmsid : matched_tms) {
        if (seen_tmsid.ixn == tmsid.ixn && seen_tmsid.idx == tmsid.idx) {
          seen_tms = true;
          break;
        }
      }
      if (seen_tms) {
        continue;}
      
      matched_tms.push_back(tmsid);
      matched_lar.push_back(larid);
      sr.nd.trkmatch.extrap.push_back(track_match); // adds successfully matched pair to StandardRecord of track matches
      sr.nd.trkmatch.nextrap += 1;
    }
  }

  std::vector<caf::SRNDTrackAssn> NDLArTMSUniqueMatchRecoFiller::Compute_match_scores(const caf::SRNDLArInt ixn, const unsigned int ixn_lar, const unsigned int n_tracks, const unsigned int ixn_tms, const unsigned int itms, const double lar_z_cutoff, const caf::SRTrack tms_trk, caf::StandardRecord &sr, const Trigger &trigger) const
  { // given a TMS track and a LAr interaction, computes the match scores between that TMS track and all LAr tracks in the interaction
    std::vector<caf::SRNDTrackAssn> potentialMatchList;

    for (unsigned int itrk = 0; itrk < n_tracks; itrk++)
    {
      caf::SRTrack trk = ixn.tracks[itrk];

      if (!Consider_LAr_track(trk,lar_z_cutoff)) {
        continue; //skips the lar track if it isn't suitable according to the function
      }
      
      std::vector<double> proj_vec = Project_track(trk,true);

      double delta_x = tms_trk.start.x - proj_vec[0];
      double delta_y = tms_trk.start.y - proj_vec[1];

      std::vector<double> angles = Angle_between_tracks(tms_trk,trk);

      double matchScore = std::numeric_limits<double>::max(); // initialize match score to max value

      double lar_time = 0;
      caf::SRVector3D start_pos;
      double delta_t = 0; // initialize time of LAr track and time difference between it and TMS track to 0

      if (single_angle) {
        double angle = angles[2]; // overall angle between LAr and TMS track
        matchScore = pow(delta_x/sigma_x,2) + pow(delta_y/sigma_y,2) + pow(angle/sigma_angle,2);
        // matchScore is a weighted sum of x and y distances between tracks (after projection) and angle between them
      }
      else {
        double angle_x = angles[0];
        double angle_y = angles[1]; // x and y components of LAr and TMS tracks
        matchScore = pow(delta_x/sigma_x,2) + pow(delta_y/sigma_y,2) + pow(angle_x/sigma_angle_x,2)+ pow(angle_y/sigma_angle_y,2);
        // same as above except there are two angles not just one
      }

      if (use_time) {
        // this handles time-based matching - using truth-level particle times for now instead of light in LAr
        // this code is bugged somehow, giving nonsense lar_time values
        // either these truth-level particles have wrong times or the IDs we're inputting into sr.mc.Particle are wrong
        std::vector<float> tOv = trk.truthOverlap;
        std::vector<caf::TrueParticleID> truIDs = trk.truth;
        int idx_max = std::distance(tOv.begin(),std::max_element(tOv.begin(),tOv.end()));
        caf::TrueParticleID partID = truIDs[idx_max]; // ID of true particle that makes up majority of track
        const auto& matchedPart = FindParticle(sr.mc,partID); // gets the particle object corresponding to the ID
        if (matchedPart != nullptr) {
	  unsigned seed	= std::chrono::high_resolution_clock::now().time_since_epoch().count();
          std::mt19937 engine(seed);
          std::normal_distribution<double> dist(0.0,10.0);
          double time_smear = dist(engine);
          lar_time = matchedPart->time - 1e9*trigger.triggerTime_s - trigger.triggerTime_ns + time_smear; // adds gaussian smear to the true time with std 10 ns
          start_pos = matchedPart->start_pos;
	  double tms_time = tms_trk.time;
          delta_t = lar_time - tms_time;
          matchScore += pow((delta_t-mean_t)/sigma_t,2); // adds the time difference term to the matchScore
        }
      }

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
      potential_match.angdispl = cos(TMath::Pi()/180.0 * angles[2]);

      caf::SRTrack joint_track = potential_match.trk;
      joint_track.start = trk.start;
      joint_track.end = tms_trk.end;
      joint_track.dir = trk.dir;
      joint_track.enddir = tms_trk.enddir;

      joint_track.time = tms_trk.time; // TODO: once we have LAr time working properly this should be switched to trk.time

      joint_track.Evis = trk.Evis + tms_trk.Evis;
      // TODO: add the rest of the joint_track attributes
      
      potentialMatchList.push_back(potential_match);
    }

    return potentialMatchList;
  }


  void
  NDLArTMSUniqueMatchRecoFiller::_FillRecoBranches(const Trigger &trigger,
                                             caf::StandardRecord &sr,
                                             const cafmaker::Params &par,
                                             const TruthMatcher *truthMatcher) const
  {
    std::vector<caf::SRNDTrackAssn> possiblePandoraMatches; // vector will store all possible matched tracks between Pandora and TMS
    std::vector<caf::SRNDTrackAssn> possibleSPINEMatches; // vector will store all possible matched tracks between SPINE and TMS

    double tms_z_cutoff = 20;
    double lar_z_cutoff = 20; // tracks must overlap last/first 20 cm of the detectors
    
    for (unsigned int ixn_tms = 0; ixn_tms < sr.nd.tms.nixn; ixn_tms++)
    {
      caf::SRTMSInt tms_int = sr.nd.tms.ixn[ixn_tms];
      unsigned int n_tms_tracks = tms_int.ntracks;
      
      for (unsigned int itms = 0; itms < n_tms_tracks; itms++)
      {
        caf::SRTrack tms_trk = tms_int.tracks[itms];

        if (!Consider_TMS_track(tms_trk,tms_z_cutoff)) {
          continue; // skips the TMS track if it isn't suitable according to the function
        }

        for (unsigned int ixn_pan = 0; ixn_pan < sr.nd.lar.npandora; ixn_pan++)
        {
          caf::SRNDLArInt pan_int = sr.nd.lar.pandora[ixn_pan];
          unsigned int n_pan_tracks = pan_int.ntracks;
          
          std::vector<caf::SRNDTrackAssn> panTrkAssns = Compute_match_scores(pan_int, ixn_pan, n_pan_tracks, ixn_tms, itms, lar_z_cutoff, tms_trk, sr, trigger);

          copy(panTrkAssns.begin(), panTrkAssns.end(), back_inserter(possiblePandoraMatches));
        }

        for (unsigned int ixn_dlp = 0; ixn_dlp < sr.nd.lar.ndlp; ixn_dlp++)
        {
          caf::SRNDLArInt dlp_int = sr.nd.lar.dlp[ixn_dlp];
          unsigned int n_dlp_tracks = dlp_int.ntracks;

          std::vector<caf::SRNDTrackAssn> dlpTrkAssns = Compute_match_scores(dlp_int, ixn_dlp, n_dlp_tracks, ixn_tms, itms, lar_z_cutoff, tms_trk, sr, trigger);

          copy(dlpTrkAssns.begin(), dlpTrkAssns.end(), back_inserter(possibleSPINEMatches));
        }
      }
    }
    if (possiblePandoraMatches.size() > 0) {
      Create_matches(possiblePandoraMatches,sr);
      }

    if (possibleSPINEMatches.size() > 0) {
      Create_matches(possibleSPINEMatches,sr);
      }
  }
  // todo: this is a placeholder
  std::deque<Trigger> NDLArTMSUniqueMatchRecoFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    return std::deque<Trigger>();
  }

}
