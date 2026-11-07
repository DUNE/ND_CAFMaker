#include "NDLArTMSUniqueMatchRecoFiller.h"
#include <cmath>
#include <set>

#include "Math/Functor.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSpline.h"

#include "duneanaobj/StandardRecord/Navigate.h"

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

  NDLArTMSUniqueMatchRecoFiller::NDLArTMSUniqueMatchRecoFiller(const double sigmaX, const double sigmaY,
                                                               const double sigmaThX, const double sigmaThY,
                                                               const bool useTime,
                                                               const double meanT, const double sigmaT,
                                                               const double fCut,
                                                               const bool useSmearTime)
    : IRecoBranchFiller("LArTMSMatcher")
  {
    fSigmaX = sigmaX;
    fSigmaY = sigmaY;
    fSigmaThetaX = sigmaThX;
    fSigmaThetaY = sigmaThY;
    fUseTime = useTime;
    fMeanT = meanT;
    fSigmaT = sigmaT;
    fScoreCut = fCut;
    fUseSmearTime = useSmearTime;
    // nothing to do
    SetConfigured(true);
  }

  std::array<float, 3> NDLArTMSUniqueMatchRecoFiller::Project_track(const caf::SRTrack &track, const bool forward) const
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

       proj_z = fTMSZLim1 - z;
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

       proj_z = z - fLArZLim2;
       proj_x = -dir_x*proj_z/dir_z + x;
       proj_y = -dir_y*proj_z/dir_z + y;
    }
    return {proj_x, proj_y, proj_z};
  }

  std::array<float, 3> NDLArTMSUniqueMatchRecoFiller::Angle_between_tracks(
    const caf::SRTrack &tms_track, const caf::SRTrack &lar_track) const
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
      return {angle_x,angle_y,angle_overall};
  }

  bool NDLArTMSUniqueMatchRecoFiller::Consider_TMS_track(const caf::SRTrack &tms_track, const double tms_z_cutoff) const
  {
    float x_start = tms_track.start.x;
    float y_start = tms_track.start.y;
    float z_start = tms_track.start.z;

    bool starts_in_FV = (x_start > fTMSXLim1) && (x_start < fTMSXLim2) &&
                        (y_start > fTMSYLim1) && (y_start < fTMSYLim2) &&
                        (z_start >= fTMSZLim1) && (z_start < fTMSZLim1 + tms_z_cutoff);
    if (!starts_in_FV)
      return false;

    // checks that direction would have allowed it to originate from LAr
    const auto projection = Project_track(tms_track,false);
    return    (projection[0] > fLArXLim1) && (projection[0] < fLArXLim2)
           && (projection[1] > fLArYLim1) && (projection[1] < fLArYLim2);
  }

  bool NDLArTMSUniqueMatchRecoFiller::Consider_LAr_track(const caf::SRTrack &lar_track, const double lar_z_cutoff) const
  {
    float x_start = lar_track.start.x;
    float y_start = lar_track.start.y;
    float z_start = lar_track.start.z;

    float x_end = lar_track.end.x;
    float y_end = lar_track.end.y;
    float z_end = lar_track.end.z;

    bool contained_and_reaches_back =    (x_start > fLArXLim1)&&(x_start < fLArXLim2)
                                      && (y_start > fLArYLim1) && (y_start < fLArYLim2)
                                      && (z_start > fLArZLim1) && (z_start < fLArZLim2)
                                      // checks track begins within fiducial volume

                                      && (x_end > fLArXLim1) && (x_end < fLArXLim2)
                                      && (y_end > fLArYLim1) && (y_end < fLArYLim2)
                                      && (z_end > fLArZLim2 - lar_z_cutoff) && (z_end < fLArZLim2);  // checks track ends close enough to back of LAr
    if (!contained_and_reaches_back)
      return false;

    const auto projection = Project_track(lar_track,true);
    return    (projection[0] > fTMSXLim1)&&(projection[0] < fTMSXLim2)
           && (projection[1] > fTMSYLim1)&&(projection[1] < fTMSYLim2); // checks that direction would allow it to hit TMS
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
      caf::SRNDTrackAssn & track_match = possibleMatches[match_idx];
      double score = track_match.matchScore;
      if (score > fScoreCut) {
        continue; // this was break before, I think that was wrong
      }
      caf::SRNDLArID larid = track_match.larid;
      bool seen_lar = false; // checks if this LAr track has been matched already
      for (auto const seen_larid : matched_lar) {
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
      caf::SRTrack & joint_track = track_match.trk;
      const caf::SRTrack * lar_track = nullptr;
      if (Pandora) {
         lar_track = &sr.nd.lar.pandora[larid.ixn].tracks[larid.idx];
      } // if the Pandora flag is true, then the LAr tracks are listed within sr.nd.lar.pandora
      else {
         lar_track = &sr.nd.lar.dlp[larid.ixn].tracks[larid.idx];
      }     // otherwise, they're listed within sr.nd.lar.dlp

      const caf::SRTrack & tms_track = sr.nd.tms.ixn[tmsid.ixn].tracks[tmsid.idx];     // this is the TMS track
      joint_track.start = lar_track->start;      // starting point of joint track is starting point of LAr track (Pandora or SPINE)
      joint_track.end = tms_track.end;          // ending point of joint track is ending point of TMS track
      joint_track.dir = lar_track->dir;          // starting direction of joint track is starting direction of LAr track (Pandora or SPINE)
      joint_track.enddir = tms_track.enddir;    // end direction of joint track is end direction of TMS track
      joint_track.time = lar_track->time;        // TODO: once we have reco LAr time working properly for both Pandora and SPINE this should be switched to lar_track.time from tms_track.time
      joint_track.Evis = lar_track->Evis + tms_track.Evis;
      joint_track.charge = tms_track.charge;
      joint_track.part = lar_track->part;        // Joint track inherits LAr track's reco particle ID

      // TODO: The following values pertain to the dead region between LAr and TMS. They may not remain accurate as geometry changes. Long-term solution is to directly query the geometry file instead
      float gap_dist = sqrt(pow(tms_track.start.x - lar_track->end.x,2)+pow(tms_track.start.y - lar_track->end.y,2)+pow(tms_track.start.z - lar_track->end.z,2));
      double LArTMSGap = fTMSZLim1 - fLArZLim2;
      double DeadLArDist = fDeadLArEnd - fLArZLim2;
      double DeadLArFrac = DeadLArDist / LArTMSGap; // What fraction of the distance between LAr and TMS is dead LAr?
      joint_track.len_cm = lar_track->len_cm + tms_track.len_gcm2/LArDen + DeadLArFrac*gap_dist; // Divide the TMS areal density by the LAr density and add the dead LAr length

      // TODO: Split this straight line distance (gap_dist) into segments as it passes through each subsequent material - only the dead LAr has been implemented so far
      joint_track.len_gcm2 = (lar_track->len_cm + DeadLArFrac*gap_dist)*LArDen + tms_track.len_gcm2;
      // TODO: add the rest of the joint_track attributes (qual, truth, truthOverlap)
      double KE_mu = Muon_LAr_KE_Reco(joint_track.len_cm); // Gives muon kinetic energy in MeV
      float M_mu = 105.658; // Muon mass in MeV
      joint_track.E = KE_mu + M_mu;

      sr.nd.trkmatch.extrap.push_back(track_match); // adds successfully matched pair to StandardRecord of track matches
      sr.nd.trkmatch.nextrap += 1;
    }
  }

  std::vector<caf::SRNDTrackAssn>
  NDLArTMSUniqueMatchRecoFiller::Compute_match_scores(const caf::SRNDLArInt &ixn,
                                                      const unsigned int ixn_lar,
                                                      const unsigned int n_tracks,
                                                      const unsigned int ixn_tms,
                                                      const unsigned int itms,
                                                      const double lar_z_cutoff,
                                                      const caf::SRTrack & tms_trk,
                                                      caf::StandardRecord & sr,
                                                      const cafmaker::Trigger &trigger,
                                                      const float time_smear,
                                                      // std::set<int> matchIDs,   // this parameter is copied on invocation and the modifications to it below are never stored anywhere
                                                      const bool Pandora) const
  { // given a TMS track and a LAr interaction, computes the match scores between that TMS track and all LAr tracks in the interaction
    std::vector<caf::SRNDTrackAssn> potentialMatchList;

    for (unsigned int itrk = 0; itrk < n_tracks; itrk++) {
      const caf::SRTrack & trk = ixn.tracks[itrk];
      bool trueMatch = false;

      if (!Consider_LAr_track(trk,lar_z_cutoff)) {
        continue; //skips the lar track if it isn't suitable according to the function
      }

      const std::array<float, 3> proj_vec = Project_track(trk,true);
      float delta_x = tms_trk.start.x - proj_vec[0];
      float delta_y = tms_trk.start.y - proj_vec[1];

      std::array<float, 3> angles = Angle_between_tracks(tms_trk,trk);
      angles[0] = std::copysign(angles[0],delta_x);
      angles[1] = std::copysign(angles[1],delta_y); // if the TMS track has a larger X or Y value than the projected LAr track, then the X or Y angle is positive. Otherwise, it's negative

      double matchScore = std::numeric_limits<double>::max(); // initialize match score to max value

      double lar_time = 0;
      caf::SRVector3D start_pos;
      double delta_t = 0; // initialize time of LAr track and time difference between it and TMS track to 0

      float angle_x = angles[0];
      float angle_y = angles[1]; // x and y components of LAr and TMS tracks
      matchScore = pow(delta_x/fSigmaX,2) + pow(delta_y/fSigmaY,2) + pow(angle_x/fSigmaThetaX,2)+ pow(angle_y/fSigmaThetaY,2);
      bool truthFail = false;
      const std::vector<float> & tOv = trk.truthOverlap;
      const std::vector<caf::TrueParticleID> & truIDs = trk.truth;
      if (tOv.empty()) {
        truthFail = true;
      }
      if (truIDs.empty()) {
        truthFail = true;
      }
      if (truIDs.size() != tOv.size()) {
        truthFail = true;
      }
      if (!truthFail) {
        int idx_max = std::distance(tOv.begin(),std::max_element(tOv.begin(),tOv.end()));
        const caf::TrueParticleID & partID = truIDs[idx_max]; // ID of true particle that makes up the majority of the track
        const auto& matchedPart = FindParticle(sr.mc,partID); // gets the particle object corresponding to the ID
        if (matchedPart == nullptr) {
          truthFail = true;
        }
      }

      if (fUseTime) {
        // this handles time-based matching
        if (!fUseSmearTime) {// this occurs when the LAr time has filled in the file
          lar_time = trk.time;
          double tms_time = tms_trk.time;
          delta_t = tms_time - lar_time;
          matchScore += pow((delta_t-fMeanT)/fSigmaT,2); // adds the time difference term to the matchScore
        }
      }

      if (!truthFail) {
        int idx_max = std::distance(tOv.begin(),std::max_element(tOv.begin(),tOv.end()));
        const caf::TrueParticleID & partID = truIDs[idx_max]; // ID of true particle that makes up the majority of the track
        const auto& matchedPart = FindParticle(sr.mc,partID); // gets the particle object corresponding to the ID
        if (fUseTime && fUseSmearTime) {// this triggers if we're using a file where the LAr time hasn't filled
          lar_time = matchedPart->time - 1e9*trigger.triggerTime_s - trigger.triggerTime_ns + time_smear; // adds gaussian smear to the true time with std 10 ns
          double tms_time = tms_trk.time;
          delta_t = tms_time - lar_time;
          matchScore += pow((delta_t-fMeanT)/fSigmaT,2); // adds the time difference term to the matchScore
        }
        if (std::isnan(matchScore)) {// should only be possible if the match score is nan, but this trick may not work with certain compilers?
          LOG.WARNING() << " NAN MATCH SCORE (score = " << matchScore << ")\n";
          LOG.WARNING() << "   tms_time = " << tms_trk.time << "; lar_time = " << lar_time << " (trk.time = " << trk.time << "); delta_t = " << delta_t << "\n";
          LOG.WARNING() << " LAr ixn " << ixn_lar << " | LAr idx " << itrk << "\n";
          LOG.WARNING() << " TMS ixn " << ixn_tms << " | TMS idx " << itms << "\n";
          LOG.WARNING() << " True particle ID " << matchedPart->G4ID << "\n";
          LOG.WARNING() << " True particle time " << matchedPart->time << "\n";
          LOG.WARNING() << " --------------------" << "\n";
          continue;
        }
        // Following code is for checking if the particle IDs for matching tracks themselves match. This allows you to identify true matches
        // TODO: When someone adds secondary info this will need to be used again: std::vector<float> tOvTMS = tms_trk.truthOverlap;
        const std::vector<caf::TrueParticleID> & truIDsTMS = tms_trk.truth;
        // if (!tOvTMS.empty() && !truIDsTMS.empty() && tOvTMS.size() == truIDsTMS.size()) {
        if (truIDsTMS.size() == 1) {
          //int idx_max_TMS = std::distance(tOvTMS.begin(),std::max_element(tOvTMS.begin(),tOvTMS.end()));
          const caf::TrueParticleID & partIDTMS = truIDsTMS[0];//idx_max_TMS];
          const auto& TMSPart = FindParticle(sr.mc,partIDTMS);
          if (TMSPart != nullptr) {
            if (matchedPart->G4ID==TMSPart->G4ID) {
              trueMatch = true;
              // matchIDs.insert(matchedPart->G4ID); // adds the ID to the set of matchIDs we're keeping track of. We already know the LAr and TMS track have the same ID due to the check above
            }
          }
        }
      }
      caf::SRTMSID tmsid;
      tmsid.ixn = ixn_tms;
      tmsid.idx = itms;
      caf::SRNDLArID larid;
      if (Pandora) {
        larid.reco = caf::kPandoraNDLAr;
      }
      else {
        larid.reco = caf::kDeepLearnPhys;
      }
      larid.ixn = ixn_lar;
      larid.idx = itrk;

      caf::SRNDTrackAssn potential_match;
      if (fUseTime) {
        potential_match.matchType = caf::NDRecoMatchType::kUniqueWithTime;
      }
      else {
        potential_match.matchType = caf::NDRecoMatchType::kUniqueNoTime;
      }
      potential_match.tmsid = std::move(tmsid);
      potential_match.larid = std::move(larid);
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

    for (unsigned int ixn_tms = 0; ixn_tms < sr.nd.tms.nixn; ixn_tms++) {
      const caf::SRTMSInt & tms_int = sr.nd.tms.ixn[ixn_tms];
      unsigned int n_tms_tracks = tms_int.ntracks;

      for (unsigned int itms = 0; itms < n_tms_tracks; itms++){
        const caf::SRTrack & tms_trk = tms_int.tracks[itms];

        if (!Consider_TMS_track(tms_trk,tms_z_cutoff)) {
          continue; // skips the TMS track if it isn't suitable according to the function
        }

        for (unsigned int ixn_pan = 0; ixn_pan < sr.nd.lar.npandora; ixn_pan++){
          const caf::SRNDLArInt & pan_int = sr.nd.lar.pandora[ixn_pan];
          unsigned int n_pan_tracks = pan_int.ntracks;

          float smearTime = rng.Gaus(0.,10.); // Used for the cheated LAr time
          std::vector<caf::SRNDTrackAssn> panTrkAssns = Compute_match_scores(pan_int, ixn_pan, n_pan_tracks, ixn_tms, itms, lar_z_cutoff, tms_trk, sr, trigger, smearTime, true);

          copy(panTrkAssns.begin(), panTrkAssns.end(), back_inserter(possiblePandoraMatches));
        }

        for (unsigned int ixn_dlp = 0; ixn_dlp < sr.nd.lar.ndlp; ixn_dlp++){
          const caf::SRNDLArInt & dlp_int = sr.nd.lar.dlp[ixn_dlp];
          unsigned int n_dlp_tracks = dlp_int.ntracks;

          float smearTime = rng.Gaus(0.,10.); // Used for the cheated LAr time
          std::vector<caf::SRNDTrackAssn> dlpTrkAssns = Compute_match_scores(dlp_int, ixn_dlp, n_dlp_tracks, ixn_tms, itms, lar_z_cutoff, tms_trk, sr, trigger, smearTime, false);

          copy(dlpTrkAssns.begin(), dlpTrkAssns.end(), back_inserter(possibleSPINEMatches));
        }
      }
    }

    if (possiblePandoraMatches.size() > 0) {
      // we won't use possiblePandoraMatches afterwards, and Create_matches() alters its copy.
      // move()ing it instead sacrifices this copy of it, but then the method doesn't wind up creating a second copy
      Create_matches(std::move(possiblePandoraMatches),true,sr); // tells the matcher that it's working with Pandora LAr tracks
      }

    if (possibleSPINEMatches.size() > 0) {
      // ditto
      Create_matches(std::move(possibleSPINEMatches),false,sr); // tells the matcher that it's not working with Pandora LAr tracks (therefore, SPINE tracks)
      }

    // Check that the script works as intended
    unsigned int num_matches = 0;
    unsigned int num_true_matches = 0;
    for (unsigned int match_no = 0; match_no < sr.nd.trkmatch.nextrap; match_no++){
      const caf::SRNDTrackAssn & match_pair = sr.nd.trkmatch.extrap[match_no];
      if (match_pair.matchType == caf::NDRecoMatchType::kUniqueWithTime || match_pair.matchType == caf::NDRecoMatchType::kUniqueNoTime){
        // One of the match pairs this algorithm found, so proceed with the counting and energy reco
        num_matches += 1;
        if (match_pair.trueMatch)
          num_true_matches += 1;

        const caf::SRNDLArID & match_trk_id = match_pair.larid;
        const caf::SRTrack & matched_lar_trk = sr.nd.lar.Reco<caf::SRTrack>(match_trk_id);
        const std::vector<float> & match_tOv = matched_lar_trk.truthOverlap;
        const std::vector<caf::TrueParticleID> & matched_truIDs = matched_lar_trk.truth;
        if (match_tOv.empty() || matched_truIDs.empty() || match_tOv.size() != matched_truIDs.size())
          continue; // no usable truth info for this LAr track, so energy reco can't be compared to truth

        int match_idx_max = std::distance(match_tOv.begin(),std::max_element(match_tOv.begin(),match_tOv.end()));
        const caf::TrueParticleID & match_partID = matched_truIDs[match_idx_max]; // ID of true particle that makes up the majority of the track
        const auto& matchedParticle = FindParticle(sr.mc,match_partID); // gets the particle object corresponding to the ID
        if (!matchedParticle)
          continue;

        float start_x = matchedParticle->start_pos.x;
        float start_y = matchedParticle->start_pos.y;
        float start_z = matchedParticle->start_pos.z;
        bool ERecoForTrueMatchesOnly = false; // if true, this only prints out the energy reco performance for true matched track pairs - may be useful for debugging
        bool PrintRecoVsTrueE = true;
        if ((ERecoForTrueMatchesOnly)&&(!match_pair.trueMatch)){
          PrintRecoVsTrueE = false; // if you tell the algorithm to only show the energy reco performance for the true matched tracks, and you DO NOT have a true matched pair, then we should not print the following
        }
        if ((start_x > fLArXLim1)&&(start_x < fLArXLim2)
            &&(start_y > fLArYLim1)&&(start_y < fLArYLim2)
            &&(start_z > fLArZLim1)&&(start_z < fLArZLim2)){
              if (PrintRecoVsTrueE) {
                float match_true_E = matchedParticle->p.E; // energy of the true particle
                int match_true_PDG = matchedParticle->pdg; // PDG of the true particle
                LOG.VERBOSE() << "TOTAL RECO ENERGY " << match_pair.trk.E << "\n";
                LOG.VERBOSE() << "LAr RECO ENERGY " << matched_lar_trk.E << "\n";

                LOG.VERBOSE() << "TOTAL TRUE ENERGY " << match_true_E << "\n";
                LOG.VERBOSE() << "TRUE PARTICLE PDG " << match_true_PDG << "\n";
              }
            else {
              LOG.VERBOSE() << "FALSE MATCH: " << "\n";
              LOG.VERBOSE() << " LAr Interaction " << match_pair.larid.ixn << " LAr Index " << match_pair.larid.idx << "\n";
              LOG.VERBOSE() << " TMS Interaction " << match_pair.tmsid.ixn << " TMS Index " << match_pair.tmsid.idx << "\n";
              LOG.VERBOSE() << " Not displaying energy reco" << "\n";
            }
          }
      }
    }
    LOG.DEBUG() << "Number of matches " << num_matches << "\n";
    LOG.DEBUG() << "Number of true matches " << num_true_matches << "\n";
    if (num_matches > 0){
      float pur = static_cast<float>(num_true_matches)/num_matches;
      LOG.DEBUG() << "Purity " << pur << "\n";
    }
  }

  // todo: this is a placeholder
  std::deque<cafmaker::Trigger> NDLArTMSUniqueMatchRecoFiller::GetTriggers(int /*triggerType*/, bool /*beamOnly*/) const
  {
    return std::deque<cafmaker::Trigger>();
  }

}
