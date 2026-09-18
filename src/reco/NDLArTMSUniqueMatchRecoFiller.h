/// Match ND-LAr and TMS tracks.
///
/// \author  Q. Weyrich <qweyrich@yorku.ca>
/// \date    July 2025

#ifndef ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H
#define ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H

#include "IRecoBranchFiller.h"
#include "MLNDLArRecoBranchFiller.h"
#include "TMSRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRTruthBranch.h"
#include "duneanaobj/StandardRecord/Navigate.h"
#include "TMath.h"
#include "TRandom3.h"
#include <set>

namespace cafmaker
{
  bool Track_match_sorter(const caf::SRNDTrackAssn trackMatch1, const caf::SRNDTrackAssn trackMatch2);

  class NDLArTMSUniqueMatchRecoFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDLArTMSUniqueMatchRecoFiller(const double sigmaX, const double sigmaY, const double sigmaThX, const double sigmaThY, const bool useTime, const double meanT, const double sigmaT, const double fCut, const bool useSmearTime);

      std::vector<double> Project_track(const caf::SRTrack track, const bool forward) const;

      std::vector<double> Angle_between_tracks(const caf::SRTrack tms_track, const caf::SRTrack lar_track) const;

      bool Consider_TMS_track(const caf::SRTrack tms_track, const double tms_z_cutoff) const;

      bool Consider_LAr_track(const caf::SRTrack lar_track, const double lar_z_cutoff) const;

      double Muon_LAr_KE_Reco(const float trk_length) const;

      void Create_matches(std::vector<caf::SRNDTrackAssn> possibleMatches, const bool Pandora, caf::StandardRecord &sr) const;

      std::vector<caf::SRNDTrackAssn> Compute_match_scores(const caf::SRNDLArInt ixn, const unsigned int ixn_lar, const unsigned int n_tracks, const unsigned int ixn_tms, const unsigned int itms, const double lar_z_cutoff, const caf::SRTrack tms_trk, caf::StandardRecord &sr, const cafmaker::Trigger &trigger, const float time_smear, std::set<int> matchIDs) const;

      std::deque<cafmaker::Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::Matcher; }

    private:
      void MatchTracks(caf::StandardRecord &sr) const;

      void _FillRecoBranches(const cafmaker::Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;
        
      std::vector<caf::SRNDTrackAssn> matchVectorSPINETMS;

      std::vector<caf::SRNDTrackAssn> matchVectorPandoraTMS;

      double sigma_x;
      double sigma_y;
      double sigma_angle_x;
      double sigma_angle_y;
      bool use_time;
      double mean_t;
      double sigma_t;
      double f_cut;

      // Dimensions of TMS and LAr fiducial volume [cm]. Best practice would be to query geometry instead of hard-coding
      double tms_x_lim1 = -352.0;
      double tms_x_lim2 = 352.0;
      double tms_y_lim1 = -386.4;
      double tms_y_lim2 = 115.9;
      double tms_z_lim1 = 1136.2;
      double tms_z_lim2 = 1831.4;

      double lar_x_lim1 = -347.848;
      double lar_x_lim2 = 347.848;
      double lar_y_lim1 = -216.671;
      double lar_y_lim2 = 82.9282;
      double lar_z_lim1 = 417.924;
      double lar_z_lim2 = 913.588;

      const double LArDen = 1.3954; // LAr density [g/cm3], from https://github.com/DUNE/dune-tms/blob/main/src/TMS_Constants.h c. Sep. 16, 2026
      const double ActiveLArEnd = 913.588; // End z-coordinate for LAr active volume [cm], from https://github.com/DUNE/dune-tms/blob/main/src/TMS_Constants.h c. Sep. 16, 2026
      const double DeadLArEnd = 937; // End z-coordinate for LAr instrumented volume [cm]
      const double TMSStart = 1117.75; // Start z-coordinate for TMS [cm], from https://github.com/DUNE/dune-tms/blob/main/src/TMS_Constants.h c. Sep. 16, 2026

  };
}

#endif //ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H
