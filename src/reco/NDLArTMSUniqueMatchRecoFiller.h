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

namespace cafmaker
{
  bool Track_match_sorter(const caf::SRNDTrackAssn trackMatch1, const caf::SRNDTrackAssn trackMatch2);

  class NDLArTMSUniqueMatchRecoFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDLArTMSUniqueMatchRecoFiller(const double sigmaX, const double sigmaY, const bool singleAngle, const double sigmaTh, const double sigmaThX, const double sigmaThY, const bool useTime, const double meanT, const double sigmaT, const double fCut);

      std::vector<double> Project_track(const caf::SRTrack track, const bool forward) const;

      std::vector<double> Angle_between_tracks(const caf::SRTrack tms_track, const caf::SRTrack lar_track) const;

      bool Consider_TMS_track(const caf::SRTrack tms_track, const double tms_z_cutoff) const;

      bool Consider_LAr_track(const caf::SRTrack lar_track, const double lar_z_cutoff) const;

      void Create_matches(std::vector<caf::SRNDTrackAssn> possibleMatches, caf::StandardRecord &sr) const;

      std::vector<caf::SRNDTrackAssn> Compute_match_scores(const caf::SRNDLArInt ixn, const unsigned int ixn_lar, const unsigned int n_tracks, const unsigned int ixn_tms, const unsigned int itms, const double lar_z_cutoff, const caf::SRTrack tms_trk, caf::StandardRecord &sr, const Trigger &trigger) const;

      std::deque<Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::Matcher; }

    private:
      void MatchTracks(caf::StandardRecord &sr) const;

      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;
        
      std::vector<caf::SRNDTrackAssn> matchVectorSPINETMS;

      std::vector<caf::SRNDTrackAssn> matchVectorPandoraTMS;

      double sigma_x;
      double sigma_y;
      bool single_angle;
      double sigma_angle;
      double sigma_angle_x;
      double sigma_angle_y;
      bool use_time;
      double mean_t;
      double sigma_t;
      double f_cut;

      // Dimensions of TMS and LAr fiducial volume [cm]
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
  };
}

#endif //ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H
