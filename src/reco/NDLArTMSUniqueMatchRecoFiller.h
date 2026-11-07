/// Match ND-LAr and TMS tracks.
///
/// \author  Q. Weyrich <qweyrich@yorku.ca>
/// \date    July 2025

#ifndef ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H
#define ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H

#include "IRecoBranchFiller.h"
#include "TMSRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRNDTrackAssn.h"

#include <array>

namespace cafmaker
{
  bool Track_match_sorter(const caf::SRNDTrackAssn trackMatch1, const caf::SRNDTrackAssn trackMatch2);

  class NDLArTMSUniqueMatchRecoFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDLArTMSUniqueMatchRecoFiller(const double sigmaX, const double sigmaY, const double sigmaThX, const double sigmaThY, const bool useTime, const double meanT, const double sigmaT, const double fCut, const bool useSmearTime);

      std::array<float, 3> Project_track(const caf::SRTrack &track, const bool forward) const;

      std::array<float, 3> Angle_between_tracks(const caf::SRTrack &tms_track, const caf::SRTrack &lar_track) const;

      bool Consider_TMS_track(const caf::SRTrack &tms_track, const double tms_z_cutoff) const;

      bool Consider_LAr_track(const caf::SRTrack &lar_track, const double lar_z_cutoff) const;

      double Muon_LAr_KE_Reco(const float trk_length) const;

      void Create_matches(std::vector<caf::SRNDTrackAssn> possibleMatches, const bool Pandora, caf::StandardRecord &sr) const;

      std::vector<caf::SRNDTrackAssn> Compute_match_scores(const caf::SRNDLArInt &ixn,
                                                           unsigned int ixn_lar,
                                                           unsigned int n_tracks,
                                                           unsigned int ixn_tms,
                                                           unsigned int itms,
                                                           double lar_z_cutoff,
                                                           const caf::SRTrack &tms_trk,
                                                           caf::StandardRecord &sr,
                                                           const Trigger &trigger,
                                                           float time_smear,
                                                           bool Pandora) const;

      std::deque<cafmaker::Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::Matcher; }

      static constexpr double LArDen = 1.3954; // LAr density [g/cm3], from https://github.com/DUNE/dune-tms/blob/main/src/TMS_Constants.h c. Sep. 16, 2026


    private:
      void _FillRecoBranches(const cafmaker::Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;

      double fSigmaX;
      double fSigmaY;
      double fSigmaThetaX;
      double fSigmaThetaY;
      bool fUseTime;
      double fMeanT;
      double fSigmaT;
      double fScoreCut;
      bool fUseSmearTime;

      // Dimensions of TMS and LAr fiducial volume [cm]. TODO: query geometry instead of hard-coding
      double fTMSXLim1 = -400.0;//previous value (c. N5p1) -352.0;
      double fTMSXLim2 = 400.0;//previous value (c. N5p1) 352.0;
      double fTMSYLim1 = -570.2;//previous value (c. N5p1) -386.4 (largest difference)
      double fTMSYLim2 = 129.8;//previous value (c. N5p1) 115.9;
      double fTMSZLim1 = 1111.25;//1117.75;//1136.2; // 1136.2 came from TMSConstants.h but seems to be wrong, 1117.75 seems to be the correct value for N5p1 but the majority of N5p3 tracks start at 1111.25
      double fTMSZLim2 = 1863.5;//previous value (c.N5p1) 1831.4;
      // These are estimates of the boundaries of TMS using the N5p3 geometry (GDML file located here: https://github.com/DUNE/ND_Production/blob/main/geometry/nd_hall_with_lar_tms_sand_stt1_v2026.03.06.gdml).
      // They were obtained by converting the GDML to root using the instructions on the README here: https://github.com/PandoraPFAOrg/LArRecoND, and then modifying line 25 of
      // /exp/dune/app/users/weyrich/investigate_geo/investigate_geo.cpp to point to that root file, before running root -l /exp/dune/app/users/weyrich/investigate_geo/investigate_geo.cpp\(\"volTMS_0\"\) to get
      // the dimensions and coordinates of the centre of the volume volTMS_0. The exception is fTMSZLim1, which comes from printing the z-coordinate of the start point of all the TMS tracks in a TMS-Reco file

      // I am not 100% sure that volTMS_0 actually describes the active volume of TMS, or if it also includes some surrounding material. The z-coordinate of the front face of volTMS_0 is listed as 1108.5, which is
      // three cm less than the actual start point of the majority of the tracks, so I hope that the two agree with each other sufficiently well. But this uncertainty, and the complicated process needed to define
      // the geometry, is another reason the script should be set up to query the geometry for us to get the dimensions automatically

      double fLArXLim1 = -347.848;
      double fLArXLim2 = 347.848;
      double fLArYLim1 = -216.671;
      double fLArYLim2 = 82.9282;
      double fLArZLim1 = 417.924;
      double fLArZLim2 = 913.588;

      const double fDeadLArEnd = 937; // End z-coordinate for LAr instrumented volume [cm], from https://github.com/DUNE/dune-tms/blob/main/src/TMS_Constants.h c. Sep. 16, 2026

  };
}

#endif //ND_CAFMAKER_NDLARTMSUNIQUEMATCHRECOFILLER_H
