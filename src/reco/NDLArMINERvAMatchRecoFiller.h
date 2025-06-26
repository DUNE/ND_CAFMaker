/// Matches NDLar tracks to Minerva tracks.
///
/// \author  N.Roy <noeroy@yorku.ca>
/// \date    Sep. 2024


#ifndef ND_CAFMAKER_NDLARMINERvAMATCHRECOFILLER_H
#define ND_CAFMAKER_NDLARMINERvAMATCHRECOFILLER_H

#include "IRecoBranchFiller.h"
#include "MLNDLArRecoBranchFiller.h"
#include "MINERvARecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"


namespace cafmaker
{
  class NDLArMINERvAMatchRecoFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDLArMINERvAMatchRecoFiller(double _z_extr, double _d_x, double _d_y, double _d_thetax, double d_theta_y);
      
      RecoFillerType FillerType() const override { return RecoFillerType::Matcher; }

      std::deque<Trigger> GetTriggers(int triggerType, bool beamOnly) const override;


    private:
      void MatchTracks(caf::StandardRecord &sr) const;
      bool Passes_cut(caf::SRTrack track_minerva, caf::SRTrack track_Lar, double &costheta, double &residual) const;
      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;

      // Matching parameters
      double z_extr;   // Extrapolated position compariton. Here it's the front of the Lar modules.
      double d_x;       // Maximum residual in x coordinate [cm];
      double d_y;       // Maximum residual in y coordinate [cm];
      double d_thetax; // Maximum Angle difference wrt to x axis [rad];
      double d_thetay; // Maximum Angle difference wrt to y axis [rad];
  };
}

#endif //ND_CAFMAKER_NDLARMINERvAMATCHRECOFILLER_H
