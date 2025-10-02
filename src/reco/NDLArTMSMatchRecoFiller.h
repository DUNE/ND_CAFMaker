/// Fill ND-LAr reco branches using DeepLearnPhysics machine learning based reconstruction.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

#ifndef ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H
#define ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H

#include "IRecoBranchFiller.h"
#include "MLNDLArRecoBranchFiller.h"
#include "TMSRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

namespace cafmaker
{
  class NDLArTMSMatchRecoFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDLArTMSMatchRecoFiller();

      std::deque<Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::Matcher; }


    private:
      void MatchTracks(caf::StandardRecord &sr) const;

      bool Spatial_Match(caf::SRTrack track_lar, caf::SRTrack track_tms, 
                          double &xthreshold, double &ythreshold, double &theta_xthreshold, 
                          double &theta_ythreshold, double &residual, double &ang_residual) const;

      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;

      double Angular_Match(double xdir1_tms, double ydir1_tms, 
                            double zdir1_tms, double xdir2_lar, 
                            double ydir2_lar, double zdir2_lar) const;
  };
}

#endif //ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H
