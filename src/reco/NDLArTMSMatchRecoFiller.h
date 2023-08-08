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

      std::deque<Trigger> GetTriggers(int triggerType) const override;

    private:
      void MatchTracks(caf::StandardRecord &sr) const;

      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;
  };
}

#endif //ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H
