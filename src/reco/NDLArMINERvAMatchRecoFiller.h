/// Fill ND-LAr reco branches using DeepLearnPhysics machine learning based reconstruction.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

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
      NDLArMINERvAMatchRecoFiller();

      std::deque<Trigger> GetTriggers(int triggerType) const override;

    private:
      void MatchTracks(caf::StandardRecord &sr) const;
      bool Passes_cut(caf::SRTrack track_minerva, caf::SRTrack track_Lar, double &costheta, double &residual) const;
      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;
  };
}

#endif //ND_CAFMAKER_NDLARMINERvAMATCHRECOFILLER_H
