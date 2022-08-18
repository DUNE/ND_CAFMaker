/// Fill ND-LAr reco branches using DeepLearnPhysics machine learning based reconstruction.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

#ifndef ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H
#define ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H

#include <memory>

#include "IRecoBranchFiller.h"

namespace cafmaker
{
  class NDLArTMSMatchRecoFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDLArTMSMatchRecoFiller();

    private:
      void MatchTracks(caf::StandardRecord &sr) const;

      void _FillRecoBranches(std::size_t evtIdx,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par) const override;
  };
}

#endif //ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H
