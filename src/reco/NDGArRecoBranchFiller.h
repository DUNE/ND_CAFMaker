/// Fill NDGAr reco branches using NDGAr reco data
///
/// \author  F. Martinez Lopez
/// \date    Oct. 2022
///

#ifndef ND_CAFMAKER_NDGArRECOBRANCHFILLER_H
#define ND_CAFMAKER_NDGArRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

#include <vector>
#include <cmath>

class TFile;
class TTree;

namespace cafmaker
{
  class NDGArRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDGArRecoBranchFiller(const std::string &NDGArRecoFilename);

    private:
      void _FillRecoBranches(std::size_t evtIdx,
			     caf::StandardRecord &sr,
			     const cafmaker::Params &par) const override;

      TFile* fNDGArRecoFile;
      TTree* NDGArRecoTree;
      
      int                fEvent;
      std::vector<float_t> * fTrackStartX;
      std::vector<int>   * fTrackPIDCheatedF;
      std::vector<float_t> * fTrackLenF;
  };

}

#endif //ND_CAFMAKER_NDGArRECOBRANCHFILLER_H
