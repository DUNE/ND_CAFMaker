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

      std::vector<float_t> * fTrackStartX=0;
      std::vector<float_t> * fTrackStartY=0;
      std::vector<float_t> * fTrackStartZ=0;
      std::vector<float_t> * fTrackEndX=0;
      std::vector<float_t> * fTrackEndY=0;
      std::vector<float_t> * fTrackEndZ=0;

      std::vector<float_t> * fTrackStartPX=0;
      std::vector<float_t> * fTrackStartPY=0;
      std::vector<float_t> * fTrackStartPZ=0;
      std::vector<float_t> * fTrackEndPX=0;
      std::vector<float_t> * fTrackEndPY=0;
      std::vector<float_t> * fTrackEndPZ=0;

      std::vector<float_t> * fTrackLenF=0;
      std::vector<float_t> * fTrackLenB=0;
      std::vector<float_t> * fTrackPF=0;
      std::vector<float_t> * fTrackPB=0;
      std::vector<float_t> * fTrackAvgIonF=0;
      std::vector<float_t> * fTrackAvgIonB=0;

      std::vector<int> * fTrackIDNumber=0;
      std::vector<int> * fTrackNClusters=0;

      std::vector<int>     * fTrackPIDF=0;
      std::vector<float_t> * fTrackPIDProbF=0;
      std::vector<int>     * fTrackPIDB=0;
      std::vector<float_t> * fTrackPIDProbB=0;

      std::vector<float_t> * fECALClusterX=0;
      std::vector<float_t> * fECALClusterY=0;
      std::vector<float_t> * fECALClusterZ=0;

      std::vector<int> * fECALClusterIDNumber=0;

      std::vector<float_t> * fECALClusterEnergy=0;
      std::vector<int>     * fECALClusterNhits=0;

      std::vector<int> * fECALAssn_ClusterID=0;
      std::vector<int> * fECALAssn_TrackID=0;
    
  };

}

#endif //ND_CAFMAKER_NDGArRECOBRANCHFILLER_H
