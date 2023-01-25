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
      std::vector<float_t> * fTrackStartY;
      std::vector<float_t> * fTrackStartZ;
      std::vector<float_t> * fTrackEndX;
      std::vector<float_t> * fTrackEndY;
      std::vector<float_t> * fTrackEndZ;

      std::vector<float_t> * fTrackStartPx;
      std::vector<float_t> * fTrackStartPy;
      std::vector<float_t> * fTrackStartPz;
      std::vector<float_t> * fTrackEndPx;
      std::vector<float_t> * fTrackEndPy;
      std::vector<float_t> * fTrackEndPz;

      std::vector<float_t> * fTrackLenF;
      std::vector<float_t> * fTrackLenB;
      std::vector<float_t> * fTrackPF;
      std::vector<float_t> * fTrackPB;
      std::vector<float_t> * fTrackAvgIonF;
      std::vector<float_t> * fTrackAvgIonB;

      std::vector<int> * fTrackIDNumber;

      std::vector<std::vector<int>>     * fTrackPIDF;
      std::vector<std::vector<float_t>> * fTrackPIDProbF;
      std::vector<std::vector<int>>     * fTrackPIDB;
      std::vector<std::vector<float_t>> * fTrackPIDProbB;

      std::vector<float_t> * fECALClusterX;
      std::vector<float_t> * fECALClusterY;
      std::vector<float_t> * fECALClusterZ;

      std::vector<int> * fECALClusterIDNumber;

      std::vector<float_t> * fECALClusterEnergy;
      std::vector<int>     * fECALClusterNhits;

      std::vector<int> * fECALAssn_ClusterID;
      std::vector<int> * fECALAssn_TrackID;
    
  };

}

#endif //ND_CAFMAKER_NDGArRECOBRANCHFILLER_H
