/// Fill ND-LAr reco branches using DeepLearnPhysics machine learning based reconstruction.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

#ifndef ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H
#define ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H

#include <memory>
#include <iostream>
#include "TF1.h"

#include "IRecoBranchFiller.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

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
                             const cafmaker::dumpTree &dt,
                             const cafmaker::Params &par) const override;
  };
  
   TF1* DrawLines(float x1, float z1, float x2, float z2, float startpoint, float endpoint)
  {
    TF1 *fun = new TF1("fun","(x - [0])*([3] - [2])/([1] - [0]) + [2]", startpoint, endpoint);
    fun->SetParameter(0,z1);
    fun->SetParameter(1,z2);
    fun->SetParameter(2,x1);
    fun->SetParameter(3,x2);

    return fun;
  }

  float x1_lar, x2_lar, y1_lar, y2_lar, z1_lar, z2_lar;
  float x1_tms, x2_tms, y1_tms, y2_tms, z1_tms, z2_tms;
  float x_LAr, x_TMS;
  float slope_LAr, slope_TMS;
  float slope_TMS_LAr;
  float ang;
  float residual;
  float Costheta;

}

#endif //ND_CAFMAKER_NDLARTMSMATCHRECOFILLER_H
