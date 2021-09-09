/// \file Params.h
///
///  Parameters extracted from command line and passed around


#ifndef ND_CAFMAKER_PARAMS_H
#define ND_CAFMAKER_PARAMS_H

namespace cafmaker
{
  // params will be extracted from command line, and passed to the reconstruction
  struct params {
    double OA_xcoord;
    bool fhc, grid, IsGasTPC;
    int seed, run, subrun, first, n, nfiles;
    double trk_muRes, LAr_muRes, ECAL_muRes;
    double em_const, em_sqrtE;
    double michelEff;
    double CC_trk_length;
    double pileup_frac, pileup_max;
    double gastpc_len, gastpc_B, gastpc_padPitch, gastpc_X0;
  };
}

#endif //ND_CAFMAKER_PARAMS_H
