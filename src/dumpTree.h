/*
 * dumpTree.h:
 *   Holder for the info from dumpTree output
 *
 */

#ifndef ND_CAFMAKER_DUMPTREE_H
#define ND_CAFMAKER_DUMPTREE_H

#include <cstdint>
#include <string>
#include <vector>

class TTree;

namespace cafmaker
{
  class dumpTree
  {
    public:
      void BindToTree(TTree * tree);

      int ievt = -1;
      int lepPdg = 0;
      int muonReco = -1;
      float lepKE = -999.;
      float muGArLen = -999.;
      float muECalLen = -999.;
      float hadTot = -999.;
      float hadCollar = -999.;
      float hadP = -999.;
      float hadN = -999.;
      float hadPip = -999.;
      float hadPim = -999.;
      float hadPi0 = -999.;
      float hadOther = -999.;
      float p3lep[3] = {-999., -999., -999.};
      float vtx[3] = {-999., -999., -999.};
      float muonExitPt[3] = {-999., -999., -999.};
      float muonExitMom[3] = {-999., -999., -999.};
      float muon_endpoint[3] = {-999., -999., -999.};

      // the arrays arrays are default-initialized to 0
      // (can only initialize to 0, not another value, with simple syntax).
      // but the size of the array to read is given by nFS,
      // and all this is meant to be overwritten when a record
      // is read out of the TTree anyway, so...
      int nFS = 0;
      int fsPdg[100] = {};
      float fsPx[100] = {};
      float fsPy[100] = {};
      float fsPz[100] = {};
      float fsE[100] = {};
      float fsTrkLen[100] = {};
      float fsTrkLenPerp[100] = {};

      std::vector <std::vector<std::vector < uint64_t>> > * geoEffThrowResults;
      std::string muon_endVolName;
  };
}

#endif //ND_CAFMAKER_DUMPTREE_H
