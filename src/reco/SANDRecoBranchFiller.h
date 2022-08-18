/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi
/// \date    Apr. 2022
///

#ifndef ND_CAFMAKER_SANDRECOBRANCHFILLER_H
#define ND_CAFMAKER_SANDRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

#include "TFile.h"
#include "TTree.h"
#include "struct.h"

namespace cafmaker
{

  class SANDRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      SANDRecoBranchFiller(const std::string &SANDRecoFilename,
			   const std::string &SANDLibFilename);

    private:
      void _FillRecoBranches(std::size_t evtIdx,
			     caf::StandardRecord &sr,
			     const cafmaker::Params &par) const override;

      TFile* fSANDRecoFile;
      TTree* fTree;
      struct event* fEvent;
  };

}

#endif //ND_CAFMAKER_SANDRECOBRANCHFILLER_H
