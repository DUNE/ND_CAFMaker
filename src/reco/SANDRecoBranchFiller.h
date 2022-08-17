///
/// Fill reco branches using SAND data.
///

#ifndef ND_CAFMAKER_SANDRECOBRANCHFILLER_H
#define ND_CAFMAKER_SANDRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

#include "TFile.h"
#include "TTree.h"
#include "struct.h" //FIXME

namespace cafmaker
{

  class SANDRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      SANDRecoBranchFiller(const std::string &SANDRecoFilename);


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
