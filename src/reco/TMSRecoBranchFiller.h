/// Fill TMS reco branches from TMS reco output.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

#ifndef ND_CAFMAKER_TMSRECOBRANCHFILLER_H
#define ND_CAFMAKER_TMSRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

#include "TFile.h"

namespace cafmaker
{
  class TMSRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      TMSRecoBranchFiller(const std::string & tmsRecoFilename);

    private:
      void _FillRecoBranches(std::size_t evtIdx,
                             caf::StandardRecord &sr,
                             const dumpTree &dt,
                             const cafmaker::Params &par) const override;

      TFile fTMSRecoFile;
  };

}
#endif //ND_CAFMAKER_TMSRECOBRANCHFILLER_H
