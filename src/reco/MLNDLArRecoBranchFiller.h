/// \file MLNDLArRecoBranchFiller.h
///
/// Fill ND-LAr reco branches using DeepLearnPhysics machine learning based reconstruction.
///


#ifndef ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H
#define ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

#include "NDLArSummaryH5.h"

namespace cafmaker
{

  /// Fill reco CAF branches using an H5 "summary" file from the ML reco.
  class MLNDLArRecoBranchFiller : public IRecoBranchFiller
  {
    public:
      MLNDLArRecoBranchFiller(const std::string &h5filename,
                              const std::string &h5dataset);

    protected:
      void _FillRecoBranches(caf::StandardRecord& sr, const cafmaker::dumpTree & dt, const cafmaker::params &par) const override;

    private:
      NDLArSummaryH5 fH5file;
  };  // class MLNDLArRecoBranchFiller

} // namespace cafmaker
#endif //ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H
