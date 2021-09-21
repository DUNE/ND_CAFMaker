/// \file MLNDLArRecoBranchFiller.h
///
/// Fill ND-LAr reco branches using DeepLearnPhysics machine learning based reconstruction.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    Sept. 2021


#ifndef ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H
#define ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H

#include <unordered_map>
#include <typeindex>

#include "IRecoBranchFiller.h"

#include "NDLArProductFiller.h"

namespace cafmaker
{

  /// Fill reco CAF branches using an H5 "summary" file from the ML reco.
  class MLNDLArRecoBranchFiller : public IRecoBranchFiller
  {
    public:
      MLNDLArRecoBranchFiller(const std::string &h5filename);

    protected:
      void _FillRecoBranches(std::size_t evtIdx,
                             caf::StandardRecord &sr,
                             const cafmaker::dumpTree &dt,
                             const cafmaker::params &par) const override;

    private:
      NDLArProductFiller<caf::SRTrack>  fTrackFiller;
      NDLArProductFiller<caf::SRShower> fShowerFiller;
  };  // class MLNDLArRecoBranchFiller

} // namespace cafmaker
#endif //ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H
