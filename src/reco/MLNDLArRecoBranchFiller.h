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

#include "NDLArDLPH5DatasetReader.h"

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
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;

    private:
      void FillTracks(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                      caf::StandardRecord & sr) const;

      void FillShowers(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                       caf::StandardRecord & sr) const;

      void FillTruth(const H5DataView<cafmaker::types::dlp::TrueParticle> & trueParticles,
                     const H5DataView<cafmaker::types::dlp::TrueInteraction> & trueInxns,
                     caf::StandardRecord &sr) const;

      NDLArDLPH5DatasetReader fDSReader;

  };  // class MLNDLArRecoBranchFiller

} // namespace cafmaker
#endif //ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H
