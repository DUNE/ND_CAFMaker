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

#include "reco/IRecoBranchFiller.h"
#include "reco/NDLArDLPH5DatasetReader.h"

namespace caf
{
  class SRTrueInteraction;
  class SRTrueParticle;
}

namespace cafmaker
{

  /// Fill reco CAF branches using an H5 "summary" file from the ML reco.
  class MLNDLArRecoBranchFiller : public IRecoBranchFiller
  {
    public:
      MLNDLArRecoBranchFiller(const std::string &h5filename);

      std::deque<Trigger> GetTriggers(int triggerType) const override;

    protected:
      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;

    private:
      void FillTracks(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                      const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueIxns,
                      const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                      const TruthMatcher * truthMatch,
                      caf::StandardRecord & sr) const;

      void FillShowers(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                       const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueIxns,
                       const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                       const TruthMatcher * truthMatch,
                       caf::StandardRecord & sr) const;

      void FillFlashes(const H5DataView<cafmaker::types::dlp::Flash> & flashes,
                       caf::StandardRecord & sr) const;
      
      void FillInteractions(const H5DataView<cafmaker::types::dlp::Interaction> &ixns,
                            const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueIxns,
                            const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                            const TruthMatcher * truthMatch,
                            caf::StandardRecord &sr) const;

      void FillParticles(const H5DataView<cafmaker::types::dlp::Particle> &particles,
                         const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueInxns,
                         const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                         const TruthMatcher * truthMatch,
                         caf::StandardRecord &sr) const;

      void FillTrueParticle(caf::SRTrueParticle & srTruePart,
                            const cafmaker::types::dlp::TrueParticle & truePartPassthrough,
                            const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles) const;

      void FillTrueInteraction(caf::SRTrueInteraction & srTrueInt,
                               const cafmaker::types::dlp::TrueInteraction & trueIntPassthrough) const;

      NDLArDLPH5DatasetReader fDSReader;
      mutable std::vector<cafmaker::Trigger> fTriggers;
      mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;    ///< the last trigger requested using _FillRecoBranches()
  };  // class MLNDLArRecoBranchFiller

} // namespace cafmaker
#endif //ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H
