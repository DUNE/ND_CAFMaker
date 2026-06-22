/// \file MLNDLArRecoBranchFiller.h
///
/// Fill ND-LAr reco branches using DeepLearnPhysics machine learning based reconstruction.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    Sept. 2021


#ifndef ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H
#define ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H

#include <unordered_map>
#include <map>
#include <typeindex>

#include "reco/IRecoBranchFiller.h"
#include "reco/NDLArDLPH5DatasetReader.h"

namespace caf
{
  class SRTrueInteraction;
  class SRTrueParticle;
  class SRRecoParticleID;
  class SRRecoParticle;
}

namespace cafmaker
{

  /// Fill reco CAF branches using an H5 "summary" file from the ML reco.
  class MLNDLArRecoBranchFiller : public IRecoBranchFiller
  {
    public:
      MLNDLArRecoBranchFiller(const std::string &h5filename);

      std::deque<Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

      bool IsBeamTrigger(int triggerType) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::BaseReco; }

      

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
      mutable std::map<int, int> fEntryMap; //Map of the filtered trigger entries stored in the caf file

      /// @brief Helper class to map from SPINE track ID to SRRecoParticle indices
      ///
      /// This mapper maintains a mapping between SPINE track IDs (from the ML reconstruction)
      /// and the corresponding indices in the CAF StandardRecord:
      /// - The interaction index in sr.common.ixn.dlp
      /// - The particle index in sr.common.ixn.dlp.part.dlp
      ///
      /// It is used to establish the linkage between low-level SPINE reconstruction
      /// and high-level SRRecoParticle objects in the CAF format.
      class MLNDLArRecoParticleMapper: public Loggable 
      {
        /// @brief Construct a new mapper with optional logging threshold
        /// @param logThresh Logging threshold (default: WARNING)
        public:
          MLNDLArRecoParticleMapper(cafmaker::Logger::THRESHOLD logThresh=cafmaker::Logger::THRESHOLD::WARNING) : Loggable("MLNDLArRecoParticleMapper", logThresh) {}

          /// @brief Get the SRRecoParticleID for a given SPINE track ID
          /// @param partID The SPINE track ID from ML reconstruction
          /// @return SRRecoParticleID containing the interaction index, collection type (kSPINE), and particle index
          ///
          /// @note This function will FATAL log and abort if the track ID is not found in the map.
          caf::SRRecoParticleID GetRecoParticleID(int64_t partID) const;

          /// @brief Get a reference to the SRRecoParticle for a given SPINE track ID
          /// @param sr The StandardRecord containing the reconstructed data
          /// @param partID The SPINE track ID from ML reconstruction
          /// @return Reference to the corresponding SRRecoParticle
          ///
          /// @note This function will FATAL log and abort if the track ID is not found in the map.
          caf::SRRecoParticle& GetRecoParticle(caf::StandardRecord & sr, int64_t partID) const;

          /// @brief Access or create a mapping for a SPINE track ID
          /// @param partID The SPINE track ID from ML reconstruction
          /// @return Reference to the (ixn_idx, prt_idx) pair
          std::pair<size_t, size_t>& operator[](int64_t partID) { return fParticleMap[partID]; }

          /// @brief Clear all particle mappings
          void Reset() { fParticleMap.clear(); }

        private:
          /// @brief Internal map from SPINE track ID to (interaction index, particle index) pairs
          std::map<int64_t, std::pair<size_t, size_t>> fParticleMap;
      };

      mutable MLNDLArRecoParticleMapper fParticleMapper; ///< helper object to map from SPINE track ID to (sr.common.ixn.dlp, sr.common.ixn.dlp.part.dlp) indices for the corresponding SRRecoParticle

      
  };  // class MLNDLArRecoBranchFiller

} // namespace cafmaker
#endif //ND_CAFMAKER_MLNDLARRECOBRANCHFILLER_H
