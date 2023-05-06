#include "MLNDLArRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "DLP_h5_classes.h"
#include "Params.h"

using namespace cafmaker::types::dlp;

namespace cafmaker
{
  // ------------------------------------------------------------------------------
  // todo: possibly build some mechanism for customizing the dataset names in the file here
  MLNDLArRecoBranchFiller::MLNDLArRecoBranchFiller(const std::string &h5filename)
    : fDSReader(h5filename,
                {{std::type_index(typeid(Particle)),         "particles"},
                 {std::type_index(typeid(Interaction)),      "interactions"},
                 {std::type_index(typeid(TrueParticle)),     "truth_particles"},
                 {std::type_index(typeid(TrueInteraction)),  "truth_interactions"},
                 {std::type_index(typeid(Event)),            "events"}})
  {
    // if we got this far, nothing bad happened trying to open the file or dataset
    SetConfigured(true);
    name = "LArML";
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::_FillRecoBranches(std::size_t evtIdx,
                                                  caf::StandardRecord &sr,
                                                  const cafmaker::Params &par) const

  {
    H5DataView<cafmaker::types::dlp::TrueParticle> trueParticles = fDSReader.GetProducts<cafmaker::types::dlp::TrueParticle>(evtIdx);
    H5DataView<cafmaker::types::dlp::TrueInteraction> trueInteractions = fDSReader.GetProducts<cafmaker::types::dlp::TrueInteraction>(evtIdx);
    FillTruth(trueParticles, trueInteractions, sr);

    H5DataView<cafmaker::types::dlp::Particle> particles = fDSReader.GetProducts<cafmaker::types::dlp::Particle>(evtIdx);

    FillTracks(particles, sr);  // todo: also need to pass along truth info...
    FillShowers(particles, sr); // todo: also need to pass along truth info...

    // todo: add these when new StandardRecord is in place
    H5DataView<cafmaker::types::dlp::Interaction> interactions = fDSReader.GetProducts<cafmaker::types::dlp::Interaction>(evtIdx);
    //FillParticles(particles, sr);
    //FillInteractions(interactions, sr);
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTruth(const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                                          const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueInxns,
                                          caf::StandardRecord &sr) const
  {
    // todo: fill in here
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTracks(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                                           caf::StandardRecord &sr) const
  {
    for (const auto & part : particles)
    {
      // only choose 'particles' that correspond to Track semantic type
      if (part.semantic_type != types::dlp::kTrack)
        continue;

      // now copy stuff from the Particle into the sr track ....
      caf::SRTrack track;
      track.Evis = part.depositions_sum;
      // etc.

      sr.nd.lar.tracks.push_back(std::move(track));
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillShowers(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                                            caf::StandardRecord &sr) const
  {
    // todo: fill in here
  }

} // namespace cafmaker
