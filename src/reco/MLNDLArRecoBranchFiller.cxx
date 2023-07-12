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
    FillTruth(N, evtIdx, trueParticles, trueInteractions, sr);

    H5DataView<cafmaker::types::dlp::Particle> particles = fDSReader.GetProducts<cafmaker::types::dlp::Particle>(evtIdx);

    FillTracks(N, evtIdx, particles, sr);  
    FillShowers(N, evtIdx, particles, sr); 

    H5DataView<cafmaker::types::dlp::Interaction> interactions = fDSReader.GetProducts<cafmaker::types::dlp::Interaction>(evtIdx);
    FillInteractions(N, evtIdx, interactions, sr);
   
    //Resize here does not work 
  //  sr.nd.lar.dlp.resize(N);
  //  sr.mc.nu.resize(N);
   // sr.common.ixn.dlp.resize(N);

    //Fill ND-LAr specificinfo in the meta branch
    sr.meta.nd_lar.enabled = true;
    sr.meta.nd_lar.run = par().runInfo().run();
    sr.meta.nd_lar.subrun = par().runInfo().subrun();
    sr.meta.nd_lar.event = evtIdx;
  
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTruth(std::size_t N, std::size_t evtIdx, const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                                          const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueInxns,
                                          caf::StandardRecord &sr) const
  {
    sr.mc.nu.resize(N);


    //Filling truth information for every interaction
    for (const auto & trueInx : trueInxns)
    {
      caf::SRTrueInteraction true_interaction;
      true_interaction.vtx.x = trueInx.vertex[0];
      true_interaction.vtx.y = trueInx.vertex[1];
      true_interaction.vtx.z = trueInx.vertex[2];
      if(trueInx.nu_current_type == 0){true_interaction.iscc = false;}
      else{true_interaction.iscc = true;} 
      true_interaction.E = trueInx.nu_energy_init;
      true_interaction.nprim = trueInx.num_primaries;
      sr.mc.nu.push_back(std::move(true_interaction));

      //Fill only true primaries of each interaction 
      for (int i =0; i < trueInx.num_primaries; i++){
        for (const auto & truePart : trueParticles)
        { 
          if(truePart.is_primary)
          {
            caf::SRTrueParticle true_particle;
            true_particle.interaction_id = trueInx.id;
            true_particle.start_pos = {truePart.start_point[0], truePart.start_point[1], truePart.start_point[2]};
            true_particle.end_pos = {truePart.end_point[0], truePart.end_point[1], truePart.end_point[2]};
            true_particle.p.E = truePart.depositions_sum;
            true_particle.p.px = truePart.momentum[0];
            true_particle.p.py = truePart.momentum[1];
            true_particle.p.pz = truePart.momentum[2];
            true_particle.interaction_id = truePart.interaction_id;
            true_particle.ancestor_id.ixn = truePart.interaction_id;
            true_particle.ancestor_id.type = caf::TrueParticleID::kPrimary;
            true_particle.ancestor_id.part = trueInx.particle_ids[i];
            sr.mc.nu[evtIdx].prim.push_back(std::move(true_particle));
          }
        }
      }
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillInteractions(std::size_t N, std::size_t evtIdx, const H5DataView<cafmaker::types::dlp::Interaction> &Inxns,
                                           caf::StandardRecord &sr) const
  {
    sr.common.ixn.dlp.resize(N);

    //filling interactions
    for (const auto & inx : Inxns)
    {

     
      caf::SRInteraction interaction;
      // fill interaction variables
      interaction.vtx  = {inx.vertex[0], inx.vertex[1], inx.vertex[2]};
      interaction.dir.lngtrk  = {0., 0., 0.};
      //this does not work
      sr.common.ixn.dlp.push_back(std::move(interaction)); //make sure the variables correspond to the event index inside dlp vector
     
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTracks(std::size_t N, std::size_t evtIdx, const H5DataView<cafmaker::types::dlp::Particle> & particles,
                                           caf::StandardRecord &sr) const
  {
    sr.nd.lar.dlp.resize(N);
    for (const auto & part : particles)
    {
      // only choose 'particles' that correspond to Track semantic type
      if (part.semantic_type != types::dlp::kTrack)
        continue;


      caf::SRTrack track;
      // fill shower variables
      track.Evis = part.depositions_sum;
      track.start = {part.start_point[0], part.start_point[1], part.start_point[2]};
      track.end = {part.end_point[0], part.end_point[1], part.end_point[2]}; 
      track.dir = {part.start_dir[0], part.start_dir[1], part.start_dir[2]};
      track.enddir = {part.end_dir[0], part.end_dir[1], part.end_dir[2]};
      track.len_cm = sqrt(pow((part.start_point[0]-part.end_point[0]),2) + pow((part.start_point[1]-part.end_point[1]),2) + pow((part.start_point[2]-part.end_point[2]),2));
       
      
      sr.nd.lar.dlp[evtIdx].tracks.push_back(std::move(track)); //make sure the variables correspond to the event index inside dlp vector
     
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillShowers(std::size_t N, std::size_t evtIdx, const H5DataView<cafmaker::types::dlp::Particle> & particles,
                                            caf::StandardRecord &sr) const
  { 
    sr.nd.lar.dlp.resize(N);
    for (const auto & part : particles)
    {
      if (part.semantic_type != types::dlp::kShower)
        continue;

      caf::SRShower shower;
      // fill shower variables
      shower.Evis = part.depositions_sum;
      shower.start = {part.start_point[0], part.start_point[1], part.start_point[2]};
      shower.direction = {part.start_dir[0], part.start_dir[1], part.start_dir[2]};
      
      sr.nd.lar.dlp[evtIdx].showers.push_back(std::move(shower)); //make sure the variables correspond to the event index inside dlp vector
     
    }
  }

} // namespace cafmaker
