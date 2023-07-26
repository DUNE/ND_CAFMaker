
//  WARNING:
//    This file was autogenerated by h5_to_cpp.py.
//    Do not edit by hand!

#ifndef CAFMAKER_TYPES_DLP_DLP_H5_CLASSES_H
#define CAFMAKER_TYPES_DLP_DLP_H5_CLASSES_H

#include <array>

#include "H5Cpp.h"
#include "readH5/BufferView.h"

namespace cafmaker::types::dlp
{

  // This generic template will be overloaded
  // for every specific type that we create below
  template <typename T>
  H5::CompType BuildCompType();

  
  enum Pid_t : int64_t
  {
    kElectron = 1,
    kMuon = 2,
    kPhoton = 0,
    kPion = 3,
    kProton = 4,
  };
  
  
  enum SemanticType_t : int64_t
  {
    kDelta = 3,
    kGhost = 5,
    kLowEnergy = 4,
    kMichel = 2,
    kShower = 0,
    kTrack = 1,
    kUnknown = 6,
  };
  
  struct Event;
  struct Interaction;
  struct Particle;
  struct TrueInteraction;
  struct TrueParticle;
  
  struct Event
  {
    hdset_reg_ref_t index;
    hdset_reg_ref_t run_info;
    hdset_reg_ref_t meta;
    hdset_reg_ref_t truth_particles;
    hdset_reg_ref_t particles;
    hdset_reg_ref_t truth_interactions;
    hdset_reg_ref_t interactions;
    
    void SyncVectors();
    
    template <typename T>
    const hdset_reg_ref_t& GetRef() const
    {
      if constexpr(std::is_same_v<T, TrueParticle>) return truth_particles;
      else if(std::is_same_v<T, Particle>) return particles;
      else if(std::is_same_v<T, TrueInteraction>) return truth_interactions;
      else if(std::is_same_v<T, Interaction>) return interactions;
    }
    
  };
  
  
  struct Interaction
  {
    int64_t crthit_id;
    uint8_t crthit_matched;
    int64_t crthit_matched_particle_id;
    int64_t flash_hypothesis;
    int64_t flash_id;
    double flash_time;
    int64_t flash_total_pE;
    uint8_t fmatched;
    int64_t id;
    int64_t image_id;
    bool is_contained;
    bool is_neutrino;
    bool is_principal_match;
    BufferView<int64_t> match;
    BufferView<float> match_overlap;
    uint8_t matched;
    int64_t nu_id;
    int64_t num_particles;
    int64_t num_primaries;
    std::array<int64_t, 6> particle_counts;
    BufferView<int64_t> particle_ids;
    std::array<int64_t, 6> primary_counts;
    int64_t size;
    char * topology;
    std::array<float, 3> vertex;
    int64_t volume_id;
    
    void SyncVectors();
    
    // note: the following 'handle' objects
    // are used internally by HDF5 to keep track
    // of the memory for variable-length buffers.
    // please use the SyncVectors() method
    // after loading data into the object
    // to fill the corresponding BufferView<>s above,
    // and then use those for access to the data.
    
    hvl_t match_handle;
    hvl_t match_overlap_handle;
    hvl_t particle_ids_handle;
  };
  
  
  struct Particle
  {
    double csda_kinetic_energy;
    double depositions_sum;
    std::array<float, 3> end_dir;
    std::array<double, 3> end_point;
    BufferView<int64_t> fragment_ids;
    int64_t id;
    int64_t image_id;
    BufferView<int64_t> index;
    int64_t interaction_id;
    bool is_contained;
    bool is_primary;
    bool is_principal_match;
    double length;
    BufferView<int64_t> match;
    BufferView<float> match_overlap;
    uint8_t matched;
    double momentum_mcs;
    int64_t nu_id;
    int64_t num_fragments;
    Pid_t pid;
    std::array<float, 5> pid_scores;
    std::array<float, 2> primary_scores;
    SemanticType_t semantic_type;
    int64_t size;
    std::array<float, 3> start_dir;
    std::array<double, 3> start_point;
    int64_t volume_id;
    
    void SyncVectors();
    
    // note: the following 'handle' objects
    // are used internally by HDF5 to keep track
    // of the memory for variable-length buffers.
    // please use the SyncVectors() method
    // after loading data into the object
    // to fill the corresponding BufferView<>s above,
    // and then use those for access to the data.
    
    hvl_t fragment_ids_handle;
    hvl_t index_handle;
    hvl_t match_handle;
    hvl_t match_overlap_handle;
  };
  
  
  struct TrueInteraction
  {
    int64_t crthit_id;
    uint8_t crthit_matched;
    int64_t crthit_matched_particle_id;
    int64_t flash_hypothesis;
    int64_t flash_id;
    double flash_time;
    int64_t flash_total_pE;
    uint8_t fmatched;
    int64_t id;
    int64_t image_id;
    bool is_contained;
    bool is_neutrino;
    bool is_principal_match;
    BufferView<int64_t> match;
    BufferView<float> match_overlap;
    uint8_t matched;
    int64_t nu_current_type;
    double nu_energy_init;
    int64_t nu_id;
    int64_t nu_interaction_mode;
    int64_t nu_interaction_type;
    int64_t num_particles;
    int64_t num_primaries;
    std::array<int64_t, 6> particle_counts;
    BufferView<int64_t> particle_ids;
    std::array<int64_t, 6> primary_counts;
    int64_t size;
    char * topology;
    BufferView<int64_t> truth_particle_counts;
    BufferView<int64_t> truth_primary_counts;
    char * truth_topology;
    std::array<double, 3> vertex;
    int64_t volume_id;
    
    void SyncVectors();
    
    // note: the following 'handle' objects
    // are used internally by HDF5 to keep track
    // of the memory for variable-length buffers.
    // please use the SyncVectors() method
    // after loading data into the object
    // to fill the corresponding BufferView<>s above,
    // and then use those for access to the data.
    
    hvl_t match_handle;
    hvl_t match_overlap_handle;
    hvl_t particle_ids_handle;
    hvl_t truth_particle_counts_handle;
    hvl_t truth_primary_counts_handle;
  };
  
  
  struct TrueParticle
  {
    char * ancestor_creation_process;
    int64_t ancestor_pdg_code;
    BufferView<double> ancestor_position;
    double ancestor_t;
    int64_t ancestor_track_id;
    BufferView<int64_t> children_counts;
    BufferView<double> children_id;
    char * creation_process;
    double csda_kinetic_energy;
    double csda_kinetic_energy_tng;
    double depositions_sum;
    BufferView<double> direction;
    double distance_travel;
    std::array<float, 3> end_dir;
    std::array<double, 3> end_point;
    std::array<double, 3> end_position;
    double energy_deposit;
    double energy_init;
    BufferView<double> first_step;
    BufferView<int64_t> fragment_ids;
    int64_t group_id;
    int64_t id;
    int64_t image_id;
    BufferView<int64_t> index;
    int64_t interaction_id;
    bool is_contained;
    bool is_primary;
    bool is_principal_match;
    BufferView<double> last_step;
    double length;
    double length_tng;
    BufferView<int64_t> match;
    BufferView<float> match_overlap;
    uint8_t matched;
    int64_t mcst_index;
    int64_t mct_index;
    std::array<double, 3> momentum;
    double momentum_mcs;
    int64_t nu_id;
    int64_t num_fragments;
    int64_t num_voxels;
    double p;
    char * parent_creation_process;
    int64_t parent_id;
    int64_t parent_pdg_code;
    BufferView<double> parent_position;
    double parent_t;
    int64_t parent_track_id;
    int64_t pdg_code;
    Pid_t pid;
    std::array<float, 5> pid_scores;
    BufferView<double> position;
    std::array<float, 2> primary_scores;
    BufferView<int64_t> sed_index;
    int64_t sed_size;
    SemanticType_t semantic_type;
    int64_t shape;
    int64_t size;
    std::array<float, 3> start_dir;
    std::array<double, 3> start_point;
    double t;
    int64_t track_id;
    BufferView<double> truth_end_dir;
    BufferView<int64_t> truth_index;
    int64_t truth_size;
    BufferView<double> truth_start_dir;
    int64_t volume_id;
    
    void SyncVectors();
    
    // note: the following 'handle' objects
    // are used internally by HDF5 to keep track
    // of the memory for variable-length buffers.
    // please use the SyncVectors() method
    // after loading data into the object
    // to fill the corresponding BufferView<>s above,
    // and then use those for access to the data.
    
    hvl_t ancestor_position_handle;
    hvl_t children_counts_handle;
    hvl_t children_id_handle;
    hvl_t direction_handle;
    hvl_t first_step_handle;
    hvl_t fragment_ids_handle;
    hvl_t index_handle;
    hvl_t last_step_handle;
    hvl_t match_handle;
    hvl_t match_overlap_handle;
    hvl_t parent_position_handle;
    hvl_t position_handle;
    hvl_t sed_index_handle;
    hvl_t truth_end_dir_handle;
    hvl_t truth_index_handle;
    hvl_t truth_start_dir_handle;
  };
  
  
  template <>
  H5::CompType BuildCompType<Event>();
  
  
  template <>
  H5::CompType BuildCompType<Interaction>();
  
  
  template <>
  H5::CompType BuildCompType<Particle>();
  
  
  template <>
  H5::CompType BuildCompType<TrueInteraction>();
  
  
  template <>
  H5::CompType BuildCompType<TrueParticle>();
  

}

#endif // CAFMAKER_TYPES_DLP_DLP_H5_CLASSES_H
