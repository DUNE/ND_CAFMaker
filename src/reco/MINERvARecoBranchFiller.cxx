#include "MINERvARecoBranchFiller.h"

namespace cafmaker
{
  MINERvARecoBranchFiller::MINERvARecoBranchFiller(const std::string &minervaRecoFilename)
  : IRecoBranchFiller("MINERvA"),
    fTriggers(),
    fLastTriggerReqd(fTriggers.end())
  {
    fMnvRecoFile = new TFile(minervaRecoFilename.c_str(), "READ");
    if (!fMnvRecoFile->IsZombie()) {
      SetConfigured(true);
      // Save pointer to input tree
      MnvRecoTree = dynamic_cast<TTree*>(fMnvRecoFile->Get("minerva"));
      if (!MnvRecoTree) {
        std::cerr << "Did not find MINERvA reco tree in input file " << minervaRecoFilename << std::endl;
        std::cerr << "Are you sure it exists?" << std::endl;

        throw;
      }

      //Meta branches
      MnvRecoTree->SetBranchAddress("offsetX", &offsetX);
      MnvRecoTree->SetBranchAddress("offsetY", &offsetY);
      MnvRecoTree->SetBranchAddress("offsetZ", &offsetZ);

      MnvRecoTree->SetBranchAddress("ev_trigger_type", &ev_trigger_type);
      MnvRecoTree->SetBranchAddress("ev_gl_gate", &ev_gl_gate);
      MnvRecoTree->SetBranchAddress("ev_gps_time_sec", &ev_gps_time_sec);
      MnvRecoTree->SetBranchAddress("ev_gps_time_usec", &ev_gps_time_usec);

      //tracks branches 
      MnvRecoTree->SetBranchAddress("n_tracks", &n_tracks);
      MnvRecoTree->SetBranchAddress("trk_index", trk_index);
      MnvRecoTree->SetBranchAddress("trk_type", trk_type);
      MnvRecoTree->SetBranchAddress("trk_patrec", trk_patrec);
      MnvRecoTree->SetBranchAddress("trk_time_slice", trk_time_slice);
      MnvRecoTree->SetBranchAddress("trk_vis_energy", trk_vis_energy);
      MnvRecoTree->SetBranchAddress("trk_theta", trk_theta);
      MnvRecoTree->SetBranchAddress("trk_phi", trk_phi);
      MnvRecoTree->SetBranchAddress("trk_hits", trk_hits);
      MnvRecoTree->SetBranchAddress("trk_dof", trk_dof);
      MnvRecoTree->SetBranchAddress("trk_chi2perDof", trk_chi2perDof);
      MnvRecoTree->SetBranchAddress("trk_fitMass", trk_fitMass);
      MnvRecoTree->SetBranchAddress("trk_nodes", trk_nodes);
      MnvRecoTree->SetBranchAddress("trk_node_X", trk_node_X);
      MnvRecoTree->SetBranchAddress("trk_node_Y", trk_node_Y);
      MnvRecoTree->SetBranchAddress("trk_node_Z", trk_node_Z);
      MnvRecoTree->SetBranchAddress("trk_node_aX", trk_node_aX);
      MnvRecoTree->SetBranchAddress("trk_node_aY", trk_node_aY);
      MnvRecoTree->SetBranchAddress("trk_node_qOverP", trk_node_qOverP);
      MnvRecoTree->SetBranchAddress("trk_node_chi2", trk_node_chi2);
      MnvRecoTree->SetBranchAddress("trk_node_cluster_idx", trk_node_cluster_idx);


      //Shower branches
      MnvRecoTree->SetBranchAddress("n_blobs_id", &n_blobs_id);
      MnvRecoTree->SetBranchAddress("blob_id_idx", blob_id_idx);
      MnvRecoTree->SetBranchAddress("blob_id_subdet", blob_id_subdet);
      MnvRecoTree->SetBranchAddress("blob_id_history", blob_id_history);
      MnvRecoTree->SetBranchAddress("blob_id_size", blob_id_size);
      MnvRecoTree->SetBranchAddress("blob_id_patrec", blob_id_patrec);
      MnvRecoTree->SetBranchAddress("blob_id_e", blob_id_e);
      MnvRecoTree->SetBranchAddress("blob_id_time", blob_id_time);
      MnvRecoTree->SetBranchAddress("blob_id_time_slice", blob_id_time_slice);
      MnvRecoTree->SetBranchAddress("blob_id_startpoint_x", blob_id_startpoint_x);
      MnvRecoTree->SetBranchAddress("blob_id_startpoint_y", blob_id_startpoint_y);
      MnvRecoTree->SetBranchAddress("blob_id_startpoint_z", blob_id_startpoint_z);
      MnvRecoTree->SetBranchAddress("blob_id_clus_idx", blob_id_clus_idx);


      //Truth branches
      MnvRecoTree->SetBranchAddress("n_mc_trajectories", &n_mc_trajectories);
      MnvRecoTree->SetBranchAddress("mc_traj_trkid", mc_traj_trkid);
      MnvRecoTree->SetBranchAddress("mc_traj_parentid", mc_traj_parentid);
      MnvRecoTree->SetBranchAddress("mc_traj_pdg", mc_traj_pdg);
      MnvRecoTree->SetBranchAddress("mc_traj_hit_e", mc_traj_hit_e);
      MnvRecoTree->SetBranchAddress("mc_traj_npoints", mc_traj_npoints);
      MnvRecoTree->SetBranchAddress("mc_traj_point_x", mc_traj_point_x);
      MnvRecoTree->SetBranchAddress("mc_traj_point_y", mc_traj_point_y);
      MnvRecoTree->SetBranchAddress("mc_traj_point_z", mc_traj_point_z);
      MnvRecoTree->SetBranchAddress("mc_traj_point_t", mc_traj_point_t);
      MnvRecoTree->SetBranchAddress("mc_traj_point_px", mc_traj_point_px);
      MnvRecoTree->SetBranchAddress("mc_traj_point_py", mc_traj_point_py);
      MnvRecoTree->SetBranchAddress("mc_traj_point_pz", mc_traj_point_pz);
      MnvRecoTree->SetBranchAddress("mc_traj_point_E", mc_traj_point_E);
      MnvRecoTree->SetBranchAddress("mc_traj_edepsim_eventid", mc_traj_edepsim_eventid);
      MnvRecoTree->SetBranchAddress("mc_traj_edepsim_trkid", mc_traj_edepsim_trkid);

      //Truth Matching variables
      MnvRecoTree->SetBranchAddress("clus_id_size", clus_id_size);
      MnvRecoTree->SetBranchAddress("clus_id_hits_idx", clus_id_hits_idx);
      MnvRecoTree->SetBranchAddress("mc_id_mchit_trkid", mc_id_mchit_trkid);
      MnvRecoTree->SetBranchAddress("mc_id_mchit_dE", mc_id_mchit_dE);

      


    } else {
      fMnvRecoFile = NULL;
      MnvRecoTree = NULL;
      std::cerr << "Did not find input MINERvA reco file you provided: " << minervaRecoFilename << std::endl;
      std::cerr << "Are you sure it exists?" << std::endl;
      throw;
    }


  }
  // ---------------------------------------------------------------------------
//Workaround to just use MINERvA passthrough true information
  void MINERvARecoBranchFiller::fillTrueInteraction(caf::StandardRecord &sr) const
  {
    std::map<long int, std::vector<caf::SRTrueParticle>> particle_map;
    for (int i = 0; i<n_mc_trajectories; i++)
    {
      //Get the true particles in the event;
      caf::SRTrueParticle part;
      part.pdg = mc_traj_pdg[i];
      part.G4ID = mc_traj_edepsim_trkid[i]; // Track id of Geant. reset for every interaction
      part.interaction_id = mc_traj_edepsim_eventid[i]; // Vector ID of edepsim Unique number for all interactionss
      part.time = mc_traj_pdg[i];
      
      part.time = mc_traj_point_t[i][0];
      part.parent = mc_traj_parentid[i];

      part.start_pos = caf::SRVector3D(mc_traj_point_x[i][0],mc_traj_point_y[i][0],mc_traj_point_z[i][0]);
      part.end_pos = caf::SRVector3D(mc_traj_point_x[i][1],mc_traj_point_y[i][1],mc_traj_point_z[i][1]);
      particle_map[mc_traj_edepsim_eventid[i]].push_back(part);
    }
    for (auto const& it : particle_map)
    {
      caf::SRTrueInteraction * ixn = nullptr;
      //Check if the id already exist in the interaction vector.
      Long_t neutrino_event_id = it.first;

      auto int_it = std::find_if(sr.mc.nu.begin(), sr.mc.nu.end(),
                  [neutrino_event_id](const caf::SRTrueInteraction& nu)
                  {
                    return nu.id == neutrino_event_id;
                  });
      if (int_it == sr.mc.nu.end()) // No matching neutrino in the record already
      {
        sr.mc.nu.emplace_back();
        sr.mc.nnu++;

        ixn = &sr.mc.nu.back();
        ixn->id = it.first;
        for (std::size_t k = 0; k< it.second.size(); k++)
        {

          ixn->sec.push_back(std::move(it.second[k]));
          ixn->nsec ++;
        }
      }
      else // Found a matching neutrino 
      {
        //Look in the primaries if the true particle was already stored by other detectors
        for (std::size_t k = 0; k< it.second.size(); k++)
        {
          Int_t edepsim_track_id = it.second[k].G4ID;
          
          auto part_it = std::find_if((*int_it).prim.begin(), (*int_it).prim.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; });
          if (part_it != (*int_it).prim.end()) break; // Particle already filled so no need to fill it
          
          part_it = std::find_if((*int_it).prefsi.begin(), (*int_it).prefsi.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; });
          if (part_it != (*int_it).prefsi.end()) break; // Particle already filled so no need to fill it
          
          part_it = std::find_if((*int_it).sec.begin(), (*int_it).sec.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; });
          if (part_it != (*int_it).sec.end()) break; // Particle already filled so no need to fill it
          
          //Not found the true particle in the list of true particle, I'm adding it as a secondary particle (primary particle should have been filled by LAr CAFs)
          (*int_it).sec.push_back(std::move(it.second[k]));
          (*int_it).nsec ++;              
        }
      }  
    }
  }

  // here we copy all the MINERvA reco into the SRMINERvA branch of the StandardRecord object.
  void MINERvARecoBranchFiller::_FillRecoBranches(const Trigger &trigger,
                                                caf::StandardRecord &sr,
                                                const cafmaker::Params &par,
                                                const TruthMatcher *truthMatch) const
  {

    // figure out where in our list of triggers this event index is.
    // we should always be looking forwards, since we expect to be traversing in that direction
    auto it_start = (fLastTriggerReqd == fTriggers.end()) ? fTriggers.cbegin() : fLastTriggerReqd;
    auto itTrig = std::find(it_start, fTriggers.cend(), trigger);
    if (itTrig == fTriggers.end())
    {
      LOG.FATAL() << "Reco branch filler '" << GetName() << "' could not find trigger with evtID == " << trigger.evtID << "!  Abort.\n";
      abort();
    }
    std::size_t idx = std::distance(fTriggers.cbegin(), itTrig);

    LOG.VERBOSE() << "    Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID << ", internal evt idx = " << idx << ".\n";


    // Get nth entry from tree
    MnvRecoTree->GetEntry(idx);  
    //Hack for the true interaction
    // fillTrueInteraction(sr);


    int max_slice = 0;
    // Fill in the track info 
    std::map<int,std::vector<caf::SRTrack>> track_map = fill_track(sr, max_slice, truthMatch);
    std::map<int,std::vector<caf::SRShower>> shower_map = fill_shower(sr, max_slice, truthMatch);

    for (int i_slice = 1; i_slice <= max_slice; i_slice++)
    {
      if (track_map[i_slice].size() == 0 && shower_map[i_slice].size() == 0) continue;
      caf::SRMINERvAInt my_int;

      my_int.ntracks = track_map[i_slice].size();
      my_int.nshowers = shower_map[i_slice].size();

      for (std::size_t i_track = 0; i_track < track_map[i_slice].size(); i_track ++) my_int.tracks.push_back(std::move(track_map[i_slice][i_track]));

      for (std::size_t i_shower = 0; i_shower < shower_map[i_slice].size(); i_shower ++) my_int.showers.push_back(std::move(shower_map[i_slice][i_shower]));
      
      sr.nd.minerva.ixn.push_back(my_int);
      sr.nd.minerva.nixn +=1;
    }

  }

  std::map<int, std::vector<caf::SRTrack>> MINERvARecoBranchFiller::fill_track(caf::StandardRecord &sr, int & max_slice, const TruthMatcher *truthMatch) const
  {
    std::map<int,std::vector<caf::SRTrack>> track_map;

    for (int i = 0; i < n_tracks; ++i) {

      int my_slice = trk_time_slice[i];
      caf::SRTrack my_track;

      // Save first and last hit in track
      // MINERvA Reco info is saved in mm whereas CAFs use CM as default -> do conversion here
      my_track.start = caf::SRVector3D(trk_node_X[i][0]/10.  - offsetX/10. ,trk_node_Y[i][0]/10. - offsetY/10., trk_node_Z[i][0]/10. - offsetZ/10.);
      my_track.end   = caf::SRVector3D(trk_node_X[i][trk_nodes[i] -1]/10.  - offsetX/10. ,trk_node_Y[i][trk_nodes[i] -1]/10. - offsetY/10., trk_node_Z[i][trk_nodes[i] -1]/10. - offsetZ/10.);

      // Track info
      my_track.len_cm  = sqrt(pow(trk_node_X[i][trk_nodes[i] -1] - trk_node_X[i][0],2)+pow(trk_node_Y[i][trk_nodes[i] -1] - trk_node_Y[i][0],2)+pow(trk_node_Z[i][trk_nodes[i] -1] - trk_node_Z[i][0],2))/10.;
      my_track.qual      = trk_chi2perDof[i];
      //We don't have different track energy definition so we fill both E and Evis for now in case only one of the two is used
      my_track.Evis    = trk_vis_energy[i]/1000.; //In GeV
      my_track.E       = trk_vis_energy[i]/1000.; //In GeV


      // Get the directions
      my_track.dir     = caf::SRVector3D(sin(trk_theta[i])*cos(trk_phi[i]),sin(trk_theta[i])*sin(trk_phi[i]),cos(trk_theta[i]));
      my_track.enddir  = caf::SRVector3D(sin(trk_theta[i])*cos(trk_phi[i]),sin(trk_theta[i])*sin(trk_phi[i]),cos(trk_theta[i]));

      //Associates the truth particle to the track
      find_truth_track(sr, my_track,i, truthMatch);
      track_map[my_slice].push_back(my_track);
      //Find the largest reconstructed time slice
      if (max_slice < my_slice) max_slice = my_slice;
    }
    return track_map;
  }

  std::map<int,std::vector<caf::SRShower>> MINERvARecoBranchFiller::fill_shower(caf::StandardRecord &sr, int & max_slice, const TruthMatcher *truthMatch) const
  {
    std::map<int,std::vector<caf::SRShower>> shower_map;

    for (int i = 0; i < n_blobs_id; ++i) 
    {
      int my_slice = blob_id_time_slice[i];
      caf::SRShower my_shower;

      // Save first and last hit in track
      // MINERvA Reco info is saved in mm whereas CAFs use CM as default -> do conversion here
      my_shower.start = caf::SRVector3D(blob_id_startpoint_x[i]/10. - offsetX/10.,blob_id_startpoint_y[i]/10. - offsetY/10., blob_id_startpoint_z[i]/10. - offsetZ/10.);
      

      //Actual direction but Centroid makes more sense 
//      double x_dir = (blob_id_centroid_x[i] - blob_id_startpoint_x[i]);
//      double y_dir = (blob_id_centroid_y[i] - blob_id_startpoint_y[i]);
//      double z_dir = (blob_id_centroid_z[i] - blob_id_startpoint_z[i]);
      
//      double norm_dir = sqrt(x_dir*x_dir + y_dir*y_dir + z_dir*z_dir);
//      x_dir /= norm_dir;
//      y_dir /= norm_dir;
//      z_dir /= norm_dir;


      //Use centroid as blob direction
      double x_dir = blob_id_centroid_x[i]/10. - offsetX/10.;
      double y_dir = blob_id_centroid_y[i]/10. - offsetY/10.;
      double z_dir = blob_id_centroid_z[i]/10. - offsetZ/10.;

      my_shower.direction = caf::SRVector3D(x_dir, y_dir, z_dir);
      my_shower.Evis = blob_id_e[i]/1000.; //Energy in GeV

      //Associates the truth particle to the shower
      find_truth_shower(sr, my_shower,i, truthMatch);
      //Fill the shower map
      shower_map[my_slice].push_back(my_shower);
      //Find the largest reconstructed time slice
      if (max_slice < my_slice) max_slice = my_slice;
    }

    return shower_map;
  }

  void MINERvARecoBranchFiller::find_truth_shower(caf::StandardRecord &sr, caf::SRShower &sh, int shower_id, const TruthMatcher *truthMatch ) const
  {
    std::map<int, double> most_trkid;
    for (int j = 0; j<blob_id_size[shower_id]; j++)
    {
        //Get the cluster associated to the node
        int id_cl = blob_id_clus_idx[shower_id][j];
        
        //Get the # of digits (hits) in the cluster
        int clus_size = clus_id_size[id_cl];

        for (int k = 0; k<clus_size; k++)
        {
            int id_hit = clus_id_hits_idx[id_cl][0];
            int traj_id = mc_id_mchit_trkid[id_hit][0]; // Get the true trajectory associated to the hit
            if (mc_id_mchit_dE[id_hit][0] > 0) 
            {
              most_trkid[traj_id] += mc_id_mchit_dE[id_hit][0]; // Looks for trajectory that contributed the most to the track 
            }
        }          
    }
   
    double max_trkid_stat = 0;
    int max_trkid = 0;

    for (auto const& tkid : most_trkid) 
    {
        if (tkid.second > max_trkid_stat) 
        {
            max_trkid_stat = tkid.second;
            max_trkid = tkid.first; 
        }
    }

    //Create the possible SRTrueParticle that correspond to the shower
    caf::SRTrueParticle part;
    part.pdg = mc_traj_pdg[max_trkid];
    part.G4ID = mc_traj_edepsim_trkid[max_trkid]; // Track id of Geant. reset for every interaction
    part.interaction_id = mc_traj_edepsim_eventid[max_trkid]; // Vector ID of edepsim Unique number for all interactionss
    part.time = mc_traj_pdg[max_trkid];
    
    part.time = mc_traj_point_t[max_trkid][0];
    part.parent = mc_traj_parentid[max_trkid];

    part.start_pos = caf::SRVector3D(mc_traj_point_x[max_trkid][0],mc_traj_point_y[max_trkid][0],mc_traj_point_z[max_trkid][0]);
    part.end_pos = caf::SRVector3D(mc_traj_point_x[max_trkid][1],mc_traj_point_y[max_trkid][1],mc_traj_point_z[max_trkid][1]);

    
    Long_t neutrino_event_id = mc_traj_edepsim_eventid[max_trkid];
    std::size_t truthVecIdx = std::distance(sr.mc.nu.begin(),
                                                  std::find_if(sr.mc.nu.begin(),
                                                               sr.mc.nu.end(),
                                                               [neutrino_event_id](const caf::SRTrueInteraction& nu)
                                                               {
                                                                 return nu.id == neutrino_event_id;
                                                               }));
//Once the true particle has been found, loop over the truthBranch to find the corresponding truth interraction
    caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, neutrino_event_id, true);
    //Find the position of the interaction corresponding to the track in the interaction vector

    sh.truth.ixn = truthVecIdx;
    Int_t edepsim_track_id = mc_traj_edepsim_trkid[max_trkid];

    //Look in the primaries if the particle wasn't filled already
    std::size_t truthPartIdx = std::distance(srTrueInt.prim.begin(), std::find_if(srTrueInt.prim.begin(), srTrueInt.prim.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
    if (truthPartIdx != srTrueInt.prim.size()) 
    {
      sh.truth.type = caf::TrueParticleID::kPrimary;
      sh.truth.part = truthPartIdx;
    }
    else {
      
      //Pre fsi trajectory?
      truthPartIdx = std::distance(srTrueInt.prefsi.begin(), std::find_if(srTrueInt.prefsi.begin(), srTrueInt.prefsi.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
      if (truthPartIdx != srTrueInt.prefsi.size())
      {
        sh.truth.type = caf::TrueParticleID::kPrimaryBeforeFSI;
        sh.truth.part = truthPartIdx;
      }
      else 
      {
        //If nothing else it's a secondary
        truthPartIdx = std::distance(srTrueInt.sec.begin(), std::find_if(srTrueInt.sec.begin(), srTrueInt.sec.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
        if (truthPartIdx != srTrueInt.sec.size())
        {
          sh.truth.type = caf::TrueParticleID::kSecondary;
          sh.truth.part = truthPartIdx;
        } 
        
        else
        {
          //If not found before, add in in the secondary
          srTrueInt.sec.push_back(std::move(part));
          srTrueInt.nsec ++;
          

          sh.truth.ixn = truthVecIdx;
          sh.truth.part = srTrueInt.nsec -1;
          sh.truth.type = caf::TrueParticleID::kSecondary;
        }
      }
    }


    

    // if (truthVecIdx == sr.mc.nu.size())
    // {
    //   caf::SRTrueInteraction * ixn = nullptr;

    //   sr.mc.nu.emplace_back();
    //   sr.mc.nnu++;
    //   //NOE newly added
    //   ixn = &sr.mc.nu.back();
    //   ixn->id = neutrino_event_id;
    //   ixn->sec.push_back(std::move(part));
    //   ixn->nsec ++;

    //   sh.truth.ixn = truthVecIdx;
    //   sh.truth.part = ixn->nsec -1;
    //   sh.truth.type = caf::TrueParticleID::kSecondary;
    // }
    // else
    // {
    //   sh.truth.ixn = truthVecIdx;
    //   Int_t edepsim_track_id = mc_traj_edepsim_trkid[max_trkid];
      
    //   std::size_t truthPartIdx = std::distance(sr.mc.nu[truthVecIdx].prim.begin(), std::find_if(sr.mc.nu[truthVecIdx].prim.begin(), sr.mc.nu[truthVecIdx].prim.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
    //   if (truthPartIdx != sr.mc.nu[truthVecIdx].prim.size()) 
    //   {
    //     sh.truth.type = caf::TrueParticleID::kPrimary;
    //     sh.truth.part = truthPartIdx;
    //   }
    //   else {
    //     truthPartIdx = std::distance(sr.mc.nu[truthVecIdx].prefsi.begin(), std::find_if(sr.mc.nu[truthVecIdx].prefsi.begin(), sr.mc.nu[truthVecIdx].prefsi.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
    //     if (truthPartIdx != sr.mc.nu[truthVecIdx].prefsi.size())
    //     {
    //       sh.truth.type = caf::TrueParticleID::kPrimaryBeforeFSI;
    //       sh.truth.part = truthPartIdx;
    //     }
    //     else 
    //     {
    //       truthPartIdx = std::distance(sr.mc.nu[truthVecIdx].sec.begin(), std::find_if(sr.mc.nu[truthVecIdx].sec.begin(), sr.mc.nu[truthVecIdx].sec.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
    //       if (truthPartIdx != sr.mc.nu[truthVecIdx].sec.size())
    //       {
    //         sh.truth.type = caf::TrueParticleID::kSecondary;
    //         sh.truth.part = truthPartIdx;
    //       } 
          
    //       else
    //       {
    //         sr.mc.nu[truthVecIdx].sec.push_back(std::move(part));
    //         sr.mc.nu[truthVecIdx].nsec ++;
            

    //         sh.truth.ixn = truthVecIdx;
    //         sh.truth.part = sr.mc.nu[truthVecIdx].nsec -1;
    //         sh.truth.type = caf::TrueParticleID::kSecondary;
    //       }
    //     }
    //   }
    // }

  }

  void MINERvARecoBranchFiller::find_truth_track(caf::StandardRecord &sr, caf::SRTrack &t, int track_id, const TruthMatcher *truthMatch) const
  {
    std::map<int, double> most_trkid;
    for (int j = 0; j<trk_nodes[track_id]; j++)
    {
        //Get the cluster associated to the node
        int id_cl = trk_node_cluster_idx[track_id][j];
        
        //Get the # of digits (hits) in the cluster
        int clus_size = clus_id_size[id_cl];

        for (int k = 0; k<clus_size; k++)
        {
            int id_hit = clus_id_hits_idx[id_cl][0];
            int traj_id = mc_id_mchit_trkid[id_hit][0]; // Get the true trajectory associated to the hit
            if (mc_id_mchit_dE[id_hit][0] > 0) 
            {
              most_trkid[traj_id] += mc_id_mchit_dE[id_hit][0]; // Looks for trajectory that contributed the most to the track 
            }
        }          
    }
   
    double max_trkid_stat = 0;
    int max_trkid = 0;

    

    for (auto const& tkid : most_trkid) 
    {
        if (tkid.second > max_trkid_stat) 
        {
            max_trkid_stat = tkid.second;
            max_trkid = tkid.first; 
        }
    }

    //find the true particle associated to this track

    caf::SRTrueParticle part;
    part.pdg = mc_traj_pdg[max_trkid];
    part.G4ID = mc_traj_edepsim_trkid[max_trkid]; // Track id of Geant. reset for every interaction
    part.interaction_id = mc_traj_edepsim_eventid[max_trkid]; // Vector ID of edepsim Unique number for all interactionss
    
    part.time = mc_traj_point_t[max_trkid][0];
    part.parent = mc_traj_parentid[max_trkid];

    part.start_pos = caf::SRVector3D(mc_traj_point_x[max_trkid][0],mc_traj_point_y[max_trkid][0],mc_traj_point_z[max_trkid][0]);
    part.end_pos = caf::SRVector3D(mc_traj_point_x[max_trkid][1],mc_traj_point_y[max_trkid][1],mc_traj_point_z[max_trkid][1]);


    //Once the true particle has been found, loop over the truthBranch to find the corresponding truth interraction
    Long_t neutrino_event_id = mc_traj_edepsim_eventid[max_trkid];
    caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, neutrino_event_id, true);
    //Find the position of the interaction corresponding to the track in the interaction vector



    std::size_t truthVecIdx = std::distance(sr.mc.nu.begin(),
                                                  std::find_if(sr.mc.nu.begin(),
                                                               sr.mc.nu.end(),
                                                               [neutrino_event_id](const caf::SRTrueInteraction& nu)
                                                               {
                                                                 return nu.id == neutrino_event_id;
                                                               }));
    t.truth.ixn = truthVecIdx;
    Int_t edepsim_track_id = mc_traj_edepsim_trkid[max_trkid];

    //Look in the primaries if the particle wasn't filled already
    std::size_t truthPartIdx = std::distance(srTrueInt.prim.begin(), std::find_if(srTrueInt.prim.begin(), srTrueInt.prim.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
    if (truthPartIdx != srTrueInt.prim.size()) 
    {
      t.truth.type = caf::TrueParticleID::kPrimary;
      t.truth.part = truthPartIdx;
    }
    else {
      
      //Pre fsi trajectory?
      truthPartIdx = std::distance(srTrueInt.prefsi.begin(), std::find_if(srTrueInt.prefsi.begin(), srTrueInt.prefsi.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
      if (truthPartIdx != srTrueInt.prefsi.size())
      {
        t.truth.type = caf::TrueParticleID::kPrimaryBeforeFSI;
        t.truth.part = truthPartIdx;
      }
      else 
      {
        //If nothing else it's a secondary
        truthPartIdx = std::distance(srTrueInt.sec.begin(), std::find_if(srTrueInt.sec.begin(), srTrueInt.sec.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
        if (truthPartIdx != srTrueInt.sec.size())
        {
          t.truth.type = caf::TrueParticleID::kSecondary;
          t.truth.part = truthPartIdx;
        } 
        
        else
        {
          //If not found before, add in in the secondary
          srTrueInt.sec.push_back(std::move(part));
          srTrueInt.nsec ++;
          

          t.truth.ixn = truthVecIdx;
          t.truth.part = srTrueInt.nsec -1;
          t.truth.type = caf::TrueParticleID::kSecondary;
        }
      }
    }


    //Hack before the new GENIE event ID
    /*
    std::size_t truthVecIdx = std::distance(sr.mc.nu.begin(),
                                                  std::find_if(sr.mc.nu.begin(),
                                                               sr.mc.nu.end(),
                                                               [neutrino_event_id](const caf::SRTrueInteraction& nu)
                                                               {
                                                                 return nu.id == neutrino_event_id;
                                                               }));

    if (truthVecIdx == sr.mc.nu.size())
    {
      caf::SRTrueInteraction * ixn = nullptr;
      sr.mc.nu.emplace_back();
      sr.mc.nnu++;
      ixn = &sr.mc.nu.back();
      ixn->id = neutrino_event_id;
      ixn->sec.push_back(std::move(part));
      ixn->nsec ++;

      sh.truth.ixn = truthVecIdx;
      sh.truth.part = ixn->nsec -1;
      sh.truth.type = caf::TrueParticleID::kSecondary;
    }
    else
    {
      sh.truth.ixn = truthVecIdx;
      Int_t edepsim_track_id = mc_traj_edepsim_trkid[max_trkid];
      
      std::size_t truthPartIdx = std::distance(sr.mc.nu[truthVecIdx].prim.begin(), std::find_if(sr.mc.nu[truthVecIdx].prim.begin(), sr.mc.nu[truthVecIdx].prim.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
      if (truthPartIdx != sr.mc.nu[truthVecIdx].prim.size()) 
      {
        sh.truth.type = caf::TrueParticleID::kPrimary;
        sh.truth.part = truthPartIdx;
      }
      else {
        truthPartIdx = std::distance(sr.mc.nu[truthVecIdx].prefsi.begin(), std::find_if(sr.mc.nu[truthVecIdx].prefsi.begin(), sr.mc.nu[truthVecIdx].prefsi.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
        if (truthPartIdx != sr.mc.nu[truthVecIdx].prefsi.size())
        {
          sh.truth.type = caf::TrueParticleID::kPrimaryBeforeFSI;
          sh.truth.part = truthPartIdx;
        }
        else 
        {
          truthPartIdx = std::distance(sr.mc.nu[truthVecIdx].sec.begin(), std::find_if(sr.mc.nu[truthVecIdx].sec.begin(), sr.mc.nu[truthVecIdx].sec.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
          if (truthPartIdx != sr.mc.nu[truthVecIdx].sec.size())
          {
            sh.truth.type = caf::TrueParticleID::kSecondary;
            sh.truth.part = truthPartIdx;
          } 
          
          else
          {
            sr.mc.nu[truthVecIdx].sec.push_back(std::move(part));
            sr.mc.nu[truthVecIdx].nsec ++;
            

            sh.truth.ixn = truthVecIdx;
            sh.truth.part = sr.mc.nu[truthVecIdx].nsec -1;
            sh.truth.type = caf::TrueParticleID::kSecondary;
          }
        }
      }
    }*/
  }

  // ------------------------------------------------------------------------------
  std::deque<Trigger> MINERvARecoBranchFiller::GetTriggers(int triggerType) const
  {
    //Just read the branches of interest

    std::deque<Trigger> triggers;
    if (fTriggers.empty())
    {
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "' from " << MnvRecoTree->GetEntries() << " MINERvA Tree:\n";
      fTriggers.reserve(MnvRecoTree->GetEntries());
      for (int entry = 0; entry < MnvRecoTree->GetEntries(); entry++)
      {


        MnvRecoTree->GetEntry(entry);



        fTriggers.emplace_back();
        Trigger & trig = fTriggers.back();

        trig.evtID = Long_t(ev_gl_gate);


        // todo: these are placeholder values until we can propagate enough info through the reco files
        trig.triggerType = ev_trigger_type;
        trig.triggerTime_s = ev_gps_time_sec;
        trig.triggerTime_ns = ev_gps_time_usec * 1000.;


        triggers.push_back(trig);

        LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
                      << ", triggerType=" << trig.triggerType
                      << ", triggerTime_s=" << trig.triggerTime_s
                      << ", triggerTime_ns=" << trig.triggerTime_ns
                      << "\n";

      }
      fLastTriggerReqd = fTriggers.end();  // since we just modified the list, any iterators have been invalidated
    }
    else
    {
      for (const Trigger & trigger : fTriggers)
      {
        if (triggerType < 0 || triggerType == fTriggers.back().triggerType)
          triggers.push_back(trigger);
      }
    }
    return triggers;
  }

} // end namespace
