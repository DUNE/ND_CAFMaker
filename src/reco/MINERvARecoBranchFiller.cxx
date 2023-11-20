#include "MINERvARecoBranchFiller.h"

namespace cafmaker
{

  MINERvARecoBranchFiller::~MINERvARecoBranchFiller() {
  delete MnvRecoTree;
  fMnvRecoFile->Close();
  delete fMnvRecoFile;
  MnvRecoTree = NULL;
  fMnvRecoFile = NULL;
  }


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

      MnvRecoTree->SetBranchAddress("ev_run", &ev_run);
      MnvRecoTree->SetBranchAddress("ev_sub_run", &ev_sub_run);
      MnvRecoTree->SetBranchAddress("ev_gate", &ev_gate);

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


    void MINERvARecoBranchFiller::FillTrueParticle(caf::SRTrueParticle & srTruePart,
                                                 int max_trkid) const
  {
    //const auto NaN = std::numeric_limits<float>::signaling_NaN();
    ValidateOrCopy(mc_traj_pdg[max_trkid], srTruePart.pdg, 0, "pdg_code");

    ValidateOrCopy(mc_traj_edepsim_trkid[max_trkid], srTruePart.G4ID, -1, "SRTrueParticle::track_id");
    ValidateOrCopy(mc_traj_parentid[max_trkid], srTruePart.parent, -1, "SRTrueParticle::parent");

    // todo: Things do not match yet the exact Genie output, need to work on the Minerva reconstruction output. 
    // For now will assume that if it's a primary and we gor the eventID and trackid right, Genie will have fill it properly
    // And if it's an important secondary shared by both detectors, MLReco will have filled it.
    /*
    ValidateOrCopy(mc_traj_point_x[max_trkid][0]/10. - offsetX/10., srTruePart.start_pos.x, NaN, "SRTrueParticle::start_pos.x");
    ValidateOrCopy(mc_traj_point_y[max_trkid][0]/10. - offsetY/10., srTruePart.start_pos.y, NaN, "SRTrueParticle::start_pos.y");
    ValidateOrCopy(mc_traj_point_z[max_trkid][0]/10. - offsetZ/10., srTruePart.start_pos.z, NaN, "SRTrueParticle::start_pos.z");

    ValidateOrCopy(mc_traj_point_x[max_trkid][1]/10. - offsetX/10., srTruePart.end_pos.x, NaN, "SRTrueParticle::end_pos.x");
    ValidateOrCopy(mc_traj_point_y[max_trkid][1]/10. - offsetY/10., srTruePart.end_pos.y, NaN, "SRTrueParticle::end_pos.y");
    ValidateOrCopy(mc_traj_point_z[max_trkid][1]/10. - offsetZ/10., srTruePart.end_pos.z, NaN, "SRTrueParticle::end_pos.z");


    ValidateOrCopy(mc_traj_point_px[max_trkid][0], srTruePart.p.px, NaN, "SRTrueParticle::end_pos.x");
    ValidateOrCopy(mc_traj_point_py[max_trkid][0], srTruePart.p.py, NaN, "SRTrueParticle::end_pos.y");
    ValidateOrCopy(mc_traj_point_pz[max_trkid][0], srTruePart.p.pz, NaN, "SRTrueParticle::end_pos.z");

    */
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

    //Fill MINERvA specific info in the meta branch
    sr.meta.minerva.enabled = true;
    sr.meta.minerva.run = ev_run;
    sr.meta.minerva.subrun = ev_sub_run;
    sr.meta.minerva.event = ev_gate;



    int max_slice = 0;
    // Fill in the track info 
    std::map<int,std::vector<caf::SRTrack>> track_map = FillTrack(sr, max_slice, truthMatch);
    std::map<int,std::vector<caf::SRShower>> shower_map = FillShower(sr, max_slice, truthMatch);

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

  std::map<int, std::vector<caf::SRTrack>> MINERvARecoBranchFiller::FillTrack(caf::StandardRecord &sr, int & max_slice, const TruthMatcher *truthMatch) const
  {
    std::map<int,std::vector<caf::SRTrack>> track_map;

    for (int i = 0; i < n_tracks; ++i) {

      int my_slice = trk_time_slice[i];
      caf::SRTrack my_track;

      // Save first and last hit in track
      // MINERvA Reco info is saved in mm whereas CAFs use cm as default -> do conversion here
      // We offset positions in MINERvA reconstruction so we reofset them here.
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
      FindTruthTrack(sr, my_track,i, truthMatch);
      track_map[my_slice].push_back(my_track);
      //Find the largest reconstructed time slice
      if (max_slice < my_slice) max_slice = my_slice;
    }
    return track_map;
  }

  std::map<int,std::vector<caf::SRShower>> MINERvARecoBranchFiller::FillShower(caf::StandardRecord &sr, int & max_slice, const TruthMatcher *truthMatch) const
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
      FindTruthShower(sr, my_shower,i, truthMatch);
      //Fill the shower map
      shower_map[my_slice].push_back(my_shower);
      //Find the largest reconstructed time slice
      if (max_slice < my_slice) max_slice = my_slice;
    }

    return shower_map;
  }

  void MINERvARecoBranchFiller::FindTruthShower(caf::StandardRecord &sr, caf::SRShower &sh, int shower_id, const TruthMatcher *truthMatch ) const
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

    
    Int_t edepsim_track_id = mc_traj_edepsim_trkid[max_trkid];
    
    //We don't store the status of the particle (primary or not) inside MNV reco, first look in the list of primaries ID if we're around 
    std::size_t truthPartIdx = std::distance(srTrueInt.prim.begin(), std::find_if(srTrueInt.prim.begin(), srTrueInt.prim.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
    bool is_primary = truthPartIdx != srTrueInt.prim.size();

    caf::SRTrueParticle & srTruePart = is_primary ? truthMatch->GetTrueParticle(sr, srTrueInt, edepsim_track_id, true, false)
                                                    : truthMatch->GetTrueParticle(sr, srTrueInt, edepsim_track_id, false, true);

    caf::TrueParticleID truePartID;
    truePartID.ixn = truthVecIdx;
    truePartID.type = is_primary ? caf::TrueParticleID::kPrimary
                                 : caf::TrueParticleID::kSecondary;
    truePartID.part = edepsim_track_id;
    sh.truth.push_back(std::move(truePartID));

    FillTrueParticle(srTruePart, max_trkid);

  }

  void MINERvARecoBranchFiller::FindTruthTrack(caf::StandardRecord &sr, caf::SRTrack &t, int track_id, const TruthMatcher *truthMatch) const
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
    Int_t edepsim_track_id = mc_traj_edepsim_trkid[max_trkid];

    //We don't store the status of the particle (primary or not) inside MNV reco, first look in the list of primaries ID if we're around
    std::size_t truthPartIdx = std::distance(srTrueInt.prim.begin(), std::find_if(srTrueInt.prim.begin(), srTrueInt.prim.end(), [edepsim_track_id](const caf::SRTrueParticle& part) { return part.G4ID == edepsim_track_id; }));
    bool is_primary = truthPartIdx != srTrueInt.prim.size();
  
    caf::SRTrueParticle & srTruePart = is_primary ? truthMatch->GetTrueParticle(sr, srTrueInt, edepsim_track_id, true, false)
                                                    : truthMatch->GetTrueParticle(sr, srTrueInt, edepsim_track_id, false, true);


    caf::TrueParticleID truePartID;
    truePartID.ixn = truthVecIdx;
    truePartID.type = is_primary ? caf::TrueParticleID::kPrimary
                                 : caf::TrueParticleID::kSecondary;
    truePartID.part = edepsim_track_id;
    t.truth.push_back(std::move(truePartID));


    FillTrueParticle(srTruePart,max_trkid);

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
      unsigned long int t0_minerva;
      for (int entry = 0; entry < MnvRecoTree->GetEntries(); entry++)
      {


        MnvRecoTree->GetEntry(entry);
        if (entry == 0) t0_minerva = ev_gps_time_sec;


        fTriggers.emplace_back();
        Trigger & trig = fTriggers.back();

        trig.evtID = Long_t(ev_gl_gate);


        // todo: these are placeholder values until we can propagate enough info through the reco files
        trig.triggerType = ev_trigger_type;
        trig.triggerTime_s = ev_gps_time_sec;
        trig.triggerTime_ns = ev_gps_time_usec * 1000.;


        //Initialize first trigger at 0 while we don't have a global time strategy
        trig.triggerTime_s -= t0_minerva;


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
