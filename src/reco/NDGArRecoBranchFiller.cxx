/// Fill NDGAr reco branches using NDGAr reco data
///
/// \author  F. Martinez Lopez
/// \date    Oct. 2022
///

#include "NDGArRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

namespace cafmaker
{
  
  NDGArRecoBranchFiller::NDGArRecoBranchFiller(const std::string &NDGArRecoFilename)
  {  
    fNDGArRecoFile = new TFile(NDGArRecoFilename.c_str(), "READ");
    if(!fNDGArRecoFile->IsZombie()){
      SetConfigured(true);

      // fTree = fNDGArRecoFile->Get<TTree>("anatree/GArAnaTree");
      NDGArRecoTree = dynamic_cast<TTree*>(fNDGArRecoFile->Get("anatree/GArAnaTree"));

      NDGArRecoTree->SetBranchAddress("Event", &fEvent);

      NDGArRecoTree->SetBranchAddress("NType", &fMCNuPDG);
      NDGArRecoTree->SetBranchAddress("MCNuPx", &fMCNuPX);
      NDGArRecoTree->SetBranchAddress("MCNuPy", &fMCNuPY);
      NDGArRecoTree->SetBranchAddress("MCNuPz", &fMCNuPZ);

      ///MC Particles info
      NDGArRecoTree->SetBranchAddress("MCTrkID", &fMCTrkID);
      NDGArRecoTree->SetBranchAddress("PDG", &fMCPDG);
      NDGArRecoTree->SetBranchAddress("MotherIndex", &fMCMotherIndex);
      NDGArRecoTree->SetBranchAddress("PDGMother", &fMCPDGMother);

      NDGArRecoTree->SetBranchAddress("MCPStartX", &fMCStartX);
      NDGArRecoTree->SetBranchAddress("MCPStartY", &fMCStartY);
      NDGArRecoTree->SetBranchAddress("MCPStartZ", &fMCStartZ);
      NDGArRecoTree->SetBranchAddress("MCPStartPX", &fMCStartPX);
      NDGArRecoTree->SetBranchAddress("MCPStartPY", &fMCStartPY);
      NDGArRecoTree->SetBranchAddress("MCPStartPZ", &fMCStartPZ);

      //Track-related info
      NDGArRecoTree->SetBranchAddress("TrackStartX", &fTrackStartX);
      NDGArRecoTree->SetBranchAddress("TrackStartY", &fTrackStartY);
      NDGArRecoTree->SetBranchAddress("TrackStartZ", &fTrackStartZ);
      NDGArRecoTree->SetBranchAddress("TrackEndX", &fTrackEndX);
      NDGArRecoTree->SetBranchAddress("TrackEndY", &fTrackEndY);
      NDGArRecoTree->SetBranchAddress("TrackEndZ", &fTrackEndZ);

      NDGArRecoTree->SetBranchAddress("TrackStartPX", &fTrackStartPX);
      NDGArRecoTree->SetBranchAddress("TrackStartPY", &fTrackStartPY);
      NDGArRecoTree->SetBranchAddress("TrackStartPZ", &fTrackStartPZ);
      NDGArRecoTree->SetBranchAddress("TrackEndPX", &fTrackEndPX);
      NDGArRecoTree->SetBranchAddress("TrackEndPY", &fTrackEndPY);
      NDGArRecoTree->SetBranchAddress("TrackEndPZ", &fTrackEndPZ);

      NDGArRecoTree->SetBranchAddress("TrackLenF", &fTrackLenF);
      NDGArRecoTree->SetBranchAddress("TrackLenB", &fTrackLenB);
      NDGArRecoTree->SetBranchAddress("TrackPF", &fTrackPF);
      NDGArRecoTree->SetBranchAddress("TrackPB", &fTrackPB);
      NDGArRecoTree->SetBranchAddress("TrackAvgIonF", &fTrackAvgIonF);
      NDGArRecoTree->SetBranchAddress("TrackAvgIonB", &fTrackAvgIonB);

      NDGArRecoTree->SetBranchAddress("TrackIDNumber", &fTrackIDNumber);
      NDGArRecoTree->SetBranchAddress("NTPCClustersOnTrack", &fTrackNClusters);

      NDGArRecoTree->SetBranchAddress("TrackPIDF", &fTrackPIDF);
      NDGArRecoTree->SetBranchAddress("TrackPIDProbF", &fTrackPIDProbF);
      NDGArRecoTree->SetBranchAddress("TrackPIDB", &fTrackPIDB);
      NDGArRecoTree->SetBranchAddress("TrackPIDProbB", &fTrackPIDProbB);

      NDGArRecoTree->SetBranchAddress("TrackMCindex", &fTrackMCindex);
      NDGArRecoTree->SetBranchAddress("TrackMCfrac", &fTrackMCfrac);

      //ECAL-related info
      NDGArRecoTree->SetBranchAddress("ClusterX", &fECALClusterX);
      NDGArRecoTree->SetBranchAddress("ClusterY", &fECALClusterY);
      NDGArRecoTree->SetBranchAddress("ClusterZ", &fECALClusterZ);

      NDGArRecoTree->SetBranchAddress("ClusterIDNumber", &fECALClusterIDNumber);

      NDGArRecoTree->SetBranchAddress("ClusterEnergy", &fECALClusterEnergy);
      NDGArRecoTree->SetBranchAddress("ClusterNhits", &fECALClusterNhits);

      NDGArRecoTree->SetBranchAddress("ClusterMCindex", &fECALClusterMCindex);
      NDGArRecoTree->SetBranchAddress("ClusterMCfrac", &fECALClusterMCfrac);

      //ECAL-track associations
      NDGArRecoTree->SetBranchAddress("ECALAssn_ClusIDNumber", &fECALAssn_ClusterID);
      NDGArRecoTree->SetBranchAddress("ECALAssn_TrackIDNumber", &fECALAssn_TrackID);

    } else {
      NDGArRecoTree = NULL;
      std::cerr << "Did not find input ND-GAr reco file you provided: " << NDGArRecoFilename << std::endl;
      std::cerr << "Are you sure it exists?" << std::endl;
      throw;
    }
  }

  void NDGArRecoBranchFiller::_FillRecoBranches(std::size_t ii, 
					       caf::StandardRecord &sr,
					       const cafmaker::Params &par) const
  {

   NDGArRecoTree->GetEntry(ii);
   std::cout << "Event number: " << fEvent << std::endl;
   //std::cout << "Nu PDG: " << fMCNuPDG->at(0) << std::endl;
   //std::cout << "Number of reco tracks: " << fTrackStartX->size() << std::endl;
   //std::cout << "Number of reco ECAL clusters: " << fECALClusterX->size() << std::endl;

   //std::cout << "Number of MC Particles: " << fMCTrkID->size() << std::endl;

   int n_tracks = fTrackStartX->size();
   sr.nd.gar.ntracks = n_tracks;

   int pid_counter = 0;
   caf::SRGArTrack track;
   for (size_t iTrack=0; iTrack<n_tracks; iTrack++){
      
      caf::SRVector3D start(fTrackStartX->at(iTrack), fTrackStartY->at(iTrack), fTrackStartZ->at(iTrack));
      track.start = start;
      
      caf::SRVector3D end(fTrackEndX->at(iTrack), fTrackEndY->at(iTrack), fTrackEndZ->at(iTrack));
      track.end = end;
      
      caf::SRVector3D dir(fTrackStartPX->at(iTrack), fTrackStartPY->at(iTrack), fTrackStartPZ->at(iTrack));
      track.dir = dir.Unit();

      caf::SRVector3D enddir(fTrackEndPX->at(iTrack), fTrackEndPY->at(iTrack), fTrackEndPZ->at(iTrack));
      track.enddir = enddir.Unit();

      track.len_cm_fwd = fTrackLenF->at(iTrack);
      track.len_cm_bkwd = fTrackLenB->at(iTrack);
      track.p_fwd = fTrackPF->at(iTrack);
      track.p_bkwd = fTrackPB->at(iTrack);
      track.dEdx_fwd = fTrackAvgIonF->at(iTrack);
      track.dEdx_bkwd = fTrackAvgIonB->at(iTrack);

      track.garsoft_trk_id = fTrackIDNumber->at(iTrack);
      track.clusters_in_track = fTrackNClusters->at(iTrack);

      for (size_t iPID=6*pid_counter; iPID<6*(pid_counter+1); ++iPID){
        track.pid_fwd.push_back(fTrackPIDF->at(iPID));
        track.pid_prob_fwd.push_back(fTrackPIDProbF->at(iPID));
        track.pid_bkwd.push_back(fTrackPIDB->at(iPID));
        track.pid_prob_bkwd.push_back(fTrackPIDProbB->at(iPID));
      }
      ++pid_counter;

      int iMCParticleTrack = fTrackMCindex->at(iTrack);
      caf::SRParticleTruth mc_true_track;
      if(iMCParticleTrack >= 0){
        mc_true_track.trkid = fMCTrkID->at(iMCParticleTrack);
        mc_true_track.pdg = fMCPDG->at(iMCParticleTrack);

        caf::SRVector3D mc_start(fMCStartX->at(iMCParticleTrack), fMCStartY->at(iMCParticleTrack), fMCStartZ->at(iMCParticleTrack));
        mc_true_track.start = mc_start;

        caf::SRVector3D mc_p(fMCStartPX->at(iMCParticleTrack), fMCStartPY->at(iMCParticleTrack), fMCStartPZ->at(iMCParticleTrack));
        mc_true_track.p.E  = mc_p.Mag();
        mc_true_track.p.px = mc_p.X();
        mc_true_track.p.py = mc_p.Y();
        mc_true_track.p.pz = mc_p.Z();

        int iMCMotherTrack = fMCMotherIndex->at(iMCParticleTrack);
        if(iMCMotherTrack >= 0){
          mc_true_track.motherpdg = fMCPDGMother->at(iMCParticleTrack);

          caf::SRVector3D mc_motherp(fMCStartPX->at(iMCMotherTrack), fMCStartPY->at(iMCMotherTrack), fMCStartPZ->at(iMCMotherTrack));
          mc_true_track.motherp.E  = mc_motherp.Mag();
          mc_true_track.motherp.px = mc_motherp.X();
          mc_true_track.motherp.py = mc_motherp.Y();
          mc_true_track.motherp.pz = mc_motherp.Z();
        } else{
          mc_true_track.motherpdg = fMCNuPDG->at(0);

          caf::SRVector3D mc_nup(fMCNuPX->at(0), fMCNuPY->at(0), fMCNuPZ->at(0));
          mc_true_track.motherp.E  = mc_nup.Mag();
          mc_true_track.motherp.px = mc_nup.X();
          mc_true_track.motherp.py = mc_nup.Y();
          mc_true_track.motherp.pz = mc_nup.Z();
        }

        track.truth_fraction = fTrackMCfrac->at(iTrack);
      } else {
        mc_true_track.trkid = -1;
        mc_true_track.pdg = -1;
        mc_true_track.motherpdg = -1;

        caf::SRVector3D mc_start(-1, -1, -1);
        mc_true_track.start = mc_start;

        caf::SRVector3D mc_p(-1, -1, -1);
        mc_true_track.p.E  = mc_p.Mag();
        mc_true_track.p.px = mc_p.X();
        mc_true_track.p.py = mc_p.Y();
        mc_true_track.p.pz = mc_p.Z();

        caf::SRVector3D mc_motherp(-1, -1, -1);
        mc_true_track.motherp.E  = mc_motherp.Mag();
        mc_true_track.motherp.px = mc_motherp.X();
        mc_true_track.motherp.py = mc_motherp.Y();
        mc_true_track.motherp.pz = mc_motherp.Z();

        track.truth_fraction = -1.;
      }

      track.truth = mc_true_track;

      sr.nd.gar.tracks.push_back(track);

   }

   int n_clusters = fECALClusterX->size();
   sr.nd.gar.nclusters = n_clusters;

   int n_assns = fECALAssn_ClusterID->size();
   caf::SRGArECAL cluster;
   for (size_t iECAL=0; iECAL<n_clusters; iECAL++){
      
      caf::SRVector3D position(fECALClusterX->at(iECAL), fECALClusterY->at(iECAL), fECALClusterZ->at(iECAL));
      cluster.position = position;

      cluster.E = fECALClusterEnergy->at(iECAL);
      cluster.hits_in_cluster = fECALClusterNhits->at(iECAL);

      cluster.garsoft_ecal_id = fECALClusterIDNumber->at(iECAL);

      for (size_t iAssn=0; iAssn<n_assns; ++iAssn){
          if (cluster.garsoft_ecal_id == fECALAssn_ClusterID->at(iAssn)){
              cluster.garsoft_trk_assn = fECALAssn_TrackID->at(iAssn);
          }
      }

      int iMCParticleCluster = fECALClusterMCindex->at(iECAL);
      caf::SRParticleTruth mc_true_cluster;
      if(iMCParticleCluster >= 0){
        mc_true_cluster.trkid = fMCTrkID->at(iMCParticleCluster);
        mc_true_cluster.pdg = fMCPDG->at(iMCParticleCluster);

        caf::SRVector3D mc_start(fMCStartX->at(iMCParticleCluster), fMCStartY->at(iMCParticleCluster), fMCStartZ->at(iMCParticleCluster));
        mc_true_cluster.start = mc_start;

        caf::SRVector3D mc_p(fMCStartPX->at(iMCParticleCluster), fMCStartPY->at(iMCParticleCluster), fMCStartPZ->at(iMCParticleCluster));
        mc_true_cluster.p.E  = mc_p.Mag();
        mc_true_cluster.p.px = mc_p.X();
        mc_true_cluster.p.py = mc_p.Y();
        mc_true_cluster.p.pz = mc_p.Z();

        int iMCMotherCluster = fMCMotherIndex->at(iMCParticleCluster);
        if(iMCMotherCluster >= 0){
          mc_true_cluster.motherpdg = fMCPDGMother->at(iMCParticleCluster);
          
          caf::SRVector3D mc_motherp(fMCStartPX->at(iMCMotherCluster), fMCStartPY->at(iMCMotherCluster), fMCStartPZ->at(iMCMotherCluster));
          mc_true_cluster.motherp.E  = mc_motherp.Mag();
          mc_true_cluster.motherp.px = mc_motherp.X();
          mc_true_cluster.motherp.py = mc_motherp.Y();
          mc_true_cluster.motherp.pz = mc_motherp.Z();
        } else{
          mc_true_cluster.motherpdg = fMCNuPDG->at(0);

          caf::SRVector3D mc_nup(fMCNuPX->at(0), fMCNuPY->at(0), fMCNuPZ->at(0));
          mc_true_cluster.motherp.E  = mc_nup.Mag();
          mc_true_cluster.motherp.px = mc_nup.X();
          mc_true_cluster.motherp.py = mc_nup.Y();
          mc_true_cluster.motherp.pz = mc_nup.Z();
        }

        cluster.truth_fraction = fECALClusterMCfrac->at(iECAL);
      } else {
        mc_true_cluster.trkid = -1;
        mc_true_cluster.pdg = -1;
        mc_true_cluster.motherpdg = -1;

        caf::SRVector3D mc_start(-1, -1, -1);
        mc_true_cluster.start = mc_start;

        caf::SRVector3D mc_p(-1, -1, -1);
        mc_true_cluster.p.E  = mc_p.Mag();
        mc_true_cluster.p.px = mc_p.X();
        mc_true_cluster.p.py = mc_p.Y();
        mc_true_cluster.p.pz = mc_p.Z();

        caf::SRVector3D mc_motherp(-1, -1, -1);
        mc_true_cluster.motherp.E  = mc_motherp.Mag();
        mc_true_cluster.motherp.px = mc_motherp.X();
        mc_true_cluster.motherp.py = mc_motherp.Y();
        mc_true_cluster.motherp.pz = mc_motherp.Z();

        cluster.truth_fraction = -1.;
      }

      cluster.truth = mc_true_cluster;

      sr.nd.gar.clusters.push_back(cluster);

   }

   //legacy variables
   //std::cout << "Event number: " << fEvent << std::endl;
   //std::cout << "Number of particles: " << fTrackStartX->size() << std::endl;

   //using forward variables only!

   // int pi_pl_mult = 0;
   // int pi_min_mult = 0;
   //
   // for (size_t i=0; i< fTrackStartX->size(); ++i){
   //    sr.nd.gar.pdg.push_back(fTrackPIDCheatedF->at(i));
   //    sr.nd.gar.trkLen.push_back(fTrackLenF->at(i));
   //    sr.nd.gar.ptrue.push_back(fTrackPF->at(i));       //filling with "reco" momentum for the moment
   //    if (fTrackPIDCheatedF->at(i) == 211){
   //      ++pi_pl_mult;
   //    } else if (fTrackPIDCheatedF->at(i) == -211){
   //      ++pi_min_mult;
   //    }
   // }
   // sr.nd.gar.gastpc_pi_pl_mult = pi_pl_mult;
   // sr.nd.gar.gastpc_pi_min_mult = pi_min_mult;
  }
}
