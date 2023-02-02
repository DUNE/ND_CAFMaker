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

      //ECAL-related info
      NDGArRecoTree->SetBranchAddress("ClusterX", &fECALClusterX);
      NDGArRecoTree->SetBranchAddress("ClusterY", &fECALClusterY);
      NDGArRecoTree->SetBranchAddress("ClusterZ", &fECALClusterZ);

      NDGArRecoTree->SetBranchAddress("ClusterIDNumber", &fECALClusterIDNumber);

      NDGArRecoTree->SetBranchAddress("ClusterEnergy", &fECALClusterEnergy);
      NDGArRecoTree->SetBranchAddress("ClusterNhits", &fECALClusterNhits);

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
   //std::cout << "Event number: " << fEvent << std::endl;
   //std::cout << "Number of reco tracks: " << fTrackStartX->size() << std::endl;
   //std::cout << "Number of reco ECAL clusters: " << fECALClusterX->size() << std::endl;

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

      track.len_cm_F = fTrackLenF->at(iTrack);
      track.len_cm_B = fTrackLenB->at(iTrack);
      track.p_F = fTrackPF->at(iTrack);
      track.p_B = fTrackPB->at(iTrack);
      track.dEdx_F = fTrackAvgIonF->at(iTrack);
      track.dEdx_B = fTrackAvgIonB->at(iTrack);

      track.trk_id = fTrackIDNumber->at(iTrack);
      track.clusters_in_track = fTrackNClusters->at(iTrack);

      for (size_t iPID=6*pid_counter; iPID<6*(pid_counter+1); ++iPID){
        track.pid_F.push_back(fTrackPIDF->at(iPID));
        track.pid_prob_F.push_back(fTrackPIDProbF->at(iPID));
        track.pid_B.push_back(fTrackPIDB->at(iPID));
        track.pid_prob_B.push_back(fTrackPIDProbB->at(iPID));
      }
      ++pid_counter;

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

      cluster.ecal_id = fECALClusterIDNumber->at(iECAL);

      for (size_t iAssn=0; iAssn<n_assns; ++iAssn){
          if (cluster.ecal_id == fECALAssn_ClusterID->at(iAssn)){
              cluster.trk_assn = fECALAssn_TrackID->at(iAssn);
          }
      }

      sr.nd.gar.clusters.push_back(cluster);

   }

   //legacy variables
   //std::cout << "Event number: " << fEvent << std::endl;
   //std::cout << "Number of particles: " << fTrackStartX->size() << std::endl;
    
   //todo: currently filling pseudo reco variables in SRGAr.h
   //once GAr reco is integrated they should be replaced with proper reco

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
