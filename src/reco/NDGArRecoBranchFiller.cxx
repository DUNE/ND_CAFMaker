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

      NDGArRecoTree->SetBranchAddress("TrackStartPx", &fTrackStartPx);
      NDGArRecoTree->SetBranchAddress("TrackStartPy", &fTrackStartPy);
      NDGArRecoTree->SetBranchAddress("TrackStartPz", &fTrackStartPz);
      NDGArRecoTree->SetBranchAddress("TrackEndPx", &fTrackEndPx);
      NDGArRecoTree->SetBranchAddress("TrackEndPy", &fTrackEndPy);
      NDGArRecoTree->SetBranchAddress("TrackEndPz", &fTrackEndPz);

      NDGArRecoTree->SetBranchAddress("TrackIDNumber", &fTrackIDNumber);

      NDGArRecoTree->SetBranchAddress("TrackLenF", &fTrackLenF);
      NDGArRecoTree->SetBranchAddress("TrackLenB", &fTrackLenB);
      NDGArRecoTree->SetBranchAddress("TrackPF", &fTrackPF);
      NDGArRecoTree->SetBranchAddress("TrackPB", &fTrackPB);
      NDGArRecoTree->SetBranchAddress("TrackAvgIonF", &fTrackAvgIonF);
      NDGArRecoTree->SetBranchAddress("TrackAvgIonB", &fTrackAvgIonB);

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
   std::cout << "Event number: " << fEvent << std::endl;
   std::cout << "Number of reco tracks: " << fTrackStartX->size() << std::endl;
   std::cout << "Number of reco ECAL clusters: " << fECALClusterX->size() << std::endl;

   int n_tracks = fTrackStartX->size();
   sr.nd.gar.ntracks = n_tracks;
   for (size_t i=0; i<n_tracks; i++){
      caf::SRGArTrack track;

      caf::SRVector3D start(fTrackStartX->at(i), fTrackStartY->at(i), fTrackStartZ->at(i));
      track.start = start;
      caf::SRVector3D end(fTrackEndX->at(i), fTrackEndY->at(i), fTrackEndZ->at(i));
      track.end = end;
      caf::SRVector3D dir(fTrackStartPx->at(i), fTrackStartPy->at(i), fTrackStartPz->at(i));
      track.start = dir;
      caf::SRVector3D enddir(fTrackEndPx->at(i), fTrackEndPy->at(i), fTrackEndPz->at(i));
      track.start = enddir;

      track.len_cm_F = fTrackLenF->at(i);
      track.len_cm_B = fTrackLenB->at(i);
      track.p_F = fTrackPF->at(i);
      track.p_B = fTrackPB->at(i);
      track.dEdx_F = fTrackAvgIonF->at(i);
      track.dEdx_B = fTrackAvgIonB->at(i);

      std::vector<int> pid_F = fTrackPIDF->at(i);
      track.pid_F = pid_F;
      std::vector<float_t> pid_prob_F = fTrackPIDProbF->at(i);
      track.pid_prob_F = pid_prob_F;
      std::vector<int> pid_B = fTrackPIDB->at(i);
      track.pid_F = pid_B;
      std::vector<float_t> pid_prob_B = fTrackPIDProbB->at(i);
      track.pid_prob_B = pid_prob_B;

      sr.nd.gar.tracks.push_back(track);
   }

   int n_clusters = fECALClusterX->size();
   sr.nd.gar.nclusters = n_clusters;

   int n_assns = fECALAssn_ClusterID->size();
   
   for (size_t i=0; i<n_clusters; i++){
      caf::SRGArECAL cluster;

      caf::SRVector3D position(fECALClusterX->at(i), fECALClusterY->at(i), fECALClusterZ->at(i));
      cluster.position = position;

      cluster.E = fECALClusterEnergy->at(i);
      cluster.hits_in_cluster = fECALClusterNhits->at(i);

      cluster.ecal_id = fECALClusterIDNumber->at(i);

      for (size_t j=0; j<n_assns; ++j){
          if (cluster.ecal_id == fECALAssn_ClusterID->at(j)){
              cluster.trk_assn = fECALAssn_TrackID->at(j);
          }
      }

   }

   //legacy variables
   //std::cout << "Event number: " << fEvent << std::endl;
   //std::cout << "Number of particles: " << fTrackStartX->size() << std::endl;
    
   //todo: currently filling pseudo reco variables in SRGAr.h
   //once GAr reco is integrated they should be replaced with proper reco

   //using forward variables only!

   int pi_pl_mult = 0;
   int pi_min_mult = 0;
  
   for (size_t i=0; i< fTrackStartX->size(); ++i){
      sr.nd.gar.pdg.push_back(fTrackPIDCheatedF->at(i));
      sr.nd.gar.trkLen.push_back(fTrackLenF->at(i));
      sr.nd.gar.ptrue.push_back(fTrackPF->at(i));       //filling with "reco" momentum for the moment
      if (fTrackPIDCheatedF->at(i) == 211){
        ++pi_pl_mult;
      } else if (fTrackPIDCheatedF->at(i) == -211){
        ++pi_min_mult;
      }
   }
   sr.nd.gar.gastpc_pi_pl_mult = pi_pl_mult;
   sr.nd.gar.gastpc_pi_min_mult = pi_min_mult;
  }
}
