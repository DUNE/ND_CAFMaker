/// Fill MINERvA reco branches from MINERvA reco output.
///
/// \author  N.Roy <noeroy@yorku.ca>
/// \date    Nov. 2022

#ifndef ND_CAFMAKER_MINERvARecoBranchFiller_H
#define ND_CAFMAKER_MINERvARecoBranchFiller_H

#include <iostream>
#include <vector>
#include <algorithm>

// The virtual base class
#include "reco/IRecoBranchFiller.h"
#include "truth/FillTruth.h"

// File handlers from ROOT
#include "TFile.h"
#include "TTree.h"

// The duneanaobj includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
// #include "duneanaobj/StandardRecord/SRVector3D.h"

namespace caf
{
  class SRTrueInteraction;
  class SRTrueParticle;
}

namespace cafmaker
{
  class MINERvARecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      MINERvARecoBranchFiller(const std::string & minervaRecoFilename);

      std::deque<Trigger> GetTriggers(int triggerType) const  override;


      ~MINERvARecoBranchFiller() {
        delete MnvRecoTree;
        fMnvRecoFile->Close();
        delete fMnvRecoFile;
        MnvRecoTree = NULL;
        fMnvRecoFile = NULL;
      }

    private:
      void _FillRecoBranches(const Trigger &trigger,
                              caf::StandardRecord &sr,
                              const cafmaker::Params &par,
                              const TruthMatcher *truthMatch= nullptr) const override;

      std::map<int, std::vector<caf::SRTrack>> fill_track(caf::StandardRecord &sr, int & max_slice, const TruthMatcher *truthMatch) const;
      std::map<int, std::vector<caf::SRShower>> fill_shower(caf::StandardRecord &sr, int & max_slice, const TruthMatcher *truthMatch) const;

      void find_truth_shower(caf::StandardRecord &sr, caf::SRShower & sh, int shower_id, const TruthMatcher *truthMatch) const;
      void find_truth_track(caf::StandardRecord &sr, caf::SRTrack & t, int track_id, const TruthMatcher *truthMatch) const;
      void fillTrueInteraction(caf::StandardRecord &sr) const;
      TFile *fMnvRecoFile;
      TTree *MnvRecoTree;
      
      Double_t        offsetX; //Minerva ref point and Edepsim not equal
      Double_t        offsetY; //Minerva ref point and Edepsim not equal
      Double_t        offsetZ; //Minerva ref point and Edepsim not equal

      Int_t           ev_trigger_type;
      Int_t           ev_gl_gate;
      Int_t           ev_gps_time_sec;
      Int_t           ev_gps_time_usec;

      //Tracks variables
      Int_t           n_tracks;
      Int_t           trk_index[100];   //[n_tracks]
      Int_t           trk_type[100];   //[n_tracks]
      Int_t           trk_patrec[100];   //[n_tracks]
      Int_t           trk_time_slice[100];   //[n_tracks]
      Double_t        trk_vis_energy[100];   //[n_tracks]
      Double_t        trk_theta[100];   //[n_tracks]
      Double_t        trk_phi[100];   //[n_tracks]
      Int_t           trk_hits[100];   //[n_tracks]
      Int_t           trk_dof[100];   //[n_tracks]
      Double_t        trk_chi2perDof[100];   //[n_tracks]
      Double_t        trk_fitMass[100];   //[n_tracks]
      Int_t           trk_nodes[100];   //[n_tracks]
      Double_t        trk_node_X[100][300];   //[n_tracks]
      Double_t        trk_node_Y[100][300];   //[n_tracks]
      Double_t        trk_node_Z[100][300];   //[n_tracks]
      Double_t        trk_node_aX[100][300];   //[n_tracks]
      Double_t        trk_node_aY[100][300];   //[n_tracks]
      Double_t        trk_node_qOverP[100][300];   //[n_tracks]
      Double_t        trk_node_chi2[100][300];   //[n_tracks]
      Int_t           trk_node_cluster_idx[100][300];   //[n_tracks]

      //Showers variables
      Int_t           n_blobs_id;
      Int_t           blob_id_idx[1000];   //[n_blobs_id]
      Int_t           blob_id_subdet[1000];   //[n_blobs_id]
      Int_t           blob_id_history[1000];   //[n_blobs_id]
      Int_t           blob_id_size[1000];   //[n_blobs_id]
      Int_t           blob_id_patrec[1000];   //[n_blobs_id]
      Double_t        blob_id_e[1000];   //[n_blobs_id]
      Double_t        blob_id_time[1000];   //[n_blobs_id]
      Int_t           blob_id_time_slice[1000];   //[n_blobs_id]
      Double_t        blob_id_startpoint_x[1000];   //[n_blobs_id]
      Double_t        blob_id_startpoint_y[1000];   //[n_blobs_id]
      Double_t        blob_id_startpoint_z[1000];   //[n_blobs_id]
      Double_t        blob_id_centroid_x[1000];   //[n_blobs_id]
      Double_t        blob_id_centroid_y[1000];   //[n_blobs_id]
      Double_t        blob_id_centroid_z[1000];   //[n_blobs_id]
      Int_t           blob_id_clus_idx[1000][1500];   //[n_blobs_id]

      //Trajectories (Truth variables)
      Int_t           n_mc_trajectories;
      Int_t           mc_traj_trkid[180];   //[n_mc_trajectories]
      Int_t           mc_traj_parentid[180];   //[n_mc_trajectories]
      Int_t           mc_traj_pdg[180];   //[n_mc_trajectories]
      Double_t        mc_traj_hit_e[180];   //[n_mc_trajectories]
      Int_t           mc_traj_npoints[180];   //[n_mc_trajectories]
      Int_t           mc_traj_edepsim_trkid[180];   //[n_mc_trajectories]
      Long64_t        mc_traj_edepsim_eventid[180];   //[n_mc_trajectories]
      Double_t        mc_traj_point_x[180][5];   //[n_mc_trajectories]
      Double_t        mc_traj_point_y[180][5];   //[n_mc_trajectories]
      Double_t        mc_traj_point_z[180][5];   //[n_mc_trajectories]
      Double_t        mc_traj_point_t[180][5];   //[n_mc_trajectories]
      Double_t        mc_traj_point_px[180][5];   //[n_mc_trajectories]
      Double_t        mc_traj_point_py[180][5];   //[n_mc_trajectories]
      Double_t        mc_traj_point_pz[180][5];   //[n_mc_trajectories]
      Double_t        mc_traj_point_E[180][5];   //[n_mc_trajectories]
   

      //Cluster list & hit list for Trajectory matching
      Int_t           clus_id_size[7500];   //[n_clusters_id]
      Int_t           clus_id_hits_idx[7500][60];   //[n_clusters_id]
      Int_t           mc_id_mchit_trkid[50000][2];   //[n_mc_id_digits]
      Double_t        mc_id_mchit_dE[50000][2];   //[n_mc_id_digits]

      mutable std::vector<cafmaker::Trigger> fTriggers;
      mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;    ///< the last trigger requested using _FillRecoBranches()
  };

}
#endif //ND_CAFMAKER_MINERvARecoBranchFiller_H
