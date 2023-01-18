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
      NDGArRecoTree->SetBranchAddress("TrackStartX", &fTrackStartX);
      NDGArRecoTree->SetBranchAddress("TrackPIDCheatedF", &fTrackPIDCheatedF);
      NDGArRecoTree->SetBranchAddress("TrackLenF", &fTrackLenF);
      NDGArRecoTree->SetBranchAddress("TrackPF", &fTrackPF);
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
   std::cout << "Number of particles: " << fTrackStartX->size() << std::endl;
    
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
