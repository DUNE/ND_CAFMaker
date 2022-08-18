/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi
/// \date    Apr. 2022
///

#include "SANDRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TSystem.h"
#include "SANDReco/SANDRecord.h" //fixme

namespace cafmaker
{
  
  SANDRecoBranchFiller::SANDRecoBranchFiller(const std::string &SANDRecoFilename,
					     const std::string &SANDLibFilename)
  {  
    
    gSystem->Load(SANDLibFilename.c_str());
    
    fSANDRecoFile = new TFile(SANDRecoFilename.c_str());
    fTree = (TTree*) fSANDRecoFile->Get("tEvent");

    fEvent = new event;
    fTree->SetBranchAddress("event", &fEvent);

    if(!fSANDRecoFile->IsZombie())
      SetConfigured(true);
  }

  void SANDRecoBranchFiller::_FillRecoBranches(std::size_t ii, 
					       caf::StandardRecord &sr,
					       const cafmaker::Params &par) const
  {

   fTree->GetEntry(ii);
    
   //todo: currently filling simple variables
   //rewrite once sr.nd.sand exists in StandardRecord

   // neutrino energy
   sr.Ev_reco = fEvent->Enureco*0.001; //GeV
  
   std::vector<particle> particle_event = fEvent->particles; 
   bool foundLepton = false;
   for ( auto it = particle_event.begin(); it != particle_event.end(); ++it){

      // primary lepton
      if( abs(((*it).pdg) == 13 || abs((*it).pdg) == 11) && (*it).primary == 1){
        sr.reco_lepton_pdg = (*it).pdg;
	sr.Elep_reco = (*it).Ereco*0.001;  //GeV
        foundLepton = true;
      }

      // other species
      // ...
	
   }
    
   // event flags
   if(!foundLepton){ //flags as nc
     sr.reco_numu = 0; sr.reco_nue=0; sr.reco_nc=1;
     sr.reco_lepton_pdg = -1;	
     sr.Elep_reco = -1;	
   }else if(sr.reco_lepton_pdg == 13){ //numu
     sr.reco_numu = 1; sr.reco_nue = 0; sr.reco_nc = 0; 
   }else if(sr.reco_lepton_pdg == 11){ //nue
     sr.reco_numu = 0; sr.reco_nue = 1; sr.reco_nc = 0;
   }


  }
}
