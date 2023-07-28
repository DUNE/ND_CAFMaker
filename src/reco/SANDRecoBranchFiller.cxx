/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi
/// \date    Apr. 2022
///

#include <iostream>

#include "SANDRecoBranchFiller.h"
#include "truth/FillTruth.h"

#ifdef ENABLE_SAND
#warning Including SANDRecoBranchFiller in build

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include "SANDReco/SANDRecord.h"

namespace cafmaker
{
  
  SANDRecoBranchFiller::SANDRecoBranchFiller(const std::string &SANDRecoFilename)
  {  
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

#else // ENABLE_SAND

#warning Not configured to build SANDRecoBranchFiller. Must set SANDRECO_INC and SANDRECO_LIB environment variables

namespace {
  void error_msg()
  {
    std::cerr << "\n\nSAND Reco support was not enabled in your build. \n"
              << " Either avoid setting `nd_cafmaker.CAFMakerSettings.SANDRecoFile` in your FCL\n"
              << " or set $SANDRECO_INC and $SANDRECO_LIB in your environment and do a clean rebuild of ND_CAFMaker...\n";

  }
}

cafmaker::SANDRecoBranchFiller::
SANDRecoBranchFiller(const std::string&)
{
  error_msg();
  abort();
}

void cafmaker::SANDRecoBranchFiller::
_FillRecoBranches(std::size_t, caf::StandardRecord &, const cafmaker::Params &, const TruthMatcher *truthMatcher) const
{
  error_msg();
  abort();
}

#endif // ENABLE_SAND
