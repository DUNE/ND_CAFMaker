/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi
/// \date    Jan. 2023
///

#include <iostream>

#include "SANDRecoBranchFiller.h"
#include "truth/FillTruth.h"

#ifdef ENABLE_SAND
#warning Including SANDRecoBranchFiller in build

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRRecoParticle.h"
#include "duneanaobj/StandardRecord/SRParticleTruth.h"


#include "TFile.h"

#include "TFile.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include "SANDRecord.h"

namespace cafmaker
{
  
  SANDRecoBranchFiller::SANDRecoBranchFiller(const std::string &SANDRecoFilename)
    : IRecoBranchFiller("SAND")
  {  
    fSANDRecoFile = new TFile(SANDRecoFilename.c_str());
    fTree = (TTree*) fSANDRecoFile->Get("tEvent");

    fEvent = new event;
    fTree->SetBranchAddress("event", &fEvent);

    if(!fSANDRecoFile->IsZombie())
      SetConfigured(true);
  }

  void SANDRecoBranchFiller::_FillRecoBranches(std::size_t N, std::size_t ii, 
					       caf::StandardRecord &sr,
					       const cafmaker::Params &par) const
  {

   fTree->GetEntry(ii);
    
   //todo: currently filling simple variables
   //rewrite once sr.nd.sand exists in StandardRecord

   // neutrino energy
   sr.Ev_reco = fEvent->Enureco*0.001; //GeV

   //sr.nd.sand.ntracks=10; 
   std::vector<particle> particle_event = fEvent->particles; 
   bool foundLepton = false;

  sr.nd.sand.nparticles=particle_event.size(); 
  caf::SRRecoParticle recop;

  for ( auto it = particle_event.begin(); it != particle_event.end(); ++it){
	recop.pdg=(*it).pdg;
	
	recop.primary=(*it).primary;
	recop.reco_trkid=(*it).tid;
	recop.mother_trkid=(*it).parent_tid;
	recop.p.E=(*it).Ereco*0.001; //in GeV
	recop.p.px=(*it).pxreco*0.001; //in GeV
	recop.p.py=(*it).pyreco*0.001; //in GeV
	recop.p.pz=(*it).pzreco*0.001; //in GeV
	caf::SRVector3D start_point((*it).xreco, (*it).yreco,(*it).zreco);  //unità?
   	recop.start=start_point;
	
	recop.trkid_best_match=(*it).tid;  //non abbiamo ancora un best match. TO DO

	caf::SRParticleTruth truep;
 	truep.trkid=(*it).tid;
	truep.pdg=(*it).pdg;
	truep.p.E=(*it).Etrue*0.001; //in GeV
	truep.p.px=(*it).pxtrue*0.001; //in GeV
        truep.p.py=(*it).pytrue*0.001; //in GeV
        truep.p.pz=(*it).pztrue*0.001; //in GeV
	caf::SRVector3D start_point_true((*it).xtrue, (*it).ytrue,(*it).ztrue);  //unità?
        truep.start=start_point_true;
	recop.particle_true=truep;

	sr.nd.sand.recoparticles.push_back(recop);
	
	

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

namespace cafmaker
{
  SANDRecoBranchFiller::SANDRecoBranchFiller(const std::string &)
    : IRecoBranchFiller("SAND")
  {
    error_msg();
    abort();
  }

  void SANDRecoBranchFiller::
  _FillRecoBranches(const Trigger &, caf::StandardRecord &, const cafmaker::Params &,
                    const TruthMatcher *truthMatcher) const
  {
    error_msg();
    abort();
  }

  // todo: this is a placeholder
  std::deque<Trigger> SANDRecoBranchFiller::GetTriggers(int triggerType) const
  {
    return std::deque<Trigger>();
  }
}

#endif // ENABLE_SAND
