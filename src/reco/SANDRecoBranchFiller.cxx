/// Fill reco branches using SAND reco data
///

#include "SANDRecoBranchFiller.h"

#include "TF1.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include <TSystem.h>

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "Params.h"
#include "struct.h"    
#include "SANDEvt.h"

namespace cafmaker
{
  const double mmu = 0.1056583745;
  void SANDRecoBranchFiller::_FillRecoBranches(std::size_t ii, //elemento del tree da guardare.. 
                                                        caf::StandardRecord &sr,
                                                        const cafmaker::Params &par) const
  {
   //std::cout<<"riempio con i dati di SAND!"<<std::endl;
   //devo aprire sand-reco da cui prendere i dati corrispondenti all'elemento size_t del tree
    std::string sandFile;
    par().cafmaker().sandRecoFile(sandFile);

    particle Recolepton;
    if(SANDEvt::Get()->GetSANDEvtLepton(Recolepton)){
	   sr.reco_lepton_pdg = Recolepton.pdg;
 	   sr.Elep_reco =Recolepton.Ereco*0.001;  //in MeV
    if(Recolepton.pdg==13){ sr.reco_numu = 1; sr.reco_nue=0; sr.reco_nc=0;}
    if(Recolepton.pdg==11){  sr.reco_numu = 0; sr.reco_nue=1; sr.reco_nc=0;}

  }else{
//Lepton not found in reco event --> it is assumed a nc event
	 sr.reco_numu = 0; sr.reco_nue=0; sr.reco_nc=1;
	 sr.reco_lepton_pdg=-1;	
    	 sr.Elep_reco=-1;	

	}
 
    //neutrino energy
    sr.Ev_reco = SANDEvt::Get()->GetSANDEvtEnureco()*0.001 ; 

    //reco energy by species
    sr.eRecoP = SANDEvt::Get()->GetSANDEvtErecopdg(2212)*0.001;   
    sr.eRecoN = SANDEvt::Get()->GetSANDEvtErecopdg(2112)*0.001;   
    sr.eRecoPip = SANDEvt::Get()->GetSANDEvtErecopdg(211)*0.001;   
    sr.eRecoPim = SANDEvt::Get()->GetSANDEvtErecopdg(-211)*0.001;   
    sr.eRecoPi0 = SANDEvt::Get()->GetSANDEvtErecopdg(111)*0.001;   
    sr.eRecoOther = SANDEvt::Get()->GetSANDEvtErecopdg(0)*0.001;

/*
   sr.reco_q = 100;
   sr.muon_contained = 100; 
   sr.muon_tracker = 100; 
   sr.muon_ecal = 100; 
   sr.muon_exit = 1;
*/
 }
}
