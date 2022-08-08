/// Get and set SAND event
///
/// \author  L. Di Noto 
/// \date    Apr. 2022

#include "SANDEvt.h"

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// Standard Record format
#include "duneanaobj/StandardRecord/StandardRecord.h"

// ND_CAFMaker
#include "dumpTree.h"
#include "CAF.h"
#include "Params.h"


SANDEvt* SANDEvt::me = 0;

SANDEvt* SANDEvt::Get() {
	if (!me) me=new SANDEvt();
	return me;

}

SANDEvt::SANDEvt(){}

SANDEvt::~SANDEvt() {}

bool SANDEvt::GetSANDEvtLepton(particle &Lepton){
 
  bool foundLepton=false;
  std::vector<particle> particle_event=fEvent->particles;
  for (std::vector<particle>::iterator it = particle_event.begin() ; it != particle_event.end(); ++it){
  	if( abs(((*it).pdg) == 13 || abs((*it).pdg)== 11)  && (*it).primary==1) { Lepton=(*it); foundLepton=true;}
	}  
  return foundLepton;
}

double SANDEvt::GetSANDEvtEnureco(){
   return fEvent->Enureco;

}

double SANDEvt::GetSANDEvtErecopdg(int pdg){

  std::vector<particle> particle_event=fEvent->particles;
  double Erecopdg = 0;	
  for (std::vector<particle>::iterator it = particle_event.begin() ; it != particle_event.end(); ++it){
  	
  	double ke = (*it).Etrue-(*it).mass; //MeV

	if( (*it).pdg == pdg && (*it).primary == 1) Erecopdg += ke;
  
	else if( pdg == 0 && (*it).primary == 1 //other
		&& std::abs((*it).pdg) != 13
		&& std::abs((*it).pdg) != 11 
		&& std::abs((*it).pdg) != 211 
		&& (*it).pdg != 111
		&& (*it).pdg != 2112
		&& (*it).pdg != 2212 ) Erecopdg += ke;
}

  return Erecopdg; 
}
