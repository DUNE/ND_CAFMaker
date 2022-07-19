/// \file FillTruthForSAND.cxx
///
/// Fill truth branches for SAND data
///
/// \author  L. Di Noto based on code by J. Wolcott <jwolcott@fnal.gov> and C. Marshall <chris.marshall@rochester.edu>
/// \date    Apr. 2022

#include "FillTruthForSAND.h"

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// GENIE
#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCEventRecord.h"

// Standard Record format
#include "duneanaobj/StandardRecord/StandardRecord.h"

// ND_CAFMaker
#include "dumpTree.h"
#include "CAF.h"
#include "Params.h"

// Fill truth info
void fillTruthForSAND(caf::StandardRecord& sr,
               TTree *intree,
               TTree * gtree,
               event * evt,
	       const genie::NtpMCEventRecord * mcrec,
               const cafmaker::Params &par,
               nusyst::response_helper& rh)
{
  // get GENIE event record
  genie::EventRecord * event = mcrec->event;
  genie::Interaction * in = event->Summary();

  // Get truth stuff out of GENIE ghep record
  sr.nuPDG = in->InitState().ProbePdg();
  sr.nuPDGunosc = in->InitState().ProbePdg(); // fill this for similarity with FD, but no oscillations
  sr.mode = in->ProcInfo().ScatteringTypeId();
  sr.Ev = in->InitState().ProbeE(genie::kRfLab);
  sr.LepPDG = in->FSPrimLeptonPdg();
  sr.isCC = (abs(sr.LepPDG) == 13 || abs(sr.LepPDG) == 11);

  TLorentzVector lepP4;
  TLorentzVector nuP4nuc = *(in->InitState().GetProbeP4(genie::kRfHitNucRest));
  TLorentzVector nuP4 = *(in->InitState().GetProbeP4(genie::kRfLab));

  sr.vtx_x =evt->x;
  sr.vtx_y =evt->y;
  sr.vtx_z =evt->z;
  sr.det_x =0;

  sr.nP = 0;
  sr.nN = 0;
  sr.nipip = 0;
  sr.nipim = 0;
  sr.nipi0 = 0;
  sr.nikp = 0;
  sr.nikm = 0;
  sr.nik0 = 0;
  sr.niem = 0;
  sr.niother = 0;
  sr.nNucleus = 0;
  sr.nUNKNOWN = 0; // there is an "other" category so this never gets used
  sr.eP = 0.;
  sr.eN = 0.;
  sr.ePip = 0.;
  sr.ePim = 0.;
  sr.ePi0 = 0.;
  sr.eOther = 0.;
  sr.eRecoP = 0.;
  sr.eRecoN = 0.;
  sr.eRecoPip = 0.;
  sr.eRecoPim = 0.;
  sr.eRecoPi0 = 0.;
  sr.eOther = 0.;
  

  std::vector<particle> particle_event=evt->particles;
 
  for (std::vector<particle>::iterator it = particle_event.begin() ; it != particle_event.end(); ++it){
   
   //std::cout<<"pdg particella "<<(*it).pdg<<std::endl; 
   double ke = 0.001*((*it).Etrue-(*it).mass);
    if( (*it).pdg == sr.LepPDG ) {
      lepP4.SetPxPyPzE( (*it).pxtrue*0.001,  (*it).pytrue*0.001,  (*it).pztrue*0.001,  (*it).Etrue*0.001); 
      sr.LepE = (*it).Etrue*0.001;
   }
     if( (*it).pdg == 2212 ) {sr.nP++; sr.eP += ke;}
     if( (*it).pdg == 2112 ) {sr.nN++; sr.eN += ke;}
     if( (*it).pdg ==  211 ) {sr.nipip++; sr.ePip += ke;}
     if( (*it).pdg == -211 ) {sr.nipim++; sr.ePim += ke;}
     if( (*it).pdg ==  111 ) {sr.nipi0++; sr.ePi0 += ke;}
     if( (*it).pdg ==  321 ) {sr.nikp++; sr.eOther += ke;}
     if( (*it).pdg == -321 ) {sr.nikm++; sr.eOther += ke;}
     if( (*it).pdg == 311 ||  (*it).pdg == -311 ||  (*it).pdg == 130 ||  (*it).pdg == 310 ) {sr.nik0++; sr.eOther += ke;}
     if( (*it).pdg ==   22 ) {sr.niem++; sr.eOther += ke;}
     if( (*it).pdg > 1000000000 ) sr.nNucleus++;
     else {sr.niother++; sr.eOther += ke;}
   }


/*
  if(ndlarcaf==true){
  for( int i = 0; i < dt.nFS; ++i ) {
    double ke = 0.001*(dt.fsE[i] - sqrt(dt.fsE[i]*dt.fsE[i] - dt.fsPx[i]*dt.fsPx[i] - dt.fsPy[i]*dt.fsPy[i] - dt.fsPz[i]*dt.fsPz[i]));
    if( dt.fsPdg[i] == sr.LepPDG ) {
      lepP4.SetPxPyPzE( dt.fsPx[i]*0.001, dt.fsPy[i]*0.001, dt.fsPz[i]*0.001, dt.fsE[i]*0.001 );
      sr.LepE = dt.fsE[i]*0.001;
    }
    else if( dt.fsPdg[i] == 2212 ) {sr.nP++; sr.eP += ke;}
    else if( dt.fsPdg[i] == 2112 ) {sr.nN++; sr.eN += ke;}
    else if( dt.fsPdg[i] ==  211 ) {sr.nipip++; sr.ePip += ke;}
    else if( dt.fsPdg[i] == -211 ) {sr.nipim++; sr.ePim += ke;}
    else if( dt.fsPdg[i] ==  111 ) {sr.nipi0++; sr.ePi0 += ke;}
    else if( dt.fsPdg[i] ==  321 ) {sr.nikp++; sr.eOther += ke;}
    else if( dt.fsPdg[i] == -321 ) {sr.nikm++; sr.eOther += ke;}
    else if( dt.fsPdg[i] == 311 || dt.fsPdg[i] == -311 || dt.fsPdg[i] == 130 || dt.fsPdg[i] == 310 ) {sr.nik0++; sr.eOther += ke;}
    else if( dt.fsPdg[i] ==   22 ) {sr.niem++; sr.eOther += ke;}
    else if( dt.fsPdg[i] > 1000000000 ) sr.nNucleus++;
    else {sr.niother++; sr.eOther += ke;}
   }
 }
*/
  // true 4-momentum transfer
  TLorentzVector q = nuP4-lepP4;

  // Q2, W, x, y frequently do not get filled in GENIE Kinematics object, so calculate manually
  sr.Q2 = -q.Mag2();
  sr.W = sqrt(0.939*0.939 + 2.*q.E()*0.939 + q.Mag2()); // "Wexp"
  sr.X = -q.Mag2()/(2*0.939*q.E());
  sr.Y = q.E()/sr.Ev;

  sr.theta_reco = -1.; // default value

  sr.NuMomX = nuP4.X();
  sr.NuMomY = nuP4.Y();
  sr.NuMomZ = nuP4.Z();
  sr.LepMomX = lepP4.X();
  sr.LepMomY = lepP4.Y();
  sr.LepMomZ = lepP4.Z();
  sr.LepE = lepP4.E();
  sr.LepNuAngle = nuP4.Angle( lepP4.Vect() );

/*
  if(ndlarcaf==true){
  // todo: come back and make this work for electrons too. what about NC?
    if (abs(sr.LepPDG) == 13)
      sr.LepEndpoint = {dt.muon_endpoint[0], dt.muon_endpoint[1], dt.muon_endpoint[2]};
  }
*/
  
//FIXME  da fare ????   endpoint dalla MCTruth? dalla traiettoria nel rivelatore??
 //  if (abs(sr.LepPDG) == 13)
     // sr.LepEndpoint = {dt.muon_endpoint[0], dt.muon_endpoint[1], dt.muon_endpoint[2]};
//  }

  // Add DUNErw weights to the CAF
  sr.total_xsSyst_cv_wgt = 1;
  // fixme: the following is disabled until DIRT-II finishes on model + uncertainty decisions
  //systtools::event_unit_response_w_cv_t resp = rh.GetEventVariationAndCVResponse(*event);
  //for( const systtools::VarAndCVResponse& it : resp ) {
  //  // Need begin/end to convert double to float
  //  sr.xsSyst_wgt.emplace_back(it.responses.begin(), it.responses.end());
  //  sr.cvwgt.push_back(it.CV_response);
  //  sr.total_xsSyst_cv_wgt *= it.CV_response;
  //}
}
