#define CAF_cxx
#ifdef CAF_cxx

#include "CAF.h"

CAF::CAF( std::string filename, bool isGas )
{
  cafFile = new TFile( filename.c_str(), "RECREATE" );
  cafMVA = new TTree( "caf", "caf" );
  cafPOT = new TTree( "meta", "meta" );
  genie = new TTree( "genieEvt", "genieEvt" );

  // initialize the GENIE record
  mcrec = NULL;

  cafMVA->Branch( "run", &run, "run/I" );
  cafMVA->Branch( "subrun", &subrun, "subrun/I" );
  cafMVA->Branch( "event", &event, "event/I" );
  cafMVA->Branch( "isFD", &isFD, "isFD/I" );
  cafMVA->Branch( "isFHC", &isFHC, "isFHC/I" );
  cafMVA->Branch( "isCC", &isCC, "isCC/I" );

  cafMVA->Branch( "nuPDG", &neutrinoPDG, "nuPDG/I" );
  cafMVA->Branch( "nuPDGunosc", &neutrinoPDGunosc, "nuPDGunosc/I");
  cafMVA->Branch( "NuMomX", &NuMomX, "NuMomX/D" );
  cafMVA->Branch( "NuMomY", &NuMomY, "NuMomY/D" );
  cafMVA->Branch( "NuMomZ", &NuMomZ, "NuMomZ/D" );
  cafMVA->Branch( "Ev", &Ev, "Ev/D" );
  cafMVA->Branch( "mode", &mode, "mode/I" );
  cafMVA->Branch( "LepPDG", &LepPDG, "LepPDG/I" );
  cafMVA->Branch( "LepMomX", &LepMomX, "LepMomX/D" );
  cafMVA->Branch( "LepMomY", &LepMomY, "LepMomY/D" );
  cafMVA->Branch( "LepMomZ", &LepMomZ, "LepMomZ/D" );
  cafMVA->Branch( "LepE", &LepE, "LepE/D" );
  cafMVA->Branch( "LepNuAngle", &LepNuAngle, "LepNuAngle/D" );
  cafMVA->Branch( "Q2", &Q2, "Q2/D" );
  cafMVA->Branch( "W", &W, "W/D" );
  cafMVA->Branch( "X", &X, "X/D" );
  cafMVA->Branch( "Y", &Y, "Y/D" );

  cafMVA->Branch( "nP", &nP, "nP/I" );
  cafMVA->Branch( "nN", &nN, "nN/I" );
  cafMVA->Branch( "nipip", &nipip, "nipip/I" );
  cafMVA->Branch( "nipim", &nipim, "nipim/I" );
  cafMVA->Branch( "nipi0", &nipi0, "nipi0/I" );
  cafMVA->Branch( "nikp", &nikp, "nikp/I" );
  cafMVA->Branch( "nikm", &nikm, "nikm/I" );
  cafMVA->Branch( "nik0", &nik0, "nik0/I" );
  cafMVA->Branch( "niem", &niem, "niem/I" );
  cafMVA->Branch( "niother", &niother, "niother/I" );
  cafMVA->Branch( "nNucleus", &nNucleus, "nNucleus/I" );
  cafMVA->Branch( "nUNKNOWN", &nUNKNOWN, "nUNKNOWN/I" );

  cafMVA->Branch("eP",        &eP,         "eP/D");
  cafMVA->Branch("eN",        &eN,         "eN/D");
  cafMVA->Branch("ePip",      &ePip,       "ePip/D");
  cafMVA->Branch("ePim",      &ePim,       "ePim/D");
  cafMVA->Branch("ePi0",      &ePi0,       "ePi0/D");
  cafMVA->Branch("eOther",    &eOther,     "eOther/D");
  cafMVA->Branch("eRecoP",        &eRecoP,         "eRecoP/D");
  cafMVA->Branch("eRecoN",        &eRecoN,         "eRecoN/D");
  cafMVA->Branch("eRecoPip",      &eRecoPip,       "eRecoPip/D");
  cafMVA->Branch("eRecoPim",      &eRecoPim,       "eRecoPim/D");
  cafMVA->Branch("eRecoPi0",      &eRecoPi0,       "eRecoPi0/D");
  cafMVA->Branch("eRecoOther",    &eRecoOther,     "eRecoOther/D");

  cafMVA->Branch( "det_x", &det_x, "det_x/D" );
  cafMVA->Branch( "vtx_x", &vtx_x, "vtx_x/D" );
  cafMVA->Branch( "vtx_y", &vtx_y, "vtx_y/D" );
  cafMVA->Branch( "vtx_z", &vtx_z, "vtx_z/D" );

  cafMVA->Branch( "Ev_reco", &Ev_reco, "Ev_reco/D" );
  cafMVA->Branch( "Elep_reco", &Elep_reco, "Elep_reco/D" );
  cafMVA->Branch( "theta_reco", &theta_reco, "theta_reco/D" );
  cafMVA->Branch( "reco_numu", &reco_numu, "reco_numu/I" );
  cafMVA->Branch( "reco_nue", &reco_nue, "reco_nue/I" );
  cafMVA->Branch( "reco_nc", &reco_nc, "reco_nc/I" );
  cafMVA->Branch( "reco_q", &reco_q, "reco_q/I" );
  cafMVA->Branch( "muon_contained", &muon_contained, "muon_contained/I" );
  cafMVA->Branch( "muon_tracker", &muon_tracker, "muon_tracker/I" );
  cafMVA->Branch( "muon_ecal", &muon_ecal, "muon_ecal/I" );
  cafMVA->Branch( "muon_exit", &muon_exit, "muon_exit/I" );
  cafMVA->Branch( "reco_lepton_pdg", &reco_lepton_pdg, "reco_lepton_pdg/I" );
  cafMVA->Branch( "Ehad_veto", &Ehad_veto, "Ehad_veto/D" );
  cafMVA->Branch( "pileup_energy", &pileup_energy, "pileup_energy/D" );

  if( isGas ) {
    cafMVA->Branch( "gastpc_pi_pl_mult", &gastpc_pi_pl_mult, "gastpc_pi_pl_mult/I" );
    cafMVA->Branch( "gastpc_pi_min_mult", &gastpc_pi_min_mult, "gastpc_pi_min_mult/I" );
    cafMVA->Branch( "nFSP", &nFSP, "nFSP/I" );
    cafMVA->Branch( "pdg", pdg, "pdg[nFSP]/I" );    
    cafMVA->Branch( "ptrue", ptrue, "ptrue[nFSP]/D" );    
    cafMVA->Branch( "trkLen", trkLen, "trkLen[nFSP]/D" );    
    cafMVA->Branch( "trkLenPerp", trkLenPerp, "trkLenPerp[nFSP]/D" );    
    cafMVA->Branch( "partEvReco", partEvReco, "partEvReco[nFSP]/D" );    
  }

  genie->Branch( "genie_record", &mcrec );

  cafPOT->Branch( "pot", &pot, "pot/D" );
  cafPOT->Branch( "run", &meta_run, "run/I" );
  cafPOT->Branch( "subrun", &meta_subrun, "subrun/I" );
  cafPOT->Branch( "version", &version, "version/I" );
}

CAF::~CAF() {}

void CAF::fill()
{
  cafMVA->Fill();
  genie->Fill();
}

void CAF::Print()
{
  printf( "Event %d:\n", event );
  printf( "   Truth: Ev = %3.3f Elep = %3.3f Q2 = %3.3f W = %3.3f x = %3.3f y = %3.3f lepton %d mode %d\n", Ev, LepE, Q2, W, X, Y, LepPDG, mode );
  printf( "    Reco: Ev = %3.3f Elep = %3.3f q %d mu/e/nc %d%d%d cont/trk/ecal/exit %d%d%d%d had veto %2.1f\n\n", Ev_reco, Elep_reco, reco_q, reco_numu, reco_nue, reco_nc, muon_contained, muon_tracker, muon_ecal, muon_exit, Ehad_veto );
}

void CAF::fillPOT()
{
  printf( "Filling metadata\n" );
  cafPOT->Fill();
}

void CAF::write()
{
  cafFile->cd();
  cafMVA->Write();
  cafPOT->Write();
  genie->Write();
  cafFile->Close();
}

void CAF::addRWbranch( int parId, std::string name, std::string wgt_var, std::vector<double> &vars )
{
  cafMVA->Branch( Form("%s_nshifts", name.c_str()), &nwgt[parId], Form("%s_nshifts/I", name.c_str()) );
  cafMVA->Branch( Form("%s_cv%s", name.c_str(), wgt_var.c_str()), &cvwgt[parId], Form("%s_cv%s/D", name.c_str(), wgt_var.c_str()) );
  cafMVA->Branch( Form("%s_%s", wgt_var.c_str(), name.c_str()), wgt[parId], Form("%s_%s[%s_nshifts]/D", wgt_var.c_str(), name.c_str(), name.c_str()) );
}

void CAF::setToBS()
{
  isFD = -1;
  isFHC = -1;
  run = -1; subrun = -1; event = -1;
  isCC = -1; neutrinoPDG = 0; neutrinoPDGunosc = 0;
  mode = 0; LepPDG = 0;
  Ev = -1.; Q2 = -1.; W = -1.; X = -1.; Y = -1.;
  NuMomX = -999.; NuMomY = -999.; NuMomZ = -999.;
  LepMomX = -999.; LepMomY = -999.; LepMomZ = -999.; 
  LepE = -999.; LepNuAngle = -999.;
  nP = 0; nN = 0; nipip = 0; nipim = 0; nipi0 = 0; nikp = 0; nikm = 0; nik0 = 0; niem = 0; niother = 0; nNucleus = 0; nUNKNOWN = 0;
  eP = 0.; eN = 0.; ePip = 0.; ePim = 0.; ePi0 = 0.; eOther = 0.;
  eRecoP = 0.; eRecoN = 0.; eRecoPip = 0.; eRecoPim = 0.; eRecoPi0 = 0.; eRecoOther = 0.;
  vtx_x = -9999.; vtx_y = -9999.; vtx_z = -9999.;
  det_x = -9999.;
  Ev_reco = 0.; Elep_reco = 0.; theta_reco = 0.;
  reco_numu = 0; reco_nue = 0; reco_nc = 0; reco_q = 0;
  muon_contained = 0; muon_tracker = 0; muon_ecal = 0; muon_exit = 0; reco_lepton_pdg = 0;
  Ehad_veto = 0.;
  pileup_energy = 0.;

  gastpc_pi_pl_mult = 0;
  gastpc_pi_min_mult = 0;
  nFSP = 0;

  // nwgt and iswgt do not change event by event and are set outside the event loop
  // do not reset them to BS here
  for( int i = 0; i < 100; ++i ) {
    if( iswgt[i] ) cvwgt[i] = 1.;
    else cvwgt[i] = 0.;
    for( int j = 0; j < 100; ++j ) {
      if( iswgt[i] ) wgt[i][j] = 1.;
      else wgt[i][j] = 0.;
    }
    trkLen[i] = 0.;
    trkLenPerp[i] = 0.;
    partEvReco[i] = 0.;
    ptrue[i] = 0.;
    pdg[i] = 0;
  }
}

#endif
