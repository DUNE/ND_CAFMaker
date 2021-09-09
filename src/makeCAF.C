#include <stdio.h>

#include "CAF.h"

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "EVGCore/EventRecord.h"
#include "nusystematics/artless/response_helper.hh"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "dumpTree.h"
#include "Params.h"
#include "reco/ParameterizedRecoBranchFiller.h"
#include "reco/IRecoBranchFiller.h"


// Fill truth info
void fillTruth(caf::StandardRecord& sr,
               const cafmaker::dumpTree & dt,
               TTree * gtree,
               const genie::NtpMCEventRecord * mcrec,
               const cafmaker::params &par,
               nusyst::response_helper& rh)
{
  sr.vtx_x = dt.vtx[0];
  sr.vtx_y = dt.vtx[1];
  sr.vtx_z = dt.vtx[2];
  sr.det_x = -100.*par.OA_xcoord;

  // get GENIE event record
  gtree->GetEntry( dt.ievt );
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

  // Add DUNErw weights to the CAF
  sr.total_xsSyst_cv_wgt = 1;
  systtools::event_unit_response_w_cv_t resp = rh.GetEventVariationAndCVResponse(*event);
  for( const systtools::VarAndCVResponse& it : resp ) {
    // Need begin/end to convert double to float
    sr.xsSyst_wgt.emplace_back(it.responses.begin(), it.responses.end());
    sr.cvwgt.push_back(it.CV_response);
    sr.total_xsSyst_cv_wgt *= it.CV_response;
  }

} // void fillTruth()


// -------------------------------------------------

// main loop function
void loop(CAF& caf, cafmaker::params &par,
          TTree * intree,
          TTree * gtree,
          const cafmaker::IRecoBranchFiller & recoFiller)
{
  //// read in dumpTree output file
  cafmaker::dumpTree dt;
  dt.BindToTree(intree);

  caf.pot = gtree->GetWeight();
  gtree->SetBranchAddress( "gmcrec", &caf.mcrec );

  // Main event loop
  int N = intree->GetEntries();
  for( int ii = par.first; ii < N; ++ii ) {

    intree->GetEntry(ii);
    if( ii % 100 == 0 ) printf( "Event %d of %d...\n", ii, N );

    // reset (the default constructor initializes its variables)
    caf.setToBS();

    caf.sr.run = par.run;
    caf.sr.subrun = par.subrun;
    caf.sr.event = ii;
    caf.sr.isFD = 0;
    caf.sr.isFHC = par.fhc;

    fillTruth(caf.sr, dt, gtree, caf.mcrec, par, caf.rh);
    recoFiller.FillRecoBranches(caf.sr, dt, par);

    caf.fill();
  }

  // set other metadata
  caf.meta_run = par.run;
  caf.meta_subrun = par.subrun;

}

// -------------------------------------------------

int main( int argc, char const *argv[] ) 
{

  if( (argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1])) ) {
    std::cout << "Help yourself by looking at the source code to see what the options are." << std::endl;
    return 0;
  }

  // Need this to store event-by-event geometric efficiency
  gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t> > >", "vector");

  // get command line options
  std::string gfile;
  std::string infile;
  std::string outfile;
  std::string fhicl_filename;

  // Make parameter object and set defaults
  cafmaker::params par;
  par.IsGasTPC = false;
  par.OA_xcoord = 0.; // on-axis by default
  par.fhc = true;
  par.grid = false;
  par.seed = 7; // a very random number
  par.run = 1; // CAFAna doesn't like run number 0
  par.subrun = 0;
  par.first = 0;
  par.trk_muRes = 0.02; // fractional muon energy resolution of HP GAr TPC
  par.LAr_muRes = 0.05; // fractional muon energy resolution of muons contained in LAr
  par.ECAL_muRes = 0.1; // fractional muon energy resolution of muons ending in ECAL
  par.em_const = 0.03; // EM energy resolution constant term: A + B/sqrt(E) (GeV)
  par.em_sqrtE = 0.1; // EM energy resolution 1/sqrt(E) term: A + B/sqrt(E) (GeV)
  par.michelEff = 0.75; // Michel finder efficiency
  par.CC_trk_length = 100.; // minimum track length for CC in cm
  par.pileup_frac = 0.1; // fraction of events with non-zero pile-up
  par.pileup_max = 0.5; // GeV
  par.gastpc_len = 6.; // track length cut in cm
  par.gastpc_B = 0.4; // B field strength in Tesla
  par.gastpc_padPitch = 0.1; // 1 mm. Actual pad pitch varies, which is going to be impossible to implement
  par.gastpc_X0 = 1300.; // cm = 13m radiation length

  int i = 0;
  while( i < argc ) {
    if( argv[i] == std::string("--infile") ) {
      infile = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--gfile") ) {
      gfile = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--outfile") ) {
      outfile = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--fhicl") ) {
      fhicl_filename = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--seed") ) {
      par.seed = atoi(argv[i+1]);
      par.run = par.seed;
      i += 2;
    } else if( argv[i] == std::string("--oa") ) {
      par.OA_xcoord = atof(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--rhc") ) {
      par.fhc = false;
      i += 1;
    } else if( argv[i] == std::string("--gastpc") ) {
      par.IsGasTPC = true;
      i += 1;
    } else i += 1; // look for next thing
  }

  CAF caf( outfile, fhicl_filename );

  TFile * tf = new TFile( infile.c_str() );
  TTree * tree = (TTree*) tf->Get( "tree" );

  TFile * gf = new TFile( gfile.c_str() );
  TTree * gtree = (TTree*) gf->Get( "gtree" );

  // for the moment this is the only reco filler we have,
  // but there are others waiting in the wings...
  TRandom3 rando;
  cafmaker::ParameterizedRecoBranchFiller recoFiller(&rando);
  loop( caf, par, tree, gtree, recoFiller );

  caf.version = 4;
  printf( "Run %d POT %g\n", caf.meta_run, caf.pot );
  caf.fillPOT();

  // Copy geometric efficiency throws TTree to CAF file
  std::cout << "Copying geometric efficiency throws TTree to output file" << std::endl;
  TTree *tGeoEfficiencyThrowsOut = (TTree*) tf->Get("geoEffThrows");
  caf.cafFile->cd();
  tGeoEfficiencyThrowsOut->CloneTree()->Write();

  std::cout << "Writing CAF" << std::endl;
  caf.write();
    
}
