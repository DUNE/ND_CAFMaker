#ifndef CAF_h
#define CAF_h

#include "TFile.h"
#include "TTree.h"
#include "Ntuple/NtpMCEventRecord.h"

class CAF {

public:
  CAF( std::string filename, bool isGas = false );
  ~CAF();
  void fill();
  void fillPOT();
  void write();
  void addRWbranch( int parId, std::string name, std::string wgt_var, std::vector<double> &vars );
  void Print();
  void setToBS();

  // Make ntuple variables public so they can be set from other file

  // configuration variables
  int isFD, isFHC;
  // event accounting
  int run, subrun, event;
  // Truth information
  int isCC, neutrinoPDG, neutrinoPDGunosc, mode, LepPDG; 
  double Ev, Q2, W, X, Y, NuMomX, NuMomY, NuMomZ, LepMomX, LepMomY, LepMomZ, LepE, LepNuAngle;
  // True particle counts
  int nP, nN, nipip, nipim, nipi0, nikp, nikm, nik0, niem, niother, nNucleus, nUNKNOWN;
  double eP, eN, ePip, ePim, ePi0, eOther;
  double eRecoP, eRecoN, eRecoPip, eRecoPim, eRecoPi0, eRecoOther;

  // vertex -- smear it?
  double vtx_x, vtx_y, vtx_z;
  double det_x;

  // Reco information CV
  double Ev_reco, Elep_reco, theta_reco;
  int reco_numu, reco_nue, reco_nc, reco_q;
  int muon_contained, muon_tracker, muon_ecal, muon_exit, reco_lepton_pdg;
  float muon_endpoint[3];
  std::string * muon_endVolName;
  double Ehad_veto;
  double lepKE;
  //DEAD_MODULE_UPDATE
  double dead_Ehad_veto;
  //INNER_UPDATE
  double inner_dead_Ehad_veto;
  bool flagr;
  double pileup_energy;

  // Gas TPC variables
  int gastpc_pi_min_mult, gastpc_pi_pl_mult;
  int nFSP;
  int pdg[100];
  double trkLen[100], trkLenPerp[100], ptrue[100], partEvReco[100];

  // reweights -- make sure big enough to hold all the variations for each knob, and all the knobs
  // the names, and what they actually mean, are determined automatically from the fhicl input file
  int nwgt[100];
  double cvwgt[100];
  double wgt[100][100];
  bool iswgt[100];

  // store the GENIE record as a branch
  genie::NtpMCEventRecord * mcrec;

  // Event-by-event geometric efficiency throw results
  std::vector< std::vector < std::vector < uint64_t > > > * geoEffThrowResults;

  // meta
  double pot;
  int meta_run, meta_subrun;
  int version;

  TFile * cafFile;
  TTree * cafMVA;
  TTree * cafPOT;
  TTree * genie;
};

#endif

