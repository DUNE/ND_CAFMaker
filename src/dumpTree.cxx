#include "dumpTree.h"

#include "TTree.h"

void cafmaker::dumpTree::BindToTree(TTree *tree)
{
  tree->SetBranchAddress( "ievt", &ievt );
  tree->SetBranchAddress( "lepPdg", &lepPdg );
  tree->SetBranchAddress( "muonReco", &muonReco );
  tree->SetBranchAddress( "lepKE", &lepKE );
  tree->SetBranchAddress( "muGArLen", &muGArLen );
  tree->SetBranchAddress( "muECalLen", &muECalLen );
  tree->SetBranchAddress( "hadTot", &hadTot );
  tree->SetBranchAddress( "hadCollar", &hadCollar );
  tree->SetBranchAddress( "hadP", &hadP );
  tree->SetBranchAddress( "hadN", &hadN );
  tree->SetBranchAddress( "hadPip", &hadPip );
  tree->SetBranchAddress( "hadPim", &hadPim );
  tree->SetBranchAddress( "hadPi0", &hadPi0 );
  tree->SetBranchAddress( "hadOther", &hadOther );
  tree->SetBranchAddress( "p3lep", p3lep );
  tree->SetBranchAddress( "vtx", vtx );
  tree->SetBranchAddress( "muonExitPt", muonExitPt );
  tree->SetBranchAddress( "muonExitMom", muonExitMom );
  tree->SetBranchAddress( "lepDeath", muon_endpoint );
  tree->SetBranchAddress( "muon_endVolName", &muon_endVolName );
  tree->SetBranchAddress( "nFS", &nFS );
  tree->SetBranchAddress( "fsPdg", fsPdg );
  tree->SetBranchAddress( "fsPx", fsPx );
  tree->SetBranchAddress( "fsPy", fsPy );
  tree->SetBranchAddress( "fsPz", fsPz );
  tree->SetBranchAddress( "fsE", fsE );
  tree->SetBranchAddress( "fsTrkLen", fsTrkLen );
//  tree->SetBranchAddress( "fsTrkLenPerp", fsTrkLenPerp );

  tree->SetBranchAddress( "geoEffThrowResults", &geoEffThrowResults );
}
