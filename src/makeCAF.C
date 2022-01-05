#include <cstdio>

#include "TRandom3.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TTree.h"

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "CAF.h"
#include "dumpTree.h"
#include "Params.h"
#include "reco/MLNDLArRecoBranchFiller.h"
#include "reco/ParameterizedRecoBranchFiller.h"
#include "reco/TMSRecoBranchFiller.h"
#include "reco/NDLArTMSMatchRecoFiller.h"
#include "truth/FillTruth.h"


// -------------------------------------------------
fhicl::Table<cafmaker::FhiclConfig> parseConfig(const std::string & configFile)
{
  fhicl::ParameterSet pset;
  cet::filepath_first_absolute_or_lookup_with_dot maker(getenv("FHICL_FILE_PATH"));
  fhicl::make_ParameterSet(configFile, maker, pset);

  // note that this insists the top-level config be named "nd_cafmaker"
  return fhicl::Table<cafmaker::FhiclConfig>{pset.get<fhicl::ParameterSet>("nd_cafmaker"), std::set<std::string>{}};  // second param is 'ignorable keys' -- we want everything validated, so, empty
}

// -------------------------------------------------
// decide which reco fillers we need based on the configuration
std::vector<std::unique_ptr<cafmaker::IRecoBranchFiller>> getRecoFillers(const cafmaker::Params &par)
{
  std::vector<std::unique_ptr<cafmaker::IRecoBranchFiller>> recoFillers;

  // first: did we do ND-LAr reco?
  std::string ndlarFile;
  if (par().cafmaker().ndlarRecoFile(ndlarFile))
    recoFillers.emplace_back(std::make_unique<cafmaker::MLNDLArRecoBranchFiller>(ndlarFile));
  else
  {
    // use the run+subrun numbers to seed the random number generator if seed is not explicitly provided
    // (yes, leak this pointer.  there's only one and sending it back to the main() is annoying)
    auto rando = new TRandom3(par().cafmaker().seed() >= 0 ?
                              par().cafmaker().seed() :
                              par().runInfo().run() * 1000 + par().runInfo().subrun());

    recoFillers.emplace_back(std::make_unique<cafmaker::ParameterizedRecoBranchFiller>(rando));
  }

  // next: did we do TMS reco?
  std::string tmsFile;
  if (par().cafmaker().tmsRecoFile(tmsFile))
    recoFillers.emplace_back(std::make_unique<cafmaker::TMSRecoBranchFiller>(ndlarFile));

  // if we did both ND-LAr and TMS, we should try to match them, too
  if (!ndlarFile.empty() && !tmsFile.empty())
     recoFillers.emplace_back(std::make_unique<cafmaker::NDLArTMSMatchRecoFiller>());

  return recoFillers;
}

// main loop function
void loop(CAF& caf,
          cafmaker::Params &par,
          TTree * intree,
          TTree * gtree,
          const std::vector<std::unique_ptr<cafmaker::IRecoBranchFiller>> & recoFillers)
{
  //// read in dumpTree output file
  cafmaker::dumpTree dt;
  dt.BindToTree(intree);

  caf.pot = gtree->GetWeight();
  gtree->SetBranchAddress( "gmcrec", &caf.mcrec );

  // Main event loop
  int N = intree->GetEntries();
  for( int ii = par().cafmaker().first(); ii < N; ++ii ) {

    intree->GetEntry(ii);
    if( ii % 100 == 0 ) printf( "Event %d of %d...\n", ii, N );

    // reset (the default constructor initializes its variables)
    caf.setToBS();

    caf.sr.run = par().runInfo().run();
    caf.sr.subrun = par().runInfo().subrun();
    caf.sr.event = ii;
    caf.sr.isFD = 0;
    caf.sr.isFHC = par().runInfo().fhc();

    // in the future this can be extended to use 'truth fillers'
    // (like the reco ones) if we find that the truth filling
    // is getting too complex for one function
    fillTruth(caf.sr, dt, gtree, caf.mcrec, par, caf.rh);

    // hand off to the correct reco filler(s).
    for (const auto & filler : recoFillers)
      filler->FillRecoBranches(ii, caf.sr, dt, par);

    caf.fill();
  }

  // set other metadata
  caf.meta_run = par().runInfo().run();
  caf.meta_subrun = par().runInfo().subrun();

}

// -------------------------------------------------

int main( int argc, char const *argv[] ) 
{

  if( (argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1])) ) {
    std::cout << "Help yourself by looking at the source code to see what the options are." << std::endl;
    return 0;
  }

  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <fhicl configuration file>" << std::endl;
    return 1;
  }

  // Need this to store event-by-event geometric efficiency
  gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t> > >", "vector");

  cafmaker::Params par = parseConfig(argv[1]);

  CAF caf( par().cafmaker().outputFile(), par().cafmaker().nusystsFcl() );

  TFile * tf = new TFile( par().cafmaker().dumpFile().c_str() );
  TTree * tree = (TTree*) tf->Get( "tree" );

  TFile * gf = new TFile( par().cafmaker().ghepFile().c_str() );
  TTree * gtree = (TTree*) gf->Get( "gtree" );

  loop( caf, par, tree, gtree, getRecoFillers(par) );

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
