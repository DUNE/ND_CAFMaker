#include <cstdio>

#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/positional_options.hpp"
#include "boost/program_options/variables_map.hpp"

#include "TRandom3.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TTree.h"

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/parse.h"

#include "CAF.h"
#include "dumpTree.h"
#include "Params.h"
#include "reco/MLNDLArRecoBranchFiller.h"
#include "reco/ParameterizedRecoBranchFiller.h"
#include "reco/TMSRecoBranchFiller.h"
#include "reco/NDLArTMSMatchRecoFiller.h"
#include "truth/FillTruth.h"

namespace progopt = boost::program_options;

// -------------------------------------------------
progopt::variables_map parseCmdLine(int argc, const char** argv)
{
  progopt::options_description genopts("General options");
  genopts.add_options()
      ("help,h", "print this help message");

  progopt::options_description fclopts("FCL overrides (for quick tests; edit your .fcl for regular usage)");
  fclopts.add_options()
      ("dump,d",     progopt::value<std::string>(), "input 'dump' file from dumpTree.py")
      ("ghep,g",     progopt::value<std::string>(), "input GENIE .ghep file")
      ("out,o",      progopt::value<std::string>(), "output CAF file")
      ("startevt",   progopt::value<int>(),         "event number to start at")
      ("numevts,n",  progopt::value<int>(),         "total number of events to process (-1 means 'all')");

  // this option needs to exist for the positional argument be assigned to it,
  // but we don't want to show it in the '--help' printout
  progopt::options_description hidden("hidden options");
  hidden.add_options()
      ("fcl",        progopt::value<std::string>(), "driver FCL");

  progopt::positional_options_description pos;
  pos.add("fcl", 1);

  progopt::options_description allopts;
  allopts.add(genopts).add(fclopts).add(hidden);
  progopt::variables_map vm;
  progopt::store(progopt::command_line_parser(argc, argv).options(allopts).positional(pos).run(), vm);
  progopt::notify(vm);

  if (vm.count("help"))
  {
    progopt::options_description opts;
    opts.add(genopts).add(fclopts);
    std::cout << "Usage: " << argv[0] << " [options] <driver.fcl> [options]" << std::endl;
    std::cout << opts << std::endl;
    exit(0);
  }

  return vm;
}

// -------------------------------------------------
fhicl::Table<cafmaker::FhiclConfig> parseConfig(const std::string & configFile, const progopt::variables_map & vm)
{
  fhicl::intermediate_table provisional;
  cet::filepath_first_absolute_or_lookup_with_dot maker(getenv("FHICL_FILE_PATH"));
  fhicl::parse_document(configFile, maker, provisional);

  if (!provisional.exists("nd_cafmaker"))
  {
    std::cerr << "Ill-formed FHICL in '" << configFile << "':" << std::endl;
    std::cerr << "Outermost table must be named 'nd_cafmaker'." << std::endl;
    exit(1);
  }

  // insert anything overridden on the command line here
  if (vm.count("dump"))
    provisional.put("nd_cafmaker.CAFMakerSettings.InputDumpFile", vm["dump"].as<std::string>());
  if (vm.count("ghep"))
    provisional.put("nd_cafmaker.CAFMakerSettings.InputGHEPFile", vm["ghep"].as<std::string>());
  if (vm.count("out"))
    provisional.put("nd_cafmaker.CAFMakerSettings.OutputFile", vm["out"].as<std::string>());

  if (vm.count("startevt"))
    provisional.put("nd_cafmaker.CAFMakerSettings.FirstEvt", vm["startevt"].as<int>());
  if (vm.count("numevts"))
    provisional.put("nd_cafmaker.CAFMakerSettings.NumEvts", vm["numevts"].as<int>());

  // now that we've updated it, convert to actual ParameterSet
  fhicl::ParameterSet pset;
  fhicl::make_ParameterSet(provisional, pset);

  // finally, convert to a table, which does the validation.
  // note that this usage insists the top-level config be named "nd_cafmaker".
  auto params = pset.get<fhicl::ParameterSet>("nd_cafmaker");
  // the second param is 'ignorable keys' -- we want everything validated, so, empty
  fhicl::Table<cafmaker::FhiclConfig> table{params, std::set<std::string>{}};

  return table;
}

// -------------------------------------------------
// update the FCL table with anything from the command line
void updateConfig(fhicl::Table<cafmaker::FhiclConfig> & table, const progopt::variables_map & vm)
{
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

// -------------------------------------------------
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
  int N = par().cafmaker().numevts() > 0 ? par().cafmaker().numevts() : intree->GetEntries() - par().cafmaker().first();
  int start = par().cafmaker().first();
  for( int ii = start; ii < start + N; ++ii ) {

    intree->GetEntry(ii);
    if( ii % 100 == 0 ) printf( "Event %d (%d of %d)...\n", ii, ii-start, N );

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

  progopt::variables_map vars = parseCmdLine(argc, argv);

  // Need this to store event-by-event geometric efficiency
  gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t> > >", "vector");

  cafmaker::Params par = parseConfig(vars["fcl"].as<std::string>(), vars);

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
