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
#include "Params.h"
#include "reco/MLNDLArRecoBranchFiller.h"
#include "reco/TMSRecoBranchFiller.h"
#include "reco/NDLArTMSMatchRecoFiller.h"
#include "reco/SANDRecoBranchFiller.h"
#include "truth/FillTruth.h"
#include "duneanaobj/StandardRecord/SREnums.h"

namespace progopt = boost::program_options;

// -------------------------------------------------
progopt::variables_map parseCmdLine(int argc, const char** argv)
{
  progopt::options_description genopts("General options");
  genopts.add_options()
      ("help,h", "print this help message");

  progopt::options_description fclopts("FCL overrides (for quick tests; edit your .fcl for regular usage)");
  fclopts.add_options()
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
  cet::filepath_first_absolute_or_lookup_with_dot maker(getenv("FHICL_FILE_PATH"));
  fhicl::intermediate_table provisional = fhicl::parse_document(configFile, maker);

  if (!provisional.exists("nd_cafmaker"))
  {
    std::cerr << "Ill-formed FHICL in '" << configFile << "':" << std::endl;
    std::cerr << "Outermost table must be named 'nd_cafmaker'." << std::endl;
    exit(1);
  }

  // insert anything overridden on the command line here
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
// decide which reco fillers we need based on the configuration
std::vector<std::unique_ptr<cafmaker::IRecoBranchFiller>> getRecoFillers(const cafmaker::Params &par)
{
  std::vector<std::unique_ptr<cafmaker::IRecoBranchFiller>> recoFillers;

  std::cout << "Filling Reco info for the following cases:\n";

  // first: we do SAND or ND-LAr reco
  std::string ndlarFile;
  std::string sandFile;
  if (par().cafmaker().ndlarRecoFile(ndlarFile))
  {
    recoFillers.emplace_back(std::make_unique<cafmaker::MLNDLArRecoBranchFiller>(ndlarFile));
    std::cout << "   ND-LAr (Deep-Learn-Physics ML)\n";
  } else if (par().cafmaker().sandRecoFile(sandFile))
  {
    recoFillers.emplace_back(std::make_unique<cafmaker::SANDRecoBranchFiller>(sandFile));
    std::cout << "   SAND\n";
  }

  // next: did we do TMS reco?
  std::string tmsFile;
  if (par().cafmaker().tmsRecoFile(tmsFile))
  {
    recoFillers.emplace_back(std::make_unique<cafmaker::TMSRecoBranchFiller>(tmsFile));
    std::cout << "   TMS\n";
  }

  // if we did both ND-LAr and TMS, we should try to match them, too
  if (!ndlarFile.empty() && !tmsFile.empty())
  {
    recoFillers.emplace_back(std::make_unique<cafmaker::NDLArTMSMatchRecoFiller>());
    std::cout << "   ND-LAr + TMS matching\n";
  }


  return recoFillers;
}

// -------------------------------------------------
// main loop function
void loop(CAF& caf,
          cafmaker::Params &par,
          TTree * gtree,
          const std::vector<std::unique_ptr<cafmaker::IRecoBranchFiller>> & recoFillers)
{

  // Enable ND_LAr detector
  if (gtree)
  {
    caf.pot = gtree->GetWeight();
    gtree->SetBranchAddress("gmcrec", &caf.mcrec);
  }

  // if this is a data file, there won't be any truth, of course
  std::unique_ptr<cafmaker::TruthMatcher> truthMatcher;
  if (gtree)
    truthMatcher = std::make_unique<cafmaker::TruthMatcher>(gtree, caf.mcrec);


  // Main event loop
  int N = par().cafmaker().numevts() > 0 ? par().cafmaker().numevts() : gtree->GetEntries() - par().cafmaker().first();
  int start = par().cafmaker().first();

  for( int ii = start; ii < start + N; ++ii ) {

    if( ii % 10000 == 0 )
      printf( "Event %d (%d of %d)...\n", ii, (ii-start)+1, N );

    // reset (the default constructor initializes its variables)
    caf.setToBS();

   //old SR variables

  //  caf.sr.meta_run = par().runInfo().run();
   // caf.sr.meta_subrun = par().runInfo().subrun();
 //   caf.sr.isFD = 0;
  //  caf.sr.isFHC = par().runInfo().fhc();
    caf.sr.beam.pulsepot = caf.pot;

    // hand off to the correct reco filler(s).
    for (const auto & filler : recoFillers)
      filler->FillRecoBranches(ii, caf.sr, par, truthMatcher.get());

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

  cafmaker::Params par = parseConfig(vars["fcl"].as<std::string>(), vars);

  CAF caf(par().cafmaker().outputFile(), par().cafmaker().nusystsFcl(), par().cafmaker().makeFlatCAF());

  TFile * gf = nullptr;
  TTree * gtree = nullptr;
  if (!par().cafmaker().ghepFile().empty())
  {
    gf = new TFile(par().cafmaker().ghepFile().c_str());   //reading genie file
    gtree = (TTree *) gf->Get("gtree");
  }

  loop( caf, par, gtree, getRecoFillers(par) );

  caf.version = 4;
  printf( "Run %d POT %g\n", caf.meta_run, caf.pot );
  caf.fillPOT();

  std::cout << "Writing CAF" << std::endl;
  caf.write();

  return 0;
}
