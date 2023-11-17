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

// GENIE
#include "Framework/Ntuple/NtpMCEventRecord.h"

#include "CAF.h"
#include "Params.h"
#include "reco/MLNDLArRecoBranchFiller.h"
#include "reco/TMSRecoBranchFiller.h"
#include "reco/NDLArTMSMatchRecoFiller.h"
#include "reco/SANDRecoBranchFiller.h"
#include "reco/MINERvARecoBranchFiller.h"
#include "truth/FillTruth.h"
#include "util/Logger.h"
#include "util/Progress.h"

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
  fhicl::ParameterSet pset = fhicl::ParameterSet::make(provisional);

  // finally, convert to a table, which does the validation.
  // note that this usage insists the top-level config be named "nd_cafmaker".
  auto params = pset.get<fhicl::ParameterSet>("nd_cafmaker");
  // the second param is 'ignorable keys' -- we want everything validated, so, empty
  fhicl::Table<cafmaker::FhiclConfig> table{params, std::set<std::string>{}};

  return table;
}

// -------------------------------------------------
// decide which reco fillers we need based on the configuration
std::vector<std::unique_ptr<cafmaker::IRecoBranchFiller>> getRecoFillers(const cafmaker::Params &par,
                                                                         cafmaker::Logger::THRESHOLD logThresh)
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

  // next: did we do the MINERvA reco?
  std::string minervaFile;
  if (par().cafmaker().minervaRecoFile(minervaFile))
  {
    recoFillers.emplace_back(std::make_unique<cafmaker::MINERvARecoBranchFiller>(minervaFile));
    std::cout<< "   MINERVA\n";
  }

  // if we did both ND-LAr and TMS, we should try to match them, too
  if (!ndlarFile.empty() && !tmsFile.empty())
  {
    recoFillers.emplace_back(std::make_unique<cafmaker::NDLArTMSMatchRecoFiller>());
    std::cout << "   ND-LAr + TMS matching\n";
  }

  // for now all the fillers get the same threshold.
  // if we decide we need to do it differently later
  // we can adjust the FCL params...
  for (std::unique_ptr<cafmaker::IRecoBranchFiller> & filler : recoFillers)
    filler->SetLogThrehsold(logThresh);

  return recoFillers;
}

// -------------------------------------------------
bool doTriggersMatch(const cafmaker::Trigger& t1, const cafmaker::Trigger& t2, unsigned int dT)
{
  return (t1.triggerTime_s == t2.triggerTime_s && abs(int(t1.triggerTime_ns - t2.triggerTime_ns)) < dT);
}

struct triggerTimeCmp
{
  bool operator()(const cafmaker::Trigger &t1, const cafmaker::Trigger &t2)
  {
    return t1.triggerTime_s < t2.triggerTime_s ||
           (t1.triggerTime_s == t2.triggerTime_s && t1.triggerTime_ns < t2.triggerTime_ns);
  }
};

struct triggerTimePtrCmp
{
  triggerTimeCmp cmp;
  bool operator()(const cafmaker::Trigger *t1, const cafmaker::Trigger *t2)
  {
    return cmp(*t1, *t2);
  }
};

// -------------------------------------------------
// return type: each element of outer vector corresponds to one group of matched triggers
std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>>
buildTriggerList(std::map<const cafmaker::IRecoBranchFiller*, std::deque<cafmaker::Trigger>> triggersByFiller,
                 unsigned int trigMatchMaxDT)
{
  // I don't want to keep typing `cafmaker::LOG_S("buildTriggerList()")` every time,
  // and the preamble to the logger resets after the first use
  auto LOG = [&]() -> const cafmaker::Logger & { return cafmaker::LOG_S("buildTriggerList()"); };

  std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>> ret;

  // don't assume input comes in sorted
  LOG().DEBUG() << "Considering triggers from the following reco branch fillers: \n";
  for (auto & fillerTrigPair : triggersByFiller)
  {
    std::sort(fillerTrigPair.second.begin(), fillerTrigPair.second.end(), triggerTimeCmp());
    LOG().DEBUG() << "   " << fillerTrigPair.first->GetName() << "\n";
  }

  while (!triggersByFiller.empty())
  {
    // look at the first element of each reco filler stream.
    LOG().VERBOSE() << "   Considering the earliest triggers in each stream:\n";
    std::vector<const cafmaker::Trigger*> firstTrigs;
    for (const auto &it: triggersByFiller)
    {
      LOG().VERBOSE() << "     " << it.first->GetName() << " --> (id = " << it.second[0].evtID
                    << ", time = " << (it.second[0].triggerTime_s + it.second[0].triggerTime_ns/1e9) << " s)\n";
      firstTrigs.push_back(&it.second[0]);
    }

    // the earliest one will be our next group seed.
    auto groupSeedIt = std::min_element(firstTrigs.begin(), firstTrigs.end(), triggerTimePtrCmp());
    LOG().VERBOSE() << "    Next trigger group seed: (" << (*groupSeedIt)->evtID << ", "
                  << (*groupSeedIt)->triggerTime_s +(*groupSeedIt)->triggerTime_ns/1e9
                  << ")\n";

    // pull that one out of its original container so we don't reconsider it in the next loop iteration
    auto seedFillerIt = triggersByFiller.begin();
    std::advance(seedFillerIt, std::distance(firstTrigs.begin(), groupSeedIt));
    ret.push_back({{seedFillerIt->first, std::move(seedFillerIt->second.front())}});  // note that we're stealing the contents of the element from its deque since we're about to pop it anyway
    seedFillerIt->second.pop_front();

    // now consider the other reco filler streams.
    // do they have any events in them that should go in this group?
    std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>> & trigGroup = ret.back();
    const cafmaker::Trigger & trigSeed = trigGroup.front().second;
    for (auto & fillerTrigPair : triggersByFiller)
    {
      // we don't want to consider the stream we're already working with.
      // (but don't continue, because we want to remove this stream from the
      //  map if it's empty, per below)
      if (fillerTrigPair.first != seedFillerIt->first)
      {
        // we will only take at most one trigger from each of the other streams.
        // since the seed was the earliest one out of all the triggers,
        // we only need to check the first one in each other stream
        if (doTriggersMatch( trigSeed, fillerTrigPair.second.front(), trigMatchMaxDT))
        {
            trigGroup.push_back({fillerTrigPair.first, std::move(fillerTrigPair.second.front())});
            fillerTrigPair.second.pop_front();
        }
      }

      // if there are no more elements in this reco filler stream,
      // remove them from consideration
      if (fillerTrigPair.second.empty())
        triggersByFiller.erase(fillerTrigPair.first);

    }
  }

  return ret;
}

// -------------------------------------------------
// main loop function
void loop(CAF &caf,
          cafmaker::Params &par,
          const std::vector<std::string> & ghepFilenames,
          const std::vector<std::unique_ptr<cafmaker::IRecoBranchFiller>> &recoFillers)
{
  // if this is a data file, there won't be any truth, of course,
  // but the TruthMatching knows not to try to do anything with a null gtree
  cafmaker::Logger::THRESHOLD thresh = cafmaker::Logger::parseStringThresh(par().cafmaker().verbosity());
  cafmaker::TruthMatcher truthMatcher(ghepFilenames, caf.mcrec,
                                      [&caf](const genie::NtpMCEventRecord* mcrec){ return caf.StoreGENIEEvent(mcrec); });
  truthMatcher.SetLogThrehsold(thresh);
  // figure out which triggers we need to loop over between the various reco fillers
  std::map<const cafmaker::IRecoBranchFiller*, std::deque<cafmaker::Trigger>> triggersByRBF;
  for (const std::unique_ptr<cafmaker::IRecoBranchFiller>& filler : recoFillers)
    triggersByRBF.insert({filler.get(), filler->GetTriggers()});
<<<<<<< HEAD
  std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>>
    groupedTriggers = buildTriggerList(triggersByRBF, par().cafmaker().trigMatchDT());

=======
  
  std::vector<std::vector<std::pair<const cafmaker::IRecoBranchFiller*, cafmaker::Trigger>>> groupedTriggers = buildTriggerList(triggersByRBF);
>>>>>>> fixed last differences with MASTER
  // sanity checks
  if (par().cafmaker().first() > static_cast<int>(groupedTriggers.size()))
  {
    std::cerr << "Requested starting event (" << par().cafmaker().first() << ") "
              << "is larger than total number of triggers (" << groupedTriggers.size() << ".\n"
              << "Do nothing ...\n";
    return;
  }
  int start = par().cafmaker().first();
  int N = par().cafmaker().numevts() > 0 ? par().cafmaker().numevts() : static_cast<int>(groupedTriggers.size()) - par().cafmaker().first();
  if (N < 1)
  {
    std::cerr << "Requested number of events (" << N << ") is non-positive!  Abort.\n";
    abort();
  }

  // Main event loop
  cafmaker::Progress progBar("Processing " + std::to_string(N - start) + " triggers");
  for( int ii = start; ii < start + N; ++ii )
  {
    // don't bother with updating the prog bar if we're going to be spamming lots of messages
    if (thresh >= cafmaker::Logger::THRESHOLD::WARNING)
      progBar.SetProgress( static_cast<double>(ii - start)/N );
    else
      cafmaker::LOG_S("loop()").INFO() << "Processing trigger: " << ii << "\n";

    // reset (the default constructor initializes its variables)
    caf.setToBS();


    // hand off to the correct reco filler(s).
    for (const auto & fillerTrigPair : groupedTriggers[ii])
    {
      cafmaker::LOG_S("loop()").INFO() << "Global trigger idx : " << ii << ", reco filler: '" << fillerTrigPair.first->GetName() << "', reco trigger eventID: " << fillerTrigPair.second.evtID << "\n";
      fillerTrigPair.first->FillRecoBranches(fillerTrigPair.second, caf.sr, par, &truthMatcher);
    }

    caf.fill();
  }
  progBar.Done();

  // set other metadata
  caf.meta_run = par().runInfo().run();
  caf.meta_subrun = par().runInfo().subrun();

}

// -------------------------------------------------

int main( int argc, char const *argv[] )
{


  progopt::variables_map vars = parseCmdLine(argc, argv);

  cafmaker::Params par = parseConfig(vars["fcl"].as<std::string>(), vars);

  cafmaker::Logger::THRESHOLD logThresh = cafmaker::Logger::parseStringThresh(par().cafmaker().verbosity());
  cafmaker::LOG_S().SetThreshold(logThresh);

  std::vector<std::string> GHEPFiles;
  par().cafmaker().GHEPFiles(GHEPFiles);  // fills the vector in if the key is found

  CAF caf(par().cafmaker().outputFile(), par().cafmaker().nusystsFcl(), par().cafmaker().makeFlatCAF(), !GHEPFiles.empty());

  loop(caf, par, GHEPFiles, getRecoFillers(par, logThresh));

  caf.version = 5;
  printf( "Run %d POT %g\n", caf.meta_run, caf.pot );
  caf.fillPOT();

  std::cout << "Writing CAF" << std::endl;
  caf.write();

  return 0;
}
