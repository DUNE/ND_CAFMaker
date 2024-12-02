/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi
/// \date    Apr. 2022
///

#include <iostream>

#include "TBranch.h"
#include "TLeaf.h"
#include "SANDRecoBranchFiller.h"
#include "truth/FillTruth.h"

// #ifdef ENABLE_SAND
#warning Including SANDRecoBranchFiller in build

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

// #include "SANDReco/SANDRecord.h"

namespace cafmaker
{

  SANDRecoBranchFiller::SANDRecoBranchFiller(const std::string &SANDRecoFilename)
      : IRecoBranchFiller("SAND"),
        fTriggers(),
        fLastTriggerReqd(fTriggers.end())
  {
    fSANDRecoFile = new TFile(SANDRecoFilename.c_str());
    // NDSANDRecoTree = (TTree*) fSANDRecoFile->Get("tEvent");
    NDSANDRecoTree = (TTree *)fSANDRecoFile->Get("tReco");
    
    // fEvent = new event;

    if (!fSANDRecoFile->IsZombie())
      SetConfigured(true);
  }

  // void SANDRecoBranchFiller::_FillRecoBranches(std::size_t N, std::size_t ii,
  //                                              caf::StandardRecord &sr,
  //                                              const cafmaker::Params &par) const
  void SANDRecoBranchFiller::_FillRecoBranches(const Trigger &trigger,
                                               caf::StandardRecord &sr,
                                               const cafmaker::Params &par, const TruthMatcher *truthMatcher) const
  {


    auto it_start = (fLastTriggerReqd == fTriggers.end()) ? fTriggers.cbegin() : fLastTriggerReqd;
    auto itTrig = std::find(it_start, fTriggers.cend(), trigger);
    if (itTrig == fTriggers.end())
    {
      LOG.FATAL() << "Reco branch filler '" << GetName() << "' could not find trigger with evtID == " << trigger.evtID << "!  Abort.\n";
      abort();
    }
    std::size_t idx = std::distance(fTriggers.cbegin(), itTrig);
    int event_num = idx;
    sr.meta.sand.event = event_num;
    LOG.VERBOSE() << "    Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID << ", internal evt idx = " << idx << ".\n";
    FillECalClusters(truthMatcher, sr, event_num);
    
    

  }

  void SANDRecoBranchFiller::FillECalClusters(const TruthMatcher * truthMatch,
                                                caf::StandardRecord &sr, int event_num) const
  {
  //NDSANDRecoTree->GetEntry(idx);

      // TBranch *cluster = NDSANDRecoTree->GetBranch("cluster");
      // if (!cluster)
      // {
      //   std::cerr << "Error: cluster branch not found in tReco tree" << std::endl;
      //   return;
      // }

      //TLeaf *cluster_energy = cluster->FindLeaf("e");
      NDSANDRecoTree->GetEntry(event_num);
      TLeaf* cluster_energy = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.e");
      if (!cluster_energy)
      {
        std::cerr << "Error: energy not found in cluster" << std::endl;
        return;
      }

      int num_clusters = cluster_energy->GetLen();
      size_t nclusters = num_clusters;
      std::cout << "EVENT: " << event_num << ",num clusters: "<< num_clusters << std::endl; 
      sr.nd.sand.ixn.resize(1);
      //std::cout << "Event " << trigger.evtID << ", cluster energy:" << std::endl; 
      sr.nd.sand.ixn[0].nshowers = nclusters;
      sr.nd.sand.ixn[0].showers.resize(num_clusters);
      std::vector<double> SANDECALClusterEnergy;
      for (int i = 0; i < num_clusters; i++)
      {
        double energy = cluster_energy->GetValue(i);
        SANDECALClusterEnergy.push_back(energy);
        //std::cout << "  Cluster " << i << ": energy(from fECALCLusternergy) = " << SANDECALClusterEnergy.at(i) << " MeV" << std::endl;
        //shower.Evis = energy; 
        std::cout << "  Cluster " << i << ": energy = " << energy << " MeV" << std::endl;
        sr.nd.sand.ixn[0].showers[i].Evis = energy;
      }
  }

  // todo: this is a placeholder
  std::deque<Trigger> SANDRecoBranchFiller::GetTriggers(int triggerType) const
  {
    std::deque<Trigger> triggers;
    size_t n_entries = NDSANDRecoTree->GetEntries();
    if (fTriggers.empty())
    {
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "\n";
      fTriggers.reserve(n_entries);

      for (size_t i = 0; i < n_entries; i++)
      {

        const int placeholderTriggerType = 0;
        // fixme: this check needs to be fixed when we have trigger type info
        if (triggerType >= 0 && triggerType != placeholderTriggerType)
        {
          LOG.VERBOSE() << "    skipping this event" << "\n";
          continue;
        }

        NDSANDRecoTree->GetEntry(i);

        fTriggers.emplace_back();
        Trigger &trig = fTriggers.back();
        trig.evtID = i;

        // todo: these are placeholder values until we can propagate enough info through the reco files
        trig.triggerType = 0;
        trig.triggerTime_s = i;
        trig.triggerTime_ns = 0.;

        triggers.push_back(trig);

        LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
                      << ", triggerType=" << trig.triggerType
                      << ", triggerTime_s=" << trig.triggerTime_s
                      << ", triggerTime_ns=" << trig.triggerTime_ns
                      << "\n";
      }
      fLastTriggerReqd = fTriggers.end(); // since we just modified the list, any iterators have been invalidated
    }
    else
    {
      for (const Trigger &trigger : fTriggers)
      {
        if (triggerType < 0 || triggerType == fTriggers.back().triggerType)
          triggers.push_back(trigger);
      }
    }
    return triggers;
  }
}

// #else // ENABLE_SAND

#warning Not configured to build SANDRecoBranchFiller. Must set SANDRECO_INC and SANDRECO_LIB environment variables

namespace
{
  void error_msg()
  {
    std::cerr << "\n\nSAND Reco support was not enabled in your build. \n"
              << " Either avoid setting `nd_cafmaker.CAFMakerSettings.SANDRecoFile` in your FCL\n"
              << " or set $SANDRECO_INC and $SANDRECO_LIB in your environment and do a clean rebuild of ND_CAFMaker...\n";
  }
}

// namespace cafmaker
// {
//   SANDRecoBranchFiller::SANDRecoBranchFiller(const std::string &)
//       : IRecoBranchFiller("SAND")
//   {
//     error_msg();
//     abort();
//   }

//   void SANDRecoBranchFiller::
//       _FillRecoBranches(const Trigger &, caf::StandardRecord &, const cafmaker::Params &,
//                         const TruthMatcher *truthMatcher) const
//   {
//     error_msg();
//     abort();
//   }

//   // todo: this is a placeholder
//   std::deque<Trigger> SANDRecoBranchFiller::GetTriggers(int triggerType) const
//   {
//     return std::deque<Trigger>();
//   }
// }

// #endif // ENABLE_SAND
