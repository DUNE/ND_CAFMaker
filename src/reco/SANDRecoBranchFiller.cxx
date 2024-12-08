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

      //TLeaf *energy = cluster->FindLeaf("e");
      NDSANDRecoTree->GetEntry(event_num);
      TLeaf* energy = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.e");
      TLeaf* position_x = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.x");
      TLeaf* position_y = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.y");
      TLeaf* position_z = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.z");
      TLeaf* var_positionx = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.varx");
      TLeaf* var_positiony = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.vary");
      TLeaf* var_positionz = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.varz");
      TLeaf* cluster_time = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.t");
      TLeaf* apex_x = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.ax");
      TLeaf* apex_y = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.ay");
      TLeaf* apex_z = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.az");
      TLeaf* dir_x = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.sx");
      TLeaf* dir_y = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.sy");
      TLeaf* dir_z = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.sz");
      TLeaf* id = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.tid");

      // if (!energy || !position_x || !position_y  || !position_z || !var_positionx || 
      // !var_positiony || !var_positionz || !cluster_time|| !apex_x || !apex_y || !apex_z ||
      //  !dir_x || !dir_y || !dir_z )//!id)
      // {
      //   std::cerr << "Error: leaf not found in cluster" << std::endl;
      //   return;
      // }

       if (!energy )
      {
        std::cerr << "Error: energy leaf not found in cluster" << std::endl;
        return;
      }

      if (!position_x || !position_y || !position_z)
      {
          std::cerr << "Error: position leaf not found in cluster" << std::endl;
          return;
      }
      
      if (!var_positionx || !var_positiony || !var_positionz)
      {
          std::cerr << "Error: var position leaf not found in cluster" << std::endl;
          return;
      }

      if (!cluster_time)
      {
          std::cerr << "Error: time leaf not found in cluster" << std::endl;
          return;
      }

      if (!apex_x || !apex_y || !apex_z)
      {
          std::cerr << "Error: apex position leaf not found in cluster" << std::endl;
          return;
      }

      if (!dir_x || !dir_y || !dir_z)
      {
          std::cerr << "Error: direction leaf not found in cluster" << std::endl;
          return;
      }
      
      if (!id)
      {
          std::cerr << "Error: id leaf not found in cluster" << std::endl;
          return;
      }

      size_t n_clusters = energy->GetLen(); //better to use a cluster size 
      std::cout << "EVENT: " << event_num << ", num clusters: "<< n_clusters << std::endl; 
      
      sr.nd.sand.ixn.resize(1);
      sr.nd.sand.ixn[0].nclusters = n_clusters;
      sr.nd.sand.ixn[0].ECALClusters.resize(n_clusters);
      
      for (int i = 0; i < n_clusters; i++)
      { 
        std::cout << "  Cluster " << i << ": energy = " << energy->GetValue(i) << " MeV" << std::endl;
        std::cout << "position x: " << position_x->GetValue(i) << ",position y: " << position_y->GetValue(i) << ", position z:" << position_z->GetValue(i) << std::endl; 
        std::cout << "var x: " << var_positionx->GetValue(i) << ", var y: " << var_positiony->GetValue(i) << ", var z:" << var_positionz->GetValue(i) << std::endl; 
        std::cout << "time: "<< cluster_time->GetValue(i) << std::endl; 
        std::cout << "apex x: " << apex_x->GetValue(i) << ", apex y: " << apex_y->GetValue(i) << ", apex z:" << apex_z->GetValue(i) << std::endl; 
        std::cout << "dir x: " << dir_x->GetValue(i) << ", dir y: " << dir_y->GetValue(i) << ", dir z:" << dir_z->GetValue(i) << std::endl; 
        std::cout << "id: "<< id->GetValue(i) << std::endl; 

        sr.nd.sand.ixn[0].ECALClusters[i].E = energy->GetValue(i);
        sr.nd.sand.ixn[0].ECALClusters[i].position.SetXYZ(position_x->GetValue(i), position_y->GetValue(i), position_z->GetValue(i));
        sr.nd.sand.ixn[0].ECALClusters[i].var_position.SetXYZ(var_positionx->GetValue(i), var_positiony->GetValue(i), var_positionz->GetValue(i));
        sr.nd.sand.ixn[0].ECALClusters[i].time = cluster_time->GetValue(i);
        sr.nd.sand.ixn[0].ECALClusters[i].start.SetXYZ(apex_x->GetValue(i), apex_y->GetValue(i), apex_z->GetValue(i));
        sr.nd.sand.ixn[0].ECALClusters[i].direction.SetXYZ(dir_x->GetValue(i), dir_y->GetValue(i), dir_z->GetValue(i));
        sr.nd.sand.ixn[0].ECALClusters[i].id = id->GetValue(i);
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
