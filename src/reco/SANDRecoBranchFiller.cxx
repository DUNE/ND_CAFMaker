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

#ifdef ENABLE_SAND
#warning Including SANDRecoBranchFiller in build

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "struct.h"

namespace cafmaker
{

  SANDRecoBranchFiller::SANDRecoBranchFiller(const std::string &SANDRecoFilename)
      : IRecoBranchFiller("SAND"),
        fTriggers(),
        fLastTriggerReqd(fTriggers.end())
  {
  
    if (SANDRecoFilename.empty()) {
      std::cerr << "ERROR: SANDRecoBranchFiller: SANDRecoFilename is empty!" << std::endl;
      SetConfigured(false);
      return;
    }
    fSANDRecoFile = new TFile(SANDRecoFilename.c_str());
    if (!fSANDRecoFile || fSANDRecoFile->IsZombie()) {
      std::cerr << "ERROR: SANDRecoBranchFiller: Could not open file " << SANDRecoFilename << std::endl;
      SetConfigured(false);
      return;
    }
    NDSANDRecoTree = (TTree *)fSANDRecoFile->Get("tReco");
    if (!NDSANDRecoTree)
    {
      std::cerr << "Error: NDSANDRecoTree is null. Tree 'tReco' not found in file." << std::endl;
      fSANDRecoFile->ls(); 
      SetConfigured(false);
      return;
    }
   
    SetConfigured(true);
  }

  
  void SANDRecoBranchFiller::_FillRecoBranches(const Trigger &trigger,
                                               caf::StandardRecord &sr,
                                               const cafmaker::Params &par, const TruthMatcher *truthMatcher) const
  {
  
    if (!NDSANDRecoTree) {
      std::cerr << "ERROR: SANDRecoBranchFiller: NDSANDRecoTree is null!" << std::endl;
      return;
    }
   
    auto it_start = (fLastTriggerReqd == fTriggers.end()) ? fTriggers.cbegin() : fLastTriggerReqd;
    auto itTrig = std::find(it_start, fTriggers.cend(), trigger);
    if (itTrig == fTriggers.end())
    {
      LOG.FATAL() << "Reco branch filler '" << GetName() << "' could not find trigger with evtID == " << trigger.evtID << "!  Abort.\n";
      abort();
    }
    
    std::size_t idx = std::distance(fTriggers.cbegin(), itTrig);
    int event_num = idx;
    NDSANDRecoTree->GetEntry(event_num);
    NDSANDEventTree->GetEntry(event_num);
    
    sr.meta.sand.event = event_num;
    LOG.VERBOSE() << "    Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID << ", internal evt idx = " << idx << ".\n";
    std::vector<cluster> cl;
    std::vector<track> tr;
    auto pcl = &cl;
    auto ptr = &tr;
    NDSANDRecoTree->SetBranchAddress("cluster", &pcl);
    NDSANDRecoTree->SetBranchAddress("track", &ptr);
    NDSANDRecoTree->GetEntry(event_num);
    FillECalClusters(truthMatcher, sr, cl);
   
    FillTracks(truthMatcher, sr, tr);
  }

  void SANDRecoBranchFiller::FillECalClusters(const TruthMatcher *truthMatch,
                                              caf::StandardRecord &sr, std::vector<cluster> &cl) const
  {
    
    size_t n_clusters = cl.size();
    for (const auto &c : cl)
    {
      std::cout << "number of clusters: " << n_clusters << std::endl;
    }

    sr.nd.sand.ixn.resize(1);
    sr.nd.sand.ixn[0].nclusters = n_clusters;
    sr.nd.sand.ixn[0].ECALClusters.resize(n_clusters);

    for (auto i = 0; i < n_clusters; i++)
    {

      sr.nd.sand.ixn[0].ECALClusters[i].E = cl[i].e;
      sr.nd.sand.ixn[0].ECALClusters[i].position.SetXYZ(cl[i].x, cl[i].y, cl[i].z);
      sr.nd.sand.ixn[0].ECALClusters[i].var_position.SetXYZ(cl[i].varx, cl[i].vary, cl[i].varz);
      sr.nd.sand.ixn[0].ECALClusters[i].time = cl[i].t; 
      sr.nd.sand.ixn[0].ECALClusters[i].start.SetXYZ(cl[i].ax, cl[i].ay, cl[i].az);
      sr.nd.sand.ixn[0].ECALClusters[i].direction.SetXYZ(cl[i].sx, cl[i].sy, cl[i].sz);
      sr.nd.sand.ixn[0].ECALClusters[i].num_cells = cl[i].reco_cells.size();
      sr.nd.sand.ixn[0].ECALClusters[i].id = cl[i].tid;
    }
    
    cl.clear();
  }

  void SANDRecoBranchFiller::FillTracks(const TruthMatcher *truthMatch,
                                        caf::StandardRecord &sr, std::vector<track> &tr) const
  {

    size_t n_tracks = tr.size();
     for (const auto &t : tr)
     {
       std::cout << "number of tracks: " <<  n_tracks << std::endl;
     }


    sr.nd.sand.ixn.resize(1);
    sr.nd.sand.ixn[0].ntracks = n_tracks;
    sr.nd.sand.ixn[0].tracks.resize(n_tracks);

    for (int i = 0; i < n_tracks; i++)
    {
      sr.nd.sand.ixn[0].tracks[i].start.SetXYZ(tr[i].x0, tr[i].y0, tr[i].z0);
      sr.nd.sand.ixn[0].tracks[i].qual = tr[i].chi2_cr;
    }
  }

  // todo: this is a placeholder
  std::deque<Trigger> SANDRecoBranchFiller::GetTriggers(int triggerType) const
  {
    std::deque<Trigger> triggers;
    size_t n_entries = NDSANDRecoTree->GetEntries();
    std::cout << "n entries: " << n_entries << std::endl;

    if (fTriggers.empty())
    {
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "\n";
      fTriggers.reserve(n_entries);

      std::cout<<"creating triggers"<<std::endl;

      for (size_t i = 0; i < n_entries; i++)
      {
        std::cout << "entry n: " << i << std::endl;
        const int placeholderTriggerType = 0;
        // fixme: this check needs to be fixed when we have trigger type info
        
        if (triggerType >= 0 && triggerType != placeholderTriggerType)
        {
          LOG.VERBOSE() << "    skipping this event" << "\n";
          std::cout << "skipping this event" << std::endl;
          continue;
        }

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
        {
          triggers.push_back(trigger);
        }
      }
    }
    std::cout<<"Trigger completed!"<<std::endl;
    return triggers;
  }
}

#else
#include "reco/SANDRecoBranchFiller.h"
#include <deque>
#include <iostream>
namespace cafmaker {
SANDRecoBranchFiller::SANDRecoBranchFiller(const std::string &)
    : IRecoBranchFiller("SAND")
{
    std::cerr << "[WARNING] SANDRecoBranchFiller: SAND code is not enabled in this build.\n"
                 "To enable SAND support, build with 'make ENABLE_SAND_DICT=1'.\n";
}
std::deque<Trigger> SANDRecoBranchFiller::GetTriggers(int) const { return {}; }
void SANDRecoBranchFiller::_FillRecoBranches(const Trigger &, caf::StandardRecord &, const cafmaker::Params &, const TruthMatcher *) const {}
}
#endif
