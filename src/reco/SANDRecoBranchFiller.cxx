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
    NDSANDRecoTree = (TTree *)fSANDRecoFile->Get("tReco");
//    NDSANDEventTree = (TTree *)fSANDRecoFile->Get("tEvent");
//    NDSANDRecoTree->AddFriend(NDSANDEventTree);
     
    // fEvent = new event;
    std::cout<<"Number of entries "<<NDSANDRecoTree->GetEntries()<<std::endl;
    

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
    NDSANDRecoTree->GetEntry(event_num);
    
    sr.meta.sand.event = event_num;
    LOG.VERBOSE() << "    Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID << ", internal evt idx = " << idx << ".\n";
    FillInteractions(truthMatcher, sr); 
    FillECalClusters(truthMatcher, sr);

  }

void SANDRecoBranchFiller::FillInteractions(const TruthMatcher * truthMatch,
                                               caf::StandardRecord &sr) const
  {
    // F says: our GArSoft samples contain one neutrino interaction
    //         per event/trigger (I think we can easily change that?)
    //         so we only need one interaction object in the common
    //         reco branch so far
    sr.common.ixn.sandreco.reserve(1);
    sr.common.ixn.nsandreco = 1;

    //         now we have one interaction per event/trigger
    //         so we simply need to create an SRInteraction and assign a
    //         dummy id of 1 to it

    caf::SRInteraction interaction;
    interaction.id  = 1;
   /* 
    TLeaf* Enureco = (TLeaf*) NDSANDEventTree->GetListOfLeaves()->FindObject("Enureco");
    TLeaf* pxnureco = (TLeaf*) NDSANDEventTree->GetListOfLeaves()->FindObject("pxnureco");
    TLeaf* pynureco = (TLeaf*) NDSANDEventTree->GetListOfLeaves()->FindObject("pynureco");
    TLeaf* pznureco = (TLeaf*) NDSANDEventTree->GetListOfLeaves()->FindObject("pznureco");
    
    TLeaf* xvertex = (TLeaf*) NDSANDEventTree->GetListOfLeaves()->FindObject("x");
    TLeaf* yvertex = (TLeaf*) NDSANDEventTree->GetListOfLeaves()->FindObject("y");
    TLeaf* zvertex = (TLeaf*) NDSANDEventTree->GetListOfLeaves()->FindObject("z");

    if (!Enureco)
      {
        std::cerr << "Error: energy not found in cluster" << std::endl;
        return;
      }
    */
    
    //interaction.vtx  = caf::SRVector3D(xvertex->GetValue(0), yvertex->GetValue(0), zvertex->GetValue(0)); 
    interaction.vtx  = caf::SRVector3D(10.0, 11.0, 12.0); 

 //   interaction.Enu.calo = Enureco->GetValue(0);
 //   interaction.Enu.momentum.SetXYZ(pxnureco->GetValue(0),pynureco->GetValue(0),pznureco->GetValue(0));
    float Enureco=12.0;
    interaction.Enu.calo = Enureco;
    interaction.Enu.momentum.SetXYZ(20.0,21.0, 22.0);
 
    int nparticles=5; //legge la size del vettore particle del tree
 

    interaction.part.npida = nparticles;  
    interaction.part.pida.resize(nparticles);  

    for(int i=0; i<nparticles; i++){
     interaction.part.pida[i].primary=false;
     interaction.part.pida[i].pdg=10.0*i;
     interaction.part.pida[i].start=caf::SRVector3D(i,i,i);
     interaction.part.pida[i].end=caf::SRVector3D(i*10,i*10,i*10);
     interaction.part.pida[i].p=caf::SRVector3D(i*0.1,i*0.1,i*0.1);
    }


    sr.common.ixn.sandreco.push_back(std::move(interaction));
}


  void SANDRecoBranchFiller::FillECalClusters(const TruthMatcher * truthMatch,
                                                caf::StandardRecord &sr) const
  {

      // TBranch *cluster = NDSANDRecoTree->GetBranch("cluster");
      // if (!cluster)
      // {
      //   std::cerr << "Error: cluster branch not found in tReco tree" << std::endl;
      //   return;
      // }

      //TLeaf *cluster_energy = cluster->FindLeaf("e");
      TLeaf* cluster_energy = (TLeaf*) NDSANDRecoTree->GetListOfLeaves()->FindObject("cluster.e");
      if (!cluster_energy)
      {
        std::cerr << "Error: energy not found in cluster" << std::endl;
        return;
      }

      int num_clusters = cluster_energy->GetLen();
      size_t nclusters = num_clusters;
      std::cout <<" num clusters: "<< num_clusters << std::endl; 
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

      std::cout<<"creating triggers"<<std::endl;

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
    std::cout<<"Trigger completed!"<<std::endl;
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
