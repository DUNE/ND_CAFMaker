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

    //Now we want to read also tEvent with GRAIN info 

    NDSANDEventTree = (TTree *)fSANDRecoFile->Get("tEvent");
    if (!NDSANDEventTree)
    {
      std::cerr << "Error: NDSANDEventTree is null. Tree 'tEvent' not found in file." << std::endl;
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

    if (!NDSANDEventTree) {
      std::cerr << "ERROR: SANDRecoBranchFiller: NDSANDEventTree is null!" << std::endl;
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
    
    sr.meta.sand.event = event_num;
    
    LOG.VERBOSE() << "    Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID << ", internal evt idx = " << idx << ".\n";
    std::vector<cluster> cl;
    std::vector<track> tr;
    std::vector<track_grain_lens> trackgrain;
    auto pcl = &cl;
    auto ptr = &tr;
    auto ptrgrain = &trackgrain;

    NDSANDRecoTree->SetBranchAddress("cluster", &pcl);
    NDSANDRecoTree->SetBranchAddress("track", &ptr);
    NDSANDRecoTree->SetBranchAddress("track_grain_lens", &ptrgrain);
    NDSANDRecoTree->GetEntry(event_num);
    event* evt = nullptr;
    NDSANDEventTree->SetBranchAddress("event", &evt);
    FillECalClusters(truthMatcher, sr, cl);
   
    /*FillTracks(truthMatcher, sr, tr);

    FillGRAINTracks(truthMatcher, sr, trackgrain);*/

    FillSANDTracks(truthMatcher, sr, tr, trackgrain);
    //new -- neutrino level info 

    NDSANDEventTree->GetEntry(event_num);  // legge l'evento corrispondente al trigger
    FillRecoNu(sr, evt, event_num);  
    FillRecoParticle(sr, evt, event_num);  
    FillTrueNu(sr, evt, event_num);
    FillTruePrimaryParticles(sr, evt, event_num);
    FillTrueSecondaryParticles(sr, evt, event_num);

  }

  //NEW FillEcalCLusters -- changes of SRSAND.H
  void SANDRecoBranchFiller::FillECalClusters(const TruthMatcher *truthMatch,
                                            caf::StandardRecord &sr,
                                            std::vector<cluster> &cl) const
{
    size_t n_clusters = cl.size();

    sr.nd.sand.ixn.resize(1);

    sr.nd.sand.ixn[0].ecal.nclusters = n_clusters;
    sr.nd.sand.ixn[0].ecal.clusters.resize(n_clusters);

    for (size_t i = 0; i < n_clusters; i++)
    {
        sr.nd.sand.ixn[0].ecal.clusters[i].E = cl[i].e;
        sr.nd.sand.ixn[0].ecal.clusters[i].position.SetXYZ(cl[i].x, cl[i].y, cl[i].z);
        sr.nd.sand.ixn[0].ecal.clusters[i].var_position.SetXYZ(cl[i].varx,
                                                               cl[i].vary,
                                                               cl[i].varz);
        sr.nd.sand.ixn[0].ecal.clusters[i].time = cl[i].t;
        sr.nd.sand.ixn[0].ecal.clusters[i].start.SetXYZ(cl[i].ax, cl[i].ay, cl[i].az);
        sr.nd.sand.ixn[0].ecal.clusters[i].direction.SetXYZ(cl[i].sx, cl[i].sy, cl[i].sz);
        sr.nd.sand.ixn[0].ecal.clusters[i].num_cells = cl[i].reco_cells.size();
        sr.nd.sand.ixn[0].ecal.clusters[i].id = cl[i].tid;
    }

    cl.clear();
}

  /*void SANDRecoBranchFiller::FillECalClusters(const TruthMatcher *truthMatch,
                                              caf::StandardRecord &sr, std::vector<cluster> &cl) const
  {
   
    size_t n_clusters = cl.size();

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

    sr.nd.sand.ixn.resize(1);
    sr.nd.sand.ixn[0].ntracks = n_tracks;
    sr.nd.sand.ixn[0].tracks.resize(n_tracks);

    for (int i = 0; i < n_tracks; i++)
    {
      sr.nd.sand.ixn[0].tracks[i].start.SetXYZ(tr[i].x0, tr[i].y0, tr[i].z0);
      sr.nd.sand.ixn[0].tracks[i].qual = tr[i].chi2_cr;
    }
    tr.clear();
  }

  void SANDRecoBranchFiller::FillGRAINTracks(const TruthMatcher *truthMatch,
                                        caf::StandardRecord &sr, std::vector<track_grain_lens> &trackgrain) const
  {

    size_t n_tracks = trackgrain.size();

    //Version with GRAIN tracks 

    size_t n_GRAINtracks= trackgrain.size();
    sr.nd.sand.ixn.resize(1);
    sr.nd.sand.ixn[0].nGRAINtracks = n_GRAINtracks;
    sr.nd.sand.ixn[0].GRAINtracks.resize(n_GRAINtracks);

    sr.nd.sand.ixn.resize(1);
    sr.nd.sand.ixn[0].ntracks = n_tracks;
    sr.nd.sand.ixn[0].tracks.resize(n_tracks);

    for (int i = 0; i < n_tracks; i++)
    {
      sr.nd.sand.ixn[0].tracks[i].start.SetXYZ(trackgrain[i].xstart, trackgrain[i].ystart, trackgrain[i].zstart);
      sr.nd.sand.ixn[0].tracks[i].end.SetXYZ(trackgrain[i].xend, trackgrain[i].yend, trackgrain[i].zend);
    }

    //Version with GRAIN tracks 
    for (int i = 0; i < n_GRAINtracks; i++)
    {
      sr.nd.sand.ixn[0].GRAINtracks[i].start.SetXYZ(trackgrain[i].xstart, trackgrain[i].ystart, trackgrain[i].zstart);
      sr.nd.sand.ixn[0].GRAINtracks[i].end.SetXYZ(trackgrain[i].xend, trackgrain[i].yend, trackgrain[i].zend);
    }

    trackgrain.clear();
  }*/

  //new -- FillSANDTracks 
  //Jeremy's suggestion

  void SANDRecoBranchFiller::FillSANDTracks(const TruthMatcher *truthMatch,
                                          caf::StandardRecord &sr,
                                          std::vector<track> &tr,
                                          std::vector<track_grain_lens> &trackgrain) const
{
    size_t n_GRAINtracks = trackgrain.size();

    sr.nd.sand.ixn.resize(1);

    // Riempio i tracks di GRAIN
    sr.nd.sand.ixn[0].grain.ntracks = n_GRAINtracks;
    sr.nd.sand.ixn[0].grain.tracks.resize(n_GRAINtracks);

    for (size_t i = 0; i < n_GRAINtracks; i++)
    {
        sr.nd.sand.ixn[0].grain.tracks[i].start.SetXYZ(trackgrain[i].xstart,
                                                       trackgrain[i].ystart,
                                                       trackgrain[i].zstart);
                                                       
        sr.nd.sand.ixn[0].grain.tracks[i].end.SetXYZ(trackgrain[i].xend,
                                                     trackgrain[i].yend,
                                                     trackgrain[i].zend);

        sr.nd.sand.ixn[0].grain.tracks[i].dir.SetXYZ(trackgrain[i].M1dirX,
                                                     trackgrain[i].M1dirY,
                                                     trackgrain[i].M1dirZ);                                             
                                                    
    }

    size_t n_tracks = tr.size();

    // Riempio i tracks di Tracker
    sr.nd.sand.ixn[0].tracker.ntracks = n_tracks;
    sr.nd.sand.ixn[0].tracker.tracks.resize(n_tracks);

    for (size_t i = 0; i < n_tracks; i++)
    {
        sr.nd.sand.ixn[0].tracker.tracks[i].start.SetXYZ(tr[i].x0,
                                                         tr[i].y0,
                                                         tr[i].z0);
        sr.nd.sand.ixn[0].tracker.tracks[i].qual = tr[i].chi2_cr;
    }

    tr.clear();
    trackgrain.clear();
}

  //new -- neutrino level info 
  
  void SANDRecoBranchFiller::FillRecoNu(caf::StandardRecord &sr, const event* evt, int event_num) const {

    sr.common.ixn.sandreco.resize(1);
    
    sr.common.ixn.sandreco[0].vtx.SetXYZ(evt->x, evt->y, evt->z);
    std::cout << "=== DEBUG FillRecoNu ===" << std::endl;
    std::cout << "Evento: posizione (" << evt->x << ", " << evt->y << ", " << evt->z << ") " << std::endl;
    std::cout << "Evento: Enureco (" << evt->Enureco << ") " << std::endl;
    std::cout << "Evento: Truth (" << event_num << ") " << std::endl;
    std::cout << "DEBUG: riempiendo Enu.part_energy_sum con Enureco = " << evt->Enureco << std::endl;
    sr.common.ixn.sandreco[0].Enu.part_energy_sum = evt->Enureco;

  }

  void SANDRecoBranchFiller::FillTrueNu(caf::StandardRecord &sr,
                                      const event* evt,
                                      int event_num) const
{
    // =========================
    // Reset branch truth neutrino
    // =========================
    //sr.mc.nu.clear();
    sr.mc.nu.resize(1);  // una sola interazione per evento

    //caf::SRTrueInteraction &nu = mc.nu[0];

    // =========================
    // Riempimento dei campi
    // =========================
    sr.mc.nu[0].id = event_num;
    sr.mc.nu[0].E  = evt->Enu;
    sr.mc.nu[0].vtx.SetXYZ(evt->x, evt->y, evt->z);
    sr.mc.nu[0].momentum.SetXYZ(evt->pxnu, evt->pynu, evt->pznu);
    std::cout << "=== DEBUG FillTrueNu ===" << std::endl;
    std::cout << "Evento: posizione (" << evt->x << ", " << evt->y << ", " << evt->z << ") " << std::endl;
    std::cout << "Evento: Enutrue (" << evt->Enu << ") " << std::endl;


}

  void SANDRecoBranchFiller::FillRecoParticle(caf::StandardRecord &sr, const event* evt, int event_num) const {

    sr.common.ixn.sandreco.resize(1);
    
    auto &sandVec = sr.common.ixn.sandreco[0].part.sandreco;  // vettore di SRRecoParticle
    sandVec.clear();  

    size_t n = evt->particles.size();
    std::cout << "DEBUG: riempiendo sandreco per evento " << event_num
              << ", numero particelle = " << n << std::endl;

    for (size_t i = 0; i < n; ++i) {
        const auto &p = evt->particles[i];
        caf::SRRecoParticle r;

        r.primary = (p.primary != 0);  // cast int → bool
        r.pdg     = p.pdg_reco;
        r.E       = p.Ereco;
        r.start.SetXYZ(p.xreco, p.yreco, p.zreco);
        r.end.SetXYZ(p.xend, p.yend, p.zend);
        r.contained = true; 
        r.parent = p.parent_tid;  
        r.truth.push_back(caf::TrueParticleID{
            0,  // index dell'interaction (tipicamente 0 se c'è una sola interazione)
            r.primary ? caf::TrueParticleID::PartType::kPrimary
                      : caf::TrueParticleID::PartType::kSecondary,
            static_cast<int>(i)  // indice della true particle corrispondente in evt->particles
        });
 
        sandVec.push_back(r);

        /*std::cout << "DEBUG: particella " << i
                  << ", pdg=" << r.pdg
                  << ", primary(int)=" << p.primary
                  << ", primary(bool)=" << r.primary
                  << ", energia=" << r.E
                  << ", contained=" << r.contained
                  << ", parent_tid=" << r.parent
                  << "\n       start=(" << r.start.X() << "," << r.start.Y() << "," << r.start.Z() << ")"
                  << "\n       end  =(" << r.end.X() << "," << r.end.Y() << "," << r.end.Z() << ")"
                  << "\n       truth[ixn,type,part]=(" 
                  << 0 << "," 
                  << (r.primary ? "kPrimary" : "kSecondary") << "," 
                  << i << ")"
                  << std::endl;*/
    }
  }

  void SANDRecoBranchFiller::FillTruePrimaryParticles(caf::StandardRecord &sr,
                                                    const event* evt,
                                                    int event_num) const {

    // =========================
    // Reset truth branch
    // =========================
    sr.mc.nu.resize(1);
    caf::SRTrueInteraction &interaction = sr.mc.nu[0];

    size_t n = evt->particles.size();
    size_t nPrim = 0;

    std::cout << "==============================================" << std::endl;
    std::cout << "DEBUG: FillTruePrimaryParticles - evento " << event_num
              << ", numero particelle totali = " << n << std::endl;

    // =========================
    // RESET q0 e modq
    // =========================
    interaction.q0 = -9999;
    interaction.modq = -9999;

    // =========================
    // NEUTRINO MOMENTUM
    // =========================
    const double pnu_x = sr.mc.nu[0].momentum.x;
    const double pnu_y = sr.mc.nu[0].momentum.y;
    const double pnu_z = sr.mc.nu[0].momentum.z;

    const double Enu = sr.mc.nu[0].E;

    for (size_t i = 0; i < n; ++i) {
        const auto &p = evt->particles[i];

        if (p.primary != 1) continue;  // salta secondarie

        if (p.pdg_true != 13) continue; // 🔴 SOLO MUONE

        caf::SRTrueParticle rtp;

        // Campi comuni
        rtp.pdg  = p.pdg_true;
        rtp.G4ID = p.tid;
        rtp.time = p.ttrue;

        rtp.p.px = p.pxtrue;
        rtp.p.py = p.pytrue;
        rtp.p.pz = p.pztrue;
        rtp.p.E  = p.Etrue;

        rtp.start_pos.SetXYZ(p.xtrue, p.ytrue, p.ztrue);
        rtp.end_pos.SetXYZ(p.xend_true, p.yend_true, p.zend_true);

        std::cout << "Primary muon found:" << std::endl;
        std::cout << "  PDG = " << rtp.pdg << ", G4ID = " << rtp.G4ID << std::endl;

        // =========================
        // Q0 = Enu - Emu
        // =========================
        const double Emu = rtp.p.E;
        interaction.q0 = Enu - Emu;

        // =========================
        // q3 = |pnu - pmu|
        // =========================
        const double dq_x = pnu_x - rtp.p.px;
        const double dq_y = pnu_y - rtp.p.py;
        const double dq_z = pnu_z - rtp.p.pz;

        interaction.modq = std::sqrt(dq_x*dq_x + dq_y*dq_y + dq_z*dq_z);

        std::cout << "  q0 = " << interaction.q0
                  << "  q3 = " << interaction.modq << std::endl;

        interaction.prim.push_back(rtp);

        nPrim++;
        break; // IMPORTANTISSIMO: un solo muone
    }

    std::cout << "DEBUG: Totale primarie salvate = " << nPrim << std::endl;
    std::cout << "==============================================" << std::endl;
}

void SANDRecoBranchFiller::FillTrueSecondaryParticles(caf::StandardRecord &sr,
                                             const event* evt,
                                             int event_num) const {

    sr.mc.nu.resize(1);

    size_t n = evt->particles.size();
    caf::SRTrueInteraction &interaction = sr.mc.nu[0];

    std::cout << "==============================================" << std::endl;
    std::cout << "DEBUG: FillTrueParticles - evento "
              << event_num << std::endl;
    std::cout << "Numero particelle totali nell'evento = "
              << n << std::endl;

    for (size_t i = 0; i < n; ++i) {

        const auto &p = evt->particles[i];

        caf::SRTrueParticle rtp;

        // =========================
        // Campi comuni
        // =========================
        rtp.pdg            = p.pdg_true;
        rtp.G4ID           = p.tid;
        rtp.time           = p.ttrue;

       // Momento
        rtp.p.px = p.pxtrue;
        rtp.p.py = p.pytrue;
        rtp.p.pz = p.pztrue;
        rtp.p.E  = p.Etrue;

        rtp.start_pos.SetXYZ(p.xtrue,p.ytrue,p.ztrue);
        rtp.end_pos.SetXYZ(p.xend_true, p.yend_true, p.zend_true);


        if (p.primary == 0) {

            rtp.parent = p.parent_tid;

            interaction.sec.push_back(rtp);

            /*std::cout << "DEBUG SECONDARIA [" << i << "]"
                      << " | pdg=" << rtp.pdg
                      << " | G4ID=" << rtp.G4ID
                      << " | parent=" << rtp.parent
                      << " | time=" << rtp.time
                      << "\n    p = ("
                      << p.pxtrue << ", "
                      << p.pytrue << ", "
                      << p.pztrue << ", "
                      << p.Etrue << ")"
                      << "\n    start = ("
                      << p.xtrue << ", "
                      << p.ytrue << ", "
                      << p.ztrue << ")"
                      << "\n    end = ("
                      << p.xend_true << ", "
                      << p.yend_true << ", "
                      << p.zend_true << ")"
                      << std::endl;*/
        }
    }

  }


  /*void SANDRecoBranchFiller::FillRecoParticle(caf::StandardRecord &sr, std::vector<particle> particles) const {

    sr.common.ixn.sandreco.resize(1);
    
    sr.common.ixn.sandreco[0].vtx.SetXYZ(evt->x, evt->y, evt->z);
    std::cout << "=== DEBUG FillRecoNu ===" << std::endl;
    std::cout << "Evento: posizione (" << evt->x << ", " << evt->y << ", " << evt->z << ") " << std::endl;
    sr.common.ixn.sandreco[0].Enu.calo = evt->Enureco;
  }*/

  // todo: this is a placeholder
  std::deque<Trigger> SANDRecoBranchFiller::GetTriggers(int triggerType, bool beamOnly) const
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
std::deque<Trigger> SANDRecoBranchFiller::GetTriggers(int, bool) const { return {}; }
void SANDRecoBranchFiller::_FillRecoBranches(const Trigger &, caf::StandardRecord &, const cafmaker::Params &, const TruthMatcher *) const {}
}
#endif


