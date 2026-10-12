/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi, S. Lanzi
/// \date    Apr. 2022
///

#include "SANDRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "truth/FillTruth.h"

#include <TFile.h>
#include <TTree.h>

#ifdef ENABLE_SANDRECO_LEGACY
#include "struct.h"
#endif

namespace cafmaker {

SANDRecoBranchFiller<SANDRecoVersion::experimental>::~SANDRecoBranchFiller() = default;

SANDRecoBranchFiller<SANDRecoVersion::experimental>::SANDRecoBranchFiller(
    const std::string& cafFilename, std::string treeName)
    : IRecoBranchFiller("SANDExperimental") {
  if (cafFilename.empty()) {
    LOG.WARNING() << GetName() << ": filename is empty\n";
    SetConfigured(false);
    return;
  }

  fFile.reset(TFile::Open(cafFilename.c_str()));
  if (!fFile || fFile->IsZombie()) {
    LOG.WARNING() << GetName() << ": could not open '" << cafFilename << "'\n";
    SetConfigured(false);
    return;
  }

  fCAFTree = fFile->Get<TTree>(treeName.c_str());
  if (!fCAFTree) {
    LOG.WARNING() << GetName() << ": tree '" << treeName << "' not found in '" << cafFilename
                  << "'\n";
    SetConfigured(false);
    return;
  }

  fCAFTree->SetBranchAddress("rec", &fSR);

  // Single pass: build evtID → entry map and trigger list
  const auto nEntries = fCAFTree->GetEntries();
  fTriggers.reserve(static_cast<std::size_t>(nEntries));
  for (Long64_t e{}; e != nEntries; ++e) {
    fCAFTree->GetEntry(e);
    if (fSR->mc.nnu == 0) {
      continue;
    }

    for (std::size_t i{}; i != fSR->mc.nnu; ++i) {
      fEvtToEntry[fSR->mc.nu[i].id] = e;
    }

    // fTriggers.push_back(Trigger{e, 0, static_cast<unsigned long>(e), 0u});
  }

  SetConfigured(true);
}

std::deque<Trigger>
SANDRecoBranchFiller<SANDRecoVersion::experimental>::GetTriggers(int triggerType,
                                                                 bool /*beamOnly*/) const {

  // This implementation is copied from the TMS BranchFiller. 
  // Without this, there was a mismatch between the triggers from different branches.
  std::deque<Trigger> triggers;

  if (fTriggers.empty())
  {
    LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "' from " << fCAFTree->GetEntries() << " TMS Reco_Tree:\n";
  
    fTriggers.reserve(fCAFTree->GetEntries());
      for (int entry = 0; entry < fCAFTree->GetEntries(); entry++)
      {
        fCAFTree->GetEntry(entry);

        const Trigger * prev_trig = fTriggers.empty() ? nullptr : &fTriggers.back();
        fTriggers.emplace_back();
        Trigger & trig = fTriggers.back();

        trig.evtID = entry;
        trig.triggerType = 1; // TODO real number?

        if (prev_trig)
        {
          trig.triggerTime_ns = prev_trig->triggerTime_ns + 2E8 ;
          trig.triggerTime_s = prev_trig->triggerTime_s + 1; // TODO: Pull the 1.2 from correct place in file
          if (trig.triggerTime_ns >= 1E9) // If we have 1s worth of ns then add 1s and remove 1s worth of ns
          {
            trig.triggerTime_s += 1;
            trig.triggerTime_ns -= 1E9;
          }
        }
        else
        {
          trig.triggerTime_ns = 0;
          trig.triggerTime_s = 0;
        }

        LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
                      << ", triggerType=" << trig.triggerType
                      << ", triggerTime_s=" << trig.triggerTime_s
                      << ", triggerTime_ns=" << trig.triggerTime_ns
                      << "\n";
      }
    }

    for (const Trigger & trigger : fTriggers)
    {
      if (triggerType < 0 || triggerType == trigger.triggerType)
      {
        triggers.push_back(trigger);
      }
    }

    return triggers;
}

void SANDRecoBranchFiller<SANDRecoVersion::experimental>::_FillRecoBranches(
    const Trigger& trigger, caf::StandardRecord& sr, const cafmaker::Params& /*par*/,
    const TruthMatcher* truthMatcher) const {

    if (!sr.mc.nu.empty() || sr.mc.nnu != 0) {
    std::ostringstream msg;
    msg << GetName() << ": stopping because sr.mc already contains true interactions (sr.mc.nu.size() = "
        << sr.mc.nu.size() << ", sr.mc.nnu = " << sr.mc.nnu << ") while the SAND reco filler is about to fill it. "
        << "The SAND filler copies sr.mc from the sandreco CAF and must be the first filler to run on a "
        << "StandardRecord: running after another filler (or the TruthMatcher) would overwrite the true "
        << "interactions they created and invalidate the genieIdx / TrueParticleID links already stored in the "
        << "reco branches. Check the order of the reco fillers and that no other filler creates true "
        << "interactions for the same trigger.";
    LOG.FATAL() << msg.str() << "\n";
    throw std::runtime_error(msg.str());
  }
	
  auto it = std::find(fTriggers.begin(), fTriggers.end(), trigger);
  if (it == fTriggers.end()) {
    LOG.FATAL() << GetName() << ": no trigger with evtID=" << trigger.evtID << "\n";
    abort();
  }
  const long entry = std::distance(fTriggers.begin(), it);
  fCAFTree->GetEntry(entry);

  sr.common.ixn.sandreco  = fSR->common.ixn.sandreco;
  sr.common.ixn.nsandreco = fSR->common.ixn.nsandreco;
  sr.nd.sand.ixn          = fSR->nd.sand.ixn;
  sr.nd.sand.nixn         = fSR->nd.sand.nixn;
  sr.nd.sand.ecal         = fSR->nd.sand.ecal;
  sr.nd.sand.tracker      = fSR->nd.sand.tracker;
  sr.nd.sand.grain        = fSR->nd.sand.grain;
  sr.mc                   = fSR->mc;

	  truthMatcher->StoreGENIEEvents(sr);
	  std::cout << "Stored Genie events" << std::endl;
}

#ifdef ENABLE_SANDRECO_LEGACY

SANDRecoBranchFiller<SANDRecoVersion::legacy>::~SANDRecoBranchFiller() = default;

SANDRecoBranchFiller<SANDRecoVersion::legacy>::SANDRecoBranchFiller(const std::string& recoFilename)
    : IRecoBranchFiller("SAND") {
  if (recoFilename.empty()) {
    LOG.WARNING() << GetName() << ": filename is empty\n";
    SetConfigured(false);
    return;
  }

  fFile.reset(TFile::Open(recoFilename.c_str()));
  if (!fFile || fFile->IsZombie()) {
    LOG.WARNING() << GetName() << ": could not open '" << recoFilename << "'\n";
    SetConfigured(false);
    return;
  }

  fRecoTree = fFile->Get<TTree>("tReco");
  if (!fRecoTree) {
    LOG.WARNING() << GetName() << ": tree 'tReco' not found\n";
    fFile->ls();
    SetConfigured(false);
    return;
  }

  // evtID = entry index (placeholder — no event ID in legacy reco format)
  const long nEntries = fRecoTree->GetEntries();
  fTriggers.reserve(nEntries);
  for (long e{}; e != nEntries; ++e) {
    Trigger t;
    t.evtID          = e;
    t.triggerType    = 0;
    t.triggerTime_s  = static_cast<unsigned long>(e);
    t.triggerTime_ns = 0;
    fTriggers.push_back(t);
  }

  SetConfigured(true);
}

std::deque<Trigger>
SANDRecoBranchFiller<SANDRecoVersion::legacy>::GetTriggers(int triggerType,
                                                           bool /*beamOnly*/) const {
  std::deque<Trigger> out;
  for (const auto& t : fTriggers) {
    if (triggerType < 0 || t.triggerType == triggerType) {
      out.push_back(t);
    }
  }
  return out;
}

void SANDRecoBranchFiller<SANDRecoVersion::legacy>::_FillRecoBranches(
    const Trigger& trigger, caf::StandardRecord& sr, const cafmaker::Params& /*par*/,
    const TruthMatcher* truthMatcher) const {
  auto it = std::find(fTriggers.begin(), fTriggers.end(), trigger);
  if (it == fTriggers.end()) {
    LOG.FATAL() << GetName() << ": no trigger with evtID=" << trigger.evtID << "\n";
    abort();
  }
  const long entry = std::distance(fTriggers.begin(), it);

  // The legacy tree stores vector<cluster>* and vector<track>* branches
  std::vector<cluster>* pcl = nullptr;
  std::vector<track>* ptr   = nullptr;
  fRecoTree->SetBranchAddress("cluster", &pcl);
  fRecoTree->SetBranchAddress("track", &ptr);
  fRecoTree->GetEntry(entry);

  if (pcl) {
    FillECalClusters(truthMatcher, sr, *pcl);
  }
  if (ptr) {
    FillTracks(truthMatcher, sr, *ptr);
  }
}

void SANDRecoBranchFiller<SANDRecoVersion::legacy>::FillECalClusters(
    const TruthMatcher* /*truthMatcher*/, caf::StandardRecord& sr, std::vector<cluster>& cl) const {
  const std::size_t n = cl.size();
  sr.nd.sand.ixn.resize(1);
  auto& ixn     = sr.nd.sand.ixn[0];
  ixn.nclusters = n;
  ixn.ECALClusters.resize(n);

  for (std::size_t i{}; i != n; ++i) {
    auto& c = ixn.ECALClusters[i];
    c.E     = cl[i].e;
    c.position.SetXYZ(cl[i].x, cl[i].y, cl[i].z);
    c.var_position.SetXYZ(cl[i].varx, cl[i].vary, cl[i].varz);
    c.time = cl[i].t;
    c.start.SetXYZ(cl[i].ax, cl[i].ay, cl[i].az);
    c.direction.SetXYZ(cl[i].sx, cl[i].sy, cl[i].sz);
    c.num_cells = cl[i].reco_cells.size();
    c.id        = cl[i].tid;
  }
}

void SANDRecoBranchFiller<SANDRecoVersion::legacy>::FillTracks(const TruthMatcher* /*truthMatcher*/,
                                                               caf::StandardRecord& sr,
                                                               std::vector<track>& tr) const {
  const std::size_t n = tr.size();
  sr.nd.sand.ixn.resize(1);
  auto& ixn   = sr.nd.sand.ixn[0];
  ixn.ntracks = n;
  ixn.tracks.resize(n);

  for (std::size_t i{}; i != n; ++i) {
    ixn.tracks[i].start.SetXYZ(tr[i].x0, tr[i].y0, tr[i].z0);
    ixn.tracks[i].qual = tr[i].chi2_cr;
  }
}

#endif // ENABLE_SANDRECO_LEGACY

} // namespace cafmaker
