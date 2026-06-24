/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi, S. Lanzi
/// \date    Apr. 2022
///

#include "SANDRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

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

    fTriggers.push_back(Trigger{fSR->mc.nu[0].id, 0, static_cast<unsigned long>(e), 0u});
  }

  SetConfigured(true);
}

std::deque<Trigger>
SANDRecoBranchFiller<SANDRecoVersion::experimental>::GetTriggers(int triggerType,
                                                                 bool /*beamOnly*/) const {
  std::deque<Trigger> out;
  for (const auto& t : fTriggers) {
    if (triggerType < 0 || t.triggerType == triggerType) {
      out.push_back(t);
    }
  }
  return out;
}

void SANDRecoBranchFiller<SANDRecoVersion::experimental>::_FillRecoBranches(
    const Trigger& trigger, caf::StandardRecord& sr, const cafmaker::Params& /*par*/,
    const TruthMatcher* /*truthMatcher*/) const {
  auto it = fEvtToEntry.find(trigger.evtID);
  if (it == fEvtToEntry.end()) {
    LOG.WARNING() << GetName() << ": no entry for evtID=" << trigger.evtID << "\n";
    return;
  }

  fCAFTree->GetEntry(it->second);

  sr.common.ixn.sandreco  = fSR->common.ixn.sandreco;
  sr.common.ixn.nsandreco = fSR->common.ixn.nsandreco;
  sr.nd.sand.ixn          = fSR->nd.sand.ixn;
  sr.nd.sand.nixn         = fSR->nd.sand.nixn;
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
