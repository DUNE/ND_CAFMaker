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

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

// Consistency check of the mc of the SAND file with sr.mc: interaction IDs and counts, particle PDG and Geant4 IDs
// in order, parent links, raw momenta and positions.
namespace {

bool Close(float a, float b) {
  return a == b || (std::isnan(a) && std::isnan(b)) || std::abs(a - b) <= 1e-6 * std::max(std::abs(a), std::abs(b));
}

bool SameVector(const caf::SRVector3D& a, const caf::SRVector3D& b) {
  return Close(a.x, b.x) && Close(a.y, b.y) && Close(a.z, b.z);
}

bool SameParticle(const caf::SRTrueParticle& a, const caf::SRTrueParticle& b) {
  return a.pdg == b.pdg && a.G4ID == b.G4ID && a.parent == b.parent && a.parentID.ixn == b.parentID.ixn
         && a.parentID.type == b.parentID.type && a.parentID.part == b.parentID.part && Close(a.p.px, b.p.px)
         && Close(a.p.py, b.p.py) && Close(a.p.pz, b.p.pz) && SameVector(a.start_pos, b.start_pos)
         && SameVector(a.end_pos, b.end_pos);
}

bool SameParticles(const std::vector<caf::SRTrueParticle>& a, const std::vector<caf::SRTrueParticle>& b) {
  return std::equal(a.begin(), a.end(), b.begin(), b.end(), SameParticle);
}

bool SameInteraction(const caf::SRTrueInteraction& a, const caf::SRTrueInteraction& b) {
  return a.id == b.id && a.nprim == b.nprim && a.nprefsi == b.nprefsi && a.nsec == b.nsec && SameVector(a.vtx, b.vtx)
         && SameVector(a.momentum, b.momentum) && SameParticles(a.prim, b.prim) && SameParticles(a.prefsi, b.prefsi)
         && SameParticles(a.sec, b.sec);
}

}  // namespace

namespace cafmaker {

SANDRecoBranchFiller::~SANDRecoBranchFiller() = default;

SANDRecoBranchFiller::SANDRecoBranchFiller(
    const std::string& cafFilename, std::string treeName)
    : IRecoBranchFiller("SAND") {
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

  const auto nEntries = fCAFTree->GetEntries();
  fTriggers.reserve(static_cast<std::size_t>(nEntries));

  SetConfigured(true);
}

std::deque<Trigger>
SANDRecoBranchFiller::GetTriggers(int triggerType,
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

void SANDRecoBranchFiller::_FillRecoBranches(
    const Trigger& trigger, caf::StandardRecord& sr, const cafmaker::Params& /*par*/,
    const TruthMatcher* /*truthMatcher*/) const {
	
  auto it = std::find(fTriggers.begin(), fTriggers.end(), trigger);
  if (it == fTriggers.end()) {
    LOG.FATAL() << GetName() << ": no trigger with evtID=" << trigger.evtID << "\n";
    abort();
  }
  const long entry = std::distance(fTriggers.begin(), it);
  CheckedGetEntry(fCAFTree, entry, "SAND reco event lookup");

  // The truth links of the SAND reco objects are indices into the SAND file's own mc and are copied as they are:
  // they are only valid if that mc is the same as sr.mc.
  const auto& sandNu = fSR->mc.nu;
  const auto& cafNu = sr.mc.nu;
  const auto [itSand, itCaf] = std::mismatch(sandNu.begin(), sandNu.end(), cafNu.begin(), cafNu.end(), SameInteraction);
  if (fSR->mc.nnu != sr.mc.nnu || itSand != sandNu.end() || itCaf != cafNu.end()) {
    std::stringstream ss;
    ss << GetName() << ": the mc of SAND entry " << entry << " differs from sr.mc at nu["
       << std::distance(sandNu.begin(), itSand) << "]\n";
    LOG.FATAL() << ss.str();
    throw std::runtime_error(ss.str());
  }

  sr.common.ixn.sandreco  = fSR->common.ixn.sandreco;
  sr.common.ixn.nsandreco = fSR->common.ixn.nsandreco;
  sr.nd.sand              = fSR->nd.sand;

}

} // namespace cafmaker
