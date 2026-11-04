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

// Comparison of the mc of the SAND file with sr.mc. Not compared: genieIdx (sandreco does not compute it).
namespace {

bool Close(float a, float b, double relTol = 1e-6) {  // default: float rounding
  return a == b || (std::isnan(a) && std::isnan(b)) || std::abs(a - b) <= relTol * std::max(std::abs(a), std::abs(b));
}

bool SameVector(const caf::SRVector3D& a, const caf::SRVector3D& b) {
  return Close(a.x, b.x) && Close(a.y, b.y) && Close(a.z, b.z);
}

bool SameID(const caf::TrueParticleID& a, const caf::TrueParticleID& b) {
  return a.ixn == b.ixn && a.type == b.type && a.part == b.part;
}

bool SameParticle(const caf::SRTrueParticle& a, const caf::SRTrueParticle& b) {
  return a.pdg == b.pdg && a.G4ID == b.G4ID && a.interaction_id == b.interaction_id && a.parent == b.parent
         && SameID(a.ancestor_id, b.ancestor_id) && SameID(a.parentID, b.parentID) && a.daughters == b.daughters
         && std::equal(a.daughtersID.begin(), a.daughtersID.end(), b.daughtersID.begin(), b.daughtersID.end(), SameID)
         && Close(a.p.px, b.p.px) && Close(a.p.py, b.p.py) && Close(a.p.pz, b.p.pz)
         && Close(a.p.E, b.p.E, 1e-3)  // sandreco uses the edep-sim mass, the CAFMaker the GENIE one
         && Close(a.time, b.time) && SameVector(a.start_pos, b.start_pos) && SameVector(a.end_pos, b.end_pos)
         && a.first_process == b.first_process && a.first_subprocess == b.first_subprocess
         && a.end_process == b.end_process && a.end_subprocess == b.end_subprocess;
}

bool SameParticles(const std::vector<caf::SRTrueParticle>& a, const std::vector<caf::SRTrueParticle>& b) {
  return std::equal(a.begin(), a.end(), b.begin(), b.end(), SameParticle);
}

bool SameDials(const std::vector<caf::SRDialWeights>& a, const std::vector<caf::SRDialWeights>& b) {
  return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](const auto& x, const auto& y) { return x.weights == y.weights; });
}

bool SameInteraction(const caf::SRTrueInteraction& a, const caf::SRTrueInteraction& b) {
  return a.id == b.id && a.pdg == b.pdg && a.pdgorig == b.pdgorig && a.iscc == b.iscc && a.mode == b.mode
         && a.targetPDG == b.targetPDG && a.hitnuc == b.hitnuc && Close(a.removalE, b.removalE) && Close(a.E, b.E)
         && SameVector(a.vtx, b.vtx) && SameVector(a.momentum, b.momentum) && a.isvtxcont == b.isvtxcont
         && Close(a.time, b.time) && Close(a.bjorkenX, b.bjorkenX) && Close(a.inelasticity, b.inelasticity)
         && Close(a.Q2, b.Q2) && Close(a.q0, b.q0) && Close(a.modq, b.modq) && Close(a.W, b.W) && Close(a.t, b.t)
         && a.ischarm == b.ischarm && a.isseaquark == b.isseaquark && a.resnum == b.resnum && Close(a.xsec, b.xsec)
         && Close(a.genweight, b.genweight) && Close(a.baseline, b.baseline) && SameVector(a.prod_vtx, b.prod_vtx)
         && SameVector(a.parent_dcy_mom, b.parent_dcy_mom) && a.parent_dcy_mode == b.parent_dcy_mode
         && a.parent_pdg == b.parent_pdg && Close(a.parent_dcy_E, b.parent_dcy_E) && Close(a.imp_weight, b.imp_weight)
         && a.generator == b.generator && a.genVersion == b.genVersion && a.nproton == b.nproton
         && a.nneutron == b.nneutron && a.npip == b.npip && a.npim == b.npim && a.npi0 == b.npi0
         && a.nprim == b.nprim && SameParticles(a.prim, b.prim) && a.nprefsi == b.nprefsi
         && SameParticles(a.prefsi, b.prefsi) && a.nsec == b.nsec && SameParticles(a.sec, b.sec)
         && Close(a.xsec_cvwgt, b.xsec_cvwgt) && SameDials(a.xsec_systs, b.xsec_systs)
         && SameDials(a.flux_systs, b.flux_systs);
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
