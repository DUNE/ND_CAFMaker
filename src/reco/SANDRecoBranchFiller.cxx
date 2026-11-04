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

  sr.common.ixn.sandreco  = fSR->common.ixn.sandreco;
  sr.common.ixn.nsandreco = fSR->common.ixn.nsandreco;
  sr.nd.sand              = fSR->nd.sand;

}

} // namespace cafmaker
