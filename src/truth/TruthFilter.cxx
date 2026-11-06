#include "truth/TruthFilter.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "util/Logger.h"

namespace
{
  struct SecInfo
  {
    int  parent     = -1;     ///< index in sec of its parent, -1 if its parent is a primary or unknown
    int  generation = 0;      ///< 1 for the daughters of a primary, parent's + 1 otherwise
    bool keep       = false;  ///< whether the secondary survives the filter
    int  newIdx     = -1;     ///< index in sec after the filter, -1 if removed
  };

  using SecInfoByIxn = std::vector<std::vector<SecInfo>>;

  // ------------------------------------------------------------
  [[noreturn]] void Fail(const std::string & msg)
  {
    cafmaker::LOG_S("FilterTrueParticles()").FATAL() << msg;
    throw std::runtime_error(msg);
  }

  // ------------------------------------------------------------
  template <typename T, typename Pred>
  void EraseIf(std::vector<T> & v, Pred pred)
  {
    v.erase(std::remove_if(v.begin(), v.end(), pred), v.end());
  }

  // ------------------------------------------------------------
  bool IsSecondaryOf(const caf::TrueParticleID & id, std::size_t ixnIdx)
  {
    return id.type == caf::TrueParticleID::kSecondary && static_cast<std::size_t>(id.ixn) == ixnIdx;
  }

  // ------------------------------------------------------------
  double KineticEnergy(const caf::SRTrueParticle & part)
  {
    const double E = part.p.E;
    const double px = part.p.px, py = part.p.py, pz = part.p.pz;
    return E - std::sqrt(std::max(E * E - (px * px + py * py + pz * pz), 0.));
  }

  // ------------------------------------------------------------
  template <typename T, typename F>
  void VisitTruthIDs(std::vector<T> & recoObjects, F & fn)
  {
    for (T & obj : recoObjects) {
      for (caf::TrueParticleID & id : obj.truth) {
        fn(id);
      }
    }
  }

  // ------------------------------------------------------------
  template <typename F>
  void ForEachRecoTruthID(caf::StandardRecord & sr, F && fn)
  {
    // reco particles of the common reco interactions
    for (std::vector<caf::SRInteraction> * ixns : {&sr.common.ixn.dlp, &sr.common.ixn.pandora, &sr.common.ixn.sandreco}) {
      for (caf::SRInteraction & ixn : *ixns) {
        VisitTruthIDs(ixn.part.dlp, fn);
        VisitTruthIDs(ixn.part.pandora, fn);
        VisitTruthIDs(ixn.part.pida, fn);
        VisitTruthIDs(ixn.part.sandreco, fn);
      }
    }

    // ND-LAr, TMS, MINERvA
    for (caf::SRNDLArInt & ixn : sr.nd.lar.dlp) {
      VisitTruthIDs(ixn.tracks, fn);
      VisitTruthIDs(ixn.showers, fn);
    }
    for (caf::SRNDLArInt & ixn : sr.nd.lar.pandora) {
      VisitTruthIDs(ixn.tracks, fn);
      VisitTruthIDs(ixn.showers, fn);
    }
    for (caf::SRTMSInt & ixn : sr.nd.tms.ixn) {
      VisitTruthIDs(ixn.tracks, fn);
    }
    for (caf::SRMINERvAInt & ixn : sr.nd.minerva.ixn) {
      VisitTruthIDs(ixn.tracks, fn);
      VisitTruthIDs(ixn.showers, fn);
    }

    // the tracks and showers built by the ND matchers
    for (caf::SRNDTrackAssn & assn : sr.nd.trkmatch.extrap) {
      for (caf::TrueParticleID & id : assn.trk.truth) {
        fn(id);
      }
    }
    for (caf::SRNDShowerAssn & assn : sr.nd.shwmatch.extrap) {
      for (caf::TrueParticleID & id : assn.shw.truth) {
        fn(id);
      }
    }

    // SAND
    VisitTruthIDs(sr.nd.sand.grain.tracks, fn);
    VisitTruthIDs(sr.nd.sand.grain.showers, fn);
    VisitTruthIDs(sr.nd.sand.tracker.tracks, fn);
    VisitTruthIDs(sr.nd.sand.tracker.showers, fn);
    VisitTruthIDs(sr.nd.sand.ecal.clusters, fn);
    for (caf::SRSANDInt & ixn : sr.nd.sand.ixn) {
      for (caf::SRSANDAssn & assn : ixn.trkmatch) {
        for (caf::TrueParticleID & id : assn.trk.truth) {
          fn(id);
        }
      }
    }
  }

  // ------------------------------------------------------------
  void MarkRecoLinked(caf::StandardRecord & sr, SecInfoByIxn & info)
  {
    ForEachRecoTruthID(sr, [&info](const caf::TrueParticleID & id) {
      if (id.type != caf::TrueParticleID::kSecondary) {
        return;
      }
      if (id.ixn < 0 || static_cast<std::size_t>(id.ixn) >= info.size()
          || id.part < 0 || static_cast<std::size_t>(id.part) >= info[id.ixn].size()) {
        std::stringstream ss;
        ss << "A reco object links to secondary " << id.part << " of true interaction " << id.ixn
           << ", which does not exist: cannot filter the true particles of this entry\n";
        Fail(ss.str());
      }
      info[id.ixn][id.part].keep = true;
    });
  }

  // ------------------------------------------------------------
  void FillParents(const caf::SRTrueInteraction & ixn, std::vector<SecInfo> & info)
  {
    for (std::size_t i = 0; i < ixn.sec.size(); ++i) {
      const caf::TrueParticleID & parentID = ixn.sec[i].parentID;
      if (parentID.type != caf::TrueParticleID::kSecondary) {
        continue;
      }
      if (parentID.part < 0 || static_cast<std::size_t>(parentID.part) >= ixn.sec.size()) {
        std::stringstream ss;
        ss << "Secondary " << i << " of true interaction " << ixn.id << " has parent secondary " << parentID.part
           << ", which does not exist: cannot filter the true particles of this entry\n";
        Fail(ss.str());
      }
      info[i].parent = parentID.part;
    }
  }

  // ------------------------------------------------------------
  void FillGenerations(const caf::SRTrueInteraction & ixn, std::vector<SecInfo> & info)
  {
    for (std::size_t i = 0; i < info.size(); ++i) {
      if (info[i].generation > 0) {
        continue;
      }
      std::vector<int> chain;
      int current = static_cast<int>(i);
      while (current >= 0 && info[current].generation == 0) {
        if (chain.size() > info.size()) {
          std::stringstream ss;
          ss << "The parents of secondary " << i << " of true interaction " << ixn.id << " form a loop: "
             << "cannot filter the true particles of this entry\n";
          Fail(ss.str());
        }
        chain.push_back(current);
        current = info[current].parent;
      }

      int gen = (current >= 0) ? info[current].generation : 0;

      for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
        info[*it].generation = ++gen;
      }
    }
  }

  // ------------------------------------------------------------
  void ApplyCuts(const std::vector<caf::SRTrueParticle> & sec, std::vector<SecInfo> & info,
                 double minKE_GeV, int maxGeneration)
  {
    for (std::size_t i = 0; i < sec.size(); ++i) {
      const double KE = KineticEnergy(sec[i]);
      const bool passKE = minKE_GeV <= 0 || std::isnan(KE) || KE >= minKE_GeV;
      const bool passGeneration = maxGeneration < 0 || info[i].generation <= maxGeneration;
      if (passKE && passGeneration) {
        info[i].keep = true;
      }
    }
  }

  // ------------------------------------------------------------
  void KeepAncestors(std::vector<SecInfo> & info)
  {
    for (const SecInfo & sec : info) {
      if (!sec.keep) {
        continue;
      }
      for (int parent = sec.parent; parent >= 0 && !info[parent].keep; parent = info[parent].parent) {
        info[parent].keep = true;
      }
    }
  }

  // ------------------------------------------------------------
  std::size_t AssignNewIndices(std::vector<SecInfo> & info)
  {
    int nKept = 0;
    for (SecInfo & sec : info) {
      if (sec.keep) {
        sec.newIdx = nKept++;
      }
    }
    return nKept;
  }

  // ------------------------------------------------------------
  void RemoveSecondaries(caf::SRTrueInteraction & ixn, const std::vector<SecInfo> & info, std::size_t nKept)
  {
    for (std::size_t i = 0; i < ixn.sec.size(); ++i) {
      const int newIdx = info[i].newIdx;
      if (newIdx >= 0 && static_cast<std::size_t>(newIdx) != i) {
        ixn.sec[newIdx] = std::move(ixn.sec[i]);
      }
    }
    ixn.sec.resize(nKept);
    ixn.nsec = nKept;
  }

  // ------------------------------------------------------------
  void RemapID(caf::TrueParticleID & id, std::size_t ixnIdx, const std::vector<SecInfo> & info)
  {
    if (IsSecondaryOf(id, ixnIdx)) {
      id.part = info.at(id.part).newIdx;
    }
  }

  // ------------------------------------------------------------
  void RemapParticleLinks(caf::SRTrueParticle & particle, std::size_t ixnIdx,
                          const std::vector<SecInfo> & info, const std::set<unsigned int> & removedG4IDs)
  {
    RemapID(particle.parentID, ixnIdx, info);

    for (caf::TrueParticleID & id : particle.daughtersID) {
      RemapID(id, ixnIdx, info);
    }
    EraseIf(particle.daughtersID, [ixnIdx](const caf::TrueParticleID & id) { return IsSecondaryOf(id, ixnIdx) && id.part < 0; });

    EraseIf(particle.daughters, [&removedG4IDs](unsigned int g4ID) { return removedG4IDs.count(g4ID) > 0; });
  }

  // ------------------------------------------------------------
  void RemapRecoLinks(caf::StandardRecord & sr, const SecInfoByIxn & info)
  {
    ForEachRecoTruthID(sr, [&info](caf::TrueParticleID & id) {
      if (id.type == caf::TrueParticleID::kSecondary) {
        id.part = info[id.ixn][id.part].newIdx;
      }
    });
  }
}

namespace cafmaker
{
  void FilterTrueParticles(caf::StandardRecord & sr, double minKE_GeV, int maxGeneration)
  {
    if (minKE_GeV <= 0 && maxGeneration < 0) {
      return;
    }

    std::vector<caf::SRTrueInteraction> & nus = sr.mc.nu;
    SecInfoByIxn info(nus.size());
    for (std::size_t ixnIdx = 0; ixnIdx < nus.size(); ++ixnIdx) {
      info[ixnIdx].resize(nus[ixnIdx].sec.size());
    }

    MarkRecoLinked(sr, info);

    std::size_t nRemoved = 0;
    for (std::size_t ixnIdx = 0; ixnIdx < nus.size(); ++ixnIdx) {
      caf::SRTrueInteraction & ixn = nus[ixnIdx];
      std::vector<SecInfo> & ixnInfo = info[ixnIdx];

      FillParents(ixn, ixnInfo);
      FillGenerations(ixn, ixnInfo);
      ApplyCuts(ixn.sec, ixnInfo, minKE_GeV, maxGeneration);
      KeepAncestors(ixnInfo);
      const std::size_t nKept = AssignNewIndices(ixnInfo);
      if (nKept == ixn.sec.size()) {
        continue;
      }

      std::set<unsigned int> removedG4IDs;
      for (std::size_t i = 0; i < ixn.sec.size(); ++i) {
        if (!ixnInfo[i].keep) {
          removedG4IDs.insert(static_cast<unsigned int>(ixn.sec[i].G4ID));
        }
      }
      nRemoved += ixn.sec.size() - nKept;

      RemoveSecondaries(ixn, ixnInfo, nKept);
      for (std::vector<caf::SRTrueParticle> * particles : {&ixn.prim, &ixn.prefsi, &ixn.sec}) {
        for (caf::SRTrueParticle & particle : *particles) {
          RemapParticleLinks(particle, ixnIdx, ixnInfo, removedG4IDs);
        }
      }
    }

    if (nRemoved > 0) {
      RemapRecoLinks(sr, info);
    }

    LOG_S("FilterTrueParticles()").VERBOSE() << "Removed " << nRemoved << " secondary true particles"
                                             << " (min kinetic energy " << minKE_GeV << " GeV, max generation "
                                             << maxGeneration << ")\n";
  }
}
