/// \file FillTruth.cxx
///
/// Fill truth branches.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>, based on code by C. Marshall <chris.marshall@rochester.edu>
/// \date    Jan. 2022

#include "FillTruth.h"

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// GENIE
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"

// Standard Record format
#include "duneanaobj/StandardRecord/StandardRecord.h"

// ND_CAFMaker
#include "CAF.h"
#include "Params.h"

/// duneanaobj not guaranteed to be the same as GENIE scattering types
caf::ScatteringMode GENIE2CAF(genie::EScatteringType sc)
{
  switch(sc)
  {
    case genie::kScQuasiElastic:
      return caf::kQE;

    case genie::kScSingleKaon:
      return caf::kSingleKaon;

    case genie::kScDeepInelastic:
      return caf::kDIS;

    case genie::kScResonant:
      return caf::kRes;

    case genie::kScCoherentProduction:
      return caf::kCoh;

    case genie::kScDiffractive:
      return caf::kDiffractive;

    case genie::kScNuElectronElastic:
      return caf::kNuElectronElastic;

    case genie::kScInverseMuDecay:
      return caf::kInvMuonDecay;

    case genie::kScAMNuGamma:
      return caf::kAMNuGamma;

    case genie::kScMEC:
      return caf::kMEC;

    case genie::kScCoherentElastic:
      return caf::kCohElastic;

    case genie::kScInverseBetaDecay:
      return caf::kInverseBetaDecay;

    case genie::kScGlashowResonance:
      return caf::kGlashowResonance;

    case genie::kScIMDAnnihilation:
      return caf::kIMDAnnihilation;

    case genie::kScPhotonCoherent:
      return caf::kPhotonCoh;

    case genie::kScPhotonResonance:
      return caf::kPhotonRes;

    case genie::kScDarkMatterElastic:
      return caf::kDarkMatterElastic;

    case genie::kScDarkMatterDeepInelastic:
      return caf::kDarkMatterDIS;

    case genie::kScDarkMatterElectron:
      return caf::kDarkMatterElectron;

    case genie::kScUnknown: [[fallthrough]];
    case genie::kScNull:
      return caf::kUnknownMode;

    default:
      std::cerr << "Unrecognized GENIE scattering mode: " << sc << "\n";
      abort();
  }

}

namespace cafmaker
{
  template <>
  void ValidateOrCopy<double, float>(const double & input, float & target, const float & unsetVal)
  {
    const auto cmp = [](const double & a, const float &b) -> bool { return std::abs(static_cast<float>(a) - b) < 1e-6; };
    const auto assgn = [](const double & a, float & b) {  b = a; };
    return ValidateOrCopy<double, float>(input, target, unsetVal, cmp, assgn);
  }


// ------------------------------------------------------------
  TruthMatcher::TruthMatcher(TTree * contGTree,
                             TTree * uncontGTree,
                             const genie::NtpMCEventRecord *gEvt)
    : fContNuGTree(contGTree), fUncontNuGTree(uncontGTree), fGEvt(gEvt)
  {}

  // --------------------------------------------------------------
  void TruthMatcher::FillInteraction(caf::SRTrueInteraction& nu, const genie::NtpMCEventRecord * gEvt)
  {

    genie::EventRecord * event = gEvt->event;
    genie::Interaction * in = event->Summary();

    nu.id = gEvt->hdr.ievent;  //todo: need to make sure this ID is the one we get all the way out the other end from det sim

    // todo: GENIE doesn't know about the detector geometry.  do we just leave these for the passthrough from G4?
//    TLorentzVector vtx = *(event->Vertex());
//    nu.vtx = vtx.Vect();
//      nu.isvtxcont =
//    nu.time = vtx.T();

    nu.pdg = in->InitState().ProbePdg();
    nu.pdgorig = in->InitState().ProbePdg(); // fill this for similarity with FD, but no oscillations

    nu.iscc = in->ProcInfo().IsWeakCC();
    nu.mode = GENIE2CAF(in->ProcInfo().ScatteringTypeId());
    nu.targetPDG = in->InitState().Tgt().Pdg();
    nu.hitnuc = in->InitState().Tgt().HitNucPdg();

    //todo: get this from Hugh G or somebody who will get it right
    // nu.removalE =

    TLorentzVector lepP4;
    TLorentzVector nuP4nuc = *(in->InitState().GetProbeP4(genie::kRfHitNucRest));
    TLorentzVector nuP4 = *(in->InitState().GetProbeP4(genie::kRfLab));
    nu.E = in->InitState().ProbeE(genie::kRfLab);
    nu.momentum = nuP4.Vect();

    // true 4-momentum transfer
    TLorentzVector q = nuP4-lepP4;

    // Q2, W, x, y frequently do not get filled in GENIE Kinematics object, so calculate manually
    const double Mnuc = 0.939;  // average nucleon mass
    nu.Q2 = -q.Mag2();
    nu.q0 = q.E();
    nu.modq = q.Vect().Mag();
    nu.W = sqrt(Mnuc*Mnuc + 2.*nu.q0*Mnuc + q.Mag2()); // "Wexp"
    nu.bjorkenX = nu.Q2/(2*Mnuc*nu.q0);
    nu.inelasticity = nu.q0/nu.E;
    if (nu.mode == caf::kCoh || nu.mode == caf::kDiffractive)
      nu.t = in->Kine().t();

    nu.ischarm = in->ExclTag().IsCharmEvent();
    nu.isseaquark = in->ProcInfo().IsDeepInelastic() && in->InitState().Tgt().HitSeaQrk();
    if (nu.mode == caf::kRes)
      nu.resnum = static_cast<int>(in->ExclTag().Resonance());

    nu.xsec = event->XSec();
    nu.genweight = event->Weight();

    // loop truth particles

    for (int j=0; j< event->GetEntries(); j++)
    {
      auto p = dynamic_cast<const genie::GHepParticle *>((*event)[j]);

      if( p->Status() != genie::EGHepStatus::kIStStableFinalState
          && p->Status() != genie::EGHepStatus::kIStHadronInTheNucleus) continue;

      caf::SRTrueParticle part;
      part.pdg = p->Pdg();
      part.interaction_id = nu.id;
      part.p = *p->P4();
      // todo: original GENIE events haven't been transferred to det coordinates
      //       or bunched in spill time structure yet.  just leave out?
//      part.time = nu.time;
//      part.start_pos = p->X4()->Vect();

      // remaining fields need to be filled in with post-G4 info

      if( p->Status() == genie::EGHepStatus::kIStStableFinalState )
      {
        // note: we leave part.id unset since it won't match with the G4 values
        // (edep-sim numbers them all sequentially from 0 throughout the whole file)
        // and the pass-through value is more useful

        nu.prim.push_back(std::move(part));
        nu.nprim++;

        if( p->Pdg() == 2212 ) nu.nproton++;
        else if( p->Pdg() == 2112 ) nu.nneutron++;
        else if( p->Pdg() ==  211 ) nu.npip++;
        else if( p->Pdg() == -211 ) nu.npim++;
        else if( p->Pdg() ==  111 ) nu.npi0++;
      }
      else // kIStHadronInTheNucleus
      {
        std::cout << "  particle " << j << " with pdg = " << p->Pdg() << " and energy = " << p->E() << " has GENIE status " << p->Status() << "\n";
        nu.prefsi.push_back(std::move(part));
        nu.nprefsi++;
      }

    }

    // todo: need to fill the flux variables in.  for 2x2, info should come from a genie::flux::GNuMIFluxPassThroughInfo object created by the flux driver.
    //       for DUNE beam, I assume there's an analogous thing?
    //nu.baseline =
  //  nu.prod_vtx       = ;              ///< Neutrino production vertex [cm; beam coordinates]
  //  nu.parent_dcy_mom =  ;        ///< Neutrino parent momentum at decay [GeV; beam coordinates]
  //  nu.parent_dcy_mode = ;  ///< Parent hadron/muon decay mode
  //  nu.parent_pdg      = ;   ///< PDG Code of parent particle ID
  //  nu.parent_dcy_E    = ; ///< Neutrino parent energy at decay [GeV]
  //  nu.imp_weight      = ; ///< Importance weight from flux file


    // Add DUNErw weights to the CAF
    nu.xsec_cvwgt = 1;

    // fixme: the following is disabled until DIRT-II finishes on model + uncertainty decisions
    //systtools::event_unit_response_w_cv_t resp = rh.GetEventVariationAndCVResponse(*event);
    //for( const systtools::VarAndCVResponse& it : resp ) {
    //  // Need begin/end to convert double to float
    //  sr.xsSyst_wgt.emplace_back(it.responses.begin(), it.responses.end());
    //  sr.cvwgt.push_back(it.CV_response);
    //  sr.total_xsSyst_cv_wgt *= it.CV_response;
    //}
  }

  // ------------------------------------------------------------
  caf::SRTrueParticle &
  TruthMatcher::GetTrueParticle(caf::StandardRecord &sr, int ixnID, int G4ID, bool isPrimary, bool createNew) const
  {
    return GetTrueParticle(sr, GetTrueInteraction(sr, ixnID, false), G4ID, isPrimary, createNew);
  }

  // ------------------------------------------------------------
  // helper class used only to avoid building and tearing down a new lambda every iteration
  namespace
  {
    struct SRPartCmp
    {
      int partID;
      bool operator()(const caf::SRTrueParticle & part) const { return part.G4ID == partID; }
    };
  }

  // ------------------------------------------------------------
  caf::SRTrueParticle &
  TruthMatcher::GetTrueParticle(caf::StandardRecord &sr, caf::SRTrueInteraction& ixn, int G4ID, bool isPrimary, bool createNew) const
  {
    static SRPartCmp srPartCmp;
    srPartCmp.partID = G4ID;
    return GetTrueParticle(sr, ixn, srPartCmp, isPrimary, createNew);
  }

  // ------------------------------------------------------------
  caf::SRTrueParticle &TruthMatcher::GetTrueParticle(caf::StandardRecord &sr,
                                                     caf::SRTrueInteraction &ixn,
                                                     std::function<bool(const caf::SRTrueParticle &)> cmp,
                                                     bool isPrimary,
                                                     bool createNew) const
  {
    caf::SRTrueParticle * part = nullptr;
    std::vector<caf::SRTrueParticle> & collection = (isPrimary) ? ixn.prim : ixn.sec;
    int & counter = (isPrimary) ? ixn.nprim : ixn.nsec;
    if ( auto itPart = std::find_if(collection.begin(), collection.end(), cmp);
         itPart == collection.end() )
    {
      if (!createNew)
        throw std::runtime_error("True particle with interaction ID " + std::to_string(ixn.id)
                                 + " was not found in the " + std::string(isPrimary ? "primary" : "secondary") + " true particle collection");
      else
        std::cout << "  made a new SRTrueParticle in " << (isPrimary ? "prim" : "sec") << " collection \n";

      collection.emplace_back();
      counter++;

      part = &collection.back();
      part->interaction_id = ixn.id;  // todo: eventually this should be ixn.track_id
    }
    else
    {
      std::cout << "    --> found previously created SRTrueParticle.  Returning that.\n";
      part = &(*itPart);
    }

    return *part;
  }

  // ------------------------------------------------------------
  caf::SRTrueInteraction &
  TruthMatcher::GetTrueInteraction(caf::StandardRecord &sr,
                                   std::function<bool(const caf::SRTrueInteraction &)> srTrueIxnCmp,
                                   std::function<bool(const genie::NtpMCEventRecord*)> genieCmp,
                                   bool createNew) const
  {
    caf::SRTrueInteraction * ixn = nullptr;

    // if we can't find a SRTrueInteraction with matching ID, we may need to make a new one
    if ( auto itIxn = std::find_if(sr.mc.nu.begin(), sr.mc.nu.end(), srTrueIxnCmp);
         itIxn == sr.mc.nu.end() )
    {
      if (!createNew)
        throw std::runtime_error("True interaction not found in this StandardRecord");

      std::cout << "    creating new SRTrueInteraction.  Trying to match to a GENIE event...\n";

      // todo: when we finally have interaction IDs that conform to the edep-sim vertexID convention,
      //       we can use those to only search the correct GENIE true (contained or uncontained).
      //       for now, though, we need to look at both.
      bool foundEv = false;
      int wrappedIdx = -1;
      std::cout << "     examining GENIE records...\n ";
      for (TTree * tree: {fContNuGTree, fUncontNuGTree})
      {
        if (!tree)
          continue;

        // if we already found the event in the first tree
        if (foundEv)
          break;

        long int lastReadEvt = tree->GetReadEvent();
        if (tree->GetReadEvent() < 0)
        {
          tree->GetEntry(0);
          lastReadEvt = -1;
        }

        // the most likely place to find the matching event is just beyond wherever we currently are,
        // so look there first, then loop back around to consider events previous to where we were
        std::cout << "       last read evt = " << lastReadEvt << ", total entries = " << tree->GetEntries() << "\n";
        for (long int evtIdx = lastReadEvt + 1; evtIdx % tree->GetEntries() != lastReadEvt; evtIdx++)
        {
          wrappedIdx = evtIdx % tree->GetEntries();
//          std::cout << " " << wrappedIdx;
          tree->GetEntry(wrappedIdx);

          if (lastReadEvt < 0)
            lastReadEvt = evtIdx;

          if (genieCmp(fGEvt))
          {
            foundEv = true;
            break;
          }
        }
        std::cout << "\n";
      } // for (tree)

      if (HaveGENIE() && !foundEv)
        throw std::runtime_error("Could not locate GENIE event record!");

      sr.mc.nu.emplace_back();
      sr.mc.nnu++;

      ixn = &sr.mc.nu.back();
      if (HaveGENIE())
      {
        std::cout << "      --> GENIE record found.  copying...\n";
        ixn->genieIdx = wrappedIdx;
        FillInteraction(*ixn, fGEvt);
      }
      else
        std::cout << "      --> no matching GENIE interaction found.  New empty SRTrueInteraction will be returned.\n";
    } // if ( didn't find a matching SRTrueInteraction )
    else
    {
      std::cout << "   Found previously created SRTrueInteraction.  Returning that.\n";
      ixn = &(*itIxn);
    }

    return *ixn;
  }

  // ------------------------------------------------------------
  // helper classes used only to avoid building and tearing down a new lambda every iteration
  namespace
  {
    struct SRCmp
    {
      int ixnID;
      bool operator()(const caf::SRTrueInteraction & ixn) const { return ixn.id == ixnID; }
    };

    struct GENIECmp
    {
      int ixnID;
      bool operator()(const genie::NtpMCEventRecord * gEvt) const { return static_cast<int>(gEvt->hdr.ievent) == ixnID; }
    };
  }

  // ------------------------------------------------------------
  caf::SRTrueInteraction & TruthMatcher::GetTrueInteraction(caf::StandardRecord &sr, int ixnID, bool createNew) const
  {
    static SRCmp srCmp;
    static GENIECmp genieCmp;
    try
    {
      srCmp.ixnID = genieCmp.ixnID = ixnID;
      return GetTrueInteraction(sr, srCmp, genieCmp, createNew);
    }
    catch (std::runtime_error & err)
    {
      // just re-throw with a nicer error message
      throw std::runtime_error("Could not locate GENIE event record with ID = " + std::to_string(ixnID));
    }
  }

}
