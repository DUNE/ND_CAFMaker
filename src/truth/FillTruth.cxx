/// \file FillTruth.cxx
///
/// Fill truth branches.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>, based on code by C. Marshall <chris.marshall@rochester.edu>
/// \date    Jan. 2022

#include "FillTruth.h"

#include <map>
#include <regex>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// GENIE
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/GHEP/GHepParticle.h"

// Standard Record format
#include "duneanaobj/StandardRecord/StandardRecord.h"

// ND_CAFMaker
#include "CAF.h"
#include "Params.h"
#include "util/FloatMath.h"

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
  void ValidateOrCopy<double, float>(const double & input, float & target, const float & unsetVal, const std::string & fieldName)
  {
    const auto cmp = [](const double & a, const float &b) -> bool { return util::AreEqual(a, b); };

    const auto assgn = [](const double & a, float & b) {  b = a; };
    return ValidateOrCopy<double, float>(input, target, unsetVal, cmp, assgn, fieldName);
  }


// ------------------------------------------------------------
  TruthMatcher::TruthMatcher(const std::vector<std::string> & ghepFilenames,
                             const genie::NtpMCEventRecord *gEvt,
                             std::function<int(const genie::NtpMCEventRecord *)> genieFillerCallback)
    : cafmaker::Loggable("TruthMatcher"),
      fGTrees(ghepFilenames, gEvt),
      fGENIEWriterCallback(std::move(genieFillerCallback))
  {}

  // --------------------------------------------------------------
  void TruthMatcher::FillInteraction(caf::SRTrueInteraction& nu, const genie::NtpMCEventRecord * gEvt)
  {

    genie::EventRecord * event = gEvt->event;
    genie::Interaction * in = event->Summary();

    // note that the nu.id and nu.genieIdx are filled in the calling function,
    // because that info is not stored inside the GENIE event proper

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

    TLorentzVector lepP4; //defaults to 0
    genie::GHepParticle* finallepton = event->FinalStatePrimaryLepton();
    if(finallepton)
      lepP4 = *(finallepton->P4());
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
    int stablePartIdx=0;
    for (int j=0; j< event->GetEntries(); j++)
    {
      auto p = dynamic_cast<const genie::GHepParticle *>((*event)[j]);

      if( p->Status() != genie::EGHepStatus::kIStStableFinalState
          && p->Status() != genie::EGHepStatus::kIStHadronInTheNucleus) continue;

      caf::SRTrueParticle part;
      part.G4ID = (p->Status() == genie::EGHepStatus::kIStStableFinalState) ? stablePartIdx++ : -1;
      part.pdg = p->Pdg();
      part.interaction_id = nu.id;
      part.p = *p->P4();
      // todo: original GENIE events haven't been transferred to det coordinates
      //       or bunched in spill time structure yet.  just leave out?
//      part.time = nu.time;
//      part.start_pos = p->X4()->Vect();

      // remaining fields need to be filled in with post-G4 info

      std::string process;
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

        process = "PRIMARY";
      }
      else // kIStHadronInTheNucleus
      {
        nu.prefsi.push_back(std::move(part));
        nu.nprefsi++;

        process = "PRE-FSI HADRON";
      }
      LOG_S("TruthMatcher::FillInteraction").DEBUG() << "  " << process << " particle "
                                                     << " (idx in GENIE vec = " << j << ", trk id = " << part.G4ID << ", pdg = " << p->Pdg() << ", energy = " << p->E() << ")"
                                                     << " has GENIE status " << p->Status() << "\n";

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
    LOG.VERBOSE() << "  Searching for true particle within interaction ID = " << ixn.id << "\n";

    caf::SRTrueParticle * part = nullptr;
    std::vector<caf::SRTrueParticle> & collection = (isPrimary) ? ixn.prim : ixn.sec;
    int & counter = (isPrimary) ? ixn.nprim : ixn.nsec;
    if ( auto itPart = std::find_if(collection.begin(), collection.end(), cmp);
         itPart == collection.end() )
    {
      if (!createNew)
        throw std::runtime_error("True particle from interaction ID " + std::to_string(ixn.id)
                                 + " was not found in the " + std::string(isPrimary ? "primary" : "secondary") + " true particle collection");
      else
        LOG.VERBOSE() << "  made a new SRTrueParticle in " << (isPrimary ? "prim" : "sec") << " collection \n";

      collection.emplace_back();
      counter++;

      part = &collection.back();
      part->interaction_id = ixn.id;
    }
    else
    {
      LOG.VERBOSE() << "    --> found previously created SRTrueParticle (interaction id = " << itPart->interaction_id << ", trk id = " << itPart->G4ID <<  ") .  Returning that.\n";
      part = &(*itPart);
    }

    return *part;
  }

  // ------------------------------------------------------------
  caf::SRTrueInteraction & TruthMatcher::GetTrueInteraction(caf::StandardRecord &sr, unsigned long ixnID, bool createNew) const
  {

    caf::SRTrueInteraction * ixn = nullptr;

    LOG.VERBOSE() << "   Searching for true interaction with interaction ID = " << ixnID << " (allowed to create new one: " << createNew << ")\n";

    // if we can't find a SRTrueInteraction with matching ID, we may need to make a new one
    if ( auto itIxn = std::find_if(sr.mc.nu.begin(), sr.mc.nu.end(),
                                   [ixnID](const caf::SRTrueInteraction & ixn) { return static_cast<unsigned long>(ixn.id) == ixnID; });
         itIxn == sr.mc.nu.end() )
    {
      if (!createNew)
      {
        LOG.FATAL() << "Could not find a true interaction with ID " << ixnID << " already in this StandardRecord, and you asked me not to make a new one\n";
        throw std::runtime_error("True interaction with interaction ID = " + std::to_string(ixnID) + " not found in this StandardRecord");
      }

      LOG.VERBOSE() << "    creating new SRTrueInteraction.  Trying to match to a GENIE event...\n";

      // todo: should this logic live somewhere else?
      unsigned int evtNum = ixnID % 1000000;
      unsigned long runNum = (ixnID - evtNum) / 1000000;
      if (HaveGENIE())
      {
        try
        {
          fGTrees.SelectEvent(runNum, evtNum);
        }
        catch (std::out_of_range & exc)
        {
          // intercept briefly to add a log message
          LOG.FATAL() << "Could not find GENIE tree for interaction with run number: " << runNum << "!  Abort.\n";
          throw exc;
        }
      }

      sr.mc.nu.emplace_back();
      sr.mc.nnu++;
      ixn = &sr.mc.nu.back();
      ixn->id = ixnID;
      if (HaveGENIE())
      {
        LOG.VERBOSE() << "      --> GENIE record found (" << fGTrees.GEvt() << ").  copying...\n";

        // this bit of info can't be extracted directly from the GENIE record,
        // so we do it here
        ixn->genieIdx = fGENIEWriterCallback(fGTrees.GEvt());  // copy the GENIE event into the CAF output GENIE tree

        FillInteraction(*ixn, fGTrees.GEvt());  // copy values from the GENIE event into the StandardRecord
      }
      else
        LOG.VERBOSE() << "      --> no matching GENIE interaction found.  New empty SRTrueInteraction will be returned.\n";
    } // if ( didn't find a matching SRTrueInteraction )
    else
    {
      LOG.VERBOSE() << "   Found previously created SRTrueInteraction.  Returning that.\n";
      ixn = &(*itIxn);
    }
    return *ixn;
  }

  // ------------------------------------------------------------
  bool TruthMatcher::HaveGENIE() const
  {
    static auto isNull = [](const std::pair<unsigned long int, const TTree*>& pair) -> bool { return !pair.second; };
    return !std::all_of(fGTrees.begin(), fGTrees.end(), isNull);
  }

  // ------------------------------------------------------------
  void TruthMatcher::SetLogThrehsold(cafmaker::Logger::THRESHOLD thresh)
  {
    Loggable::SetLogThrehsold(thresh);
    fGTrees.SetLogThrehsold(thresh);
  }

  // ------------------------------------------------------------
  TruthMatcher::GTreeContainer::GTreeContainer(const vector<std::string> &filenames, const genie::NtpMCEventRecord * gEvt)
    : cafmaker::Loggable("GTreeContainer"), fGEvt(gEvt)
  {
    for (const auto & fname : filenames)
    {
      if (fname.empty())
        continue;

      TTree * tree = nullptr;
      genie::NtpMCTreeHeader * hdr = nullptr;
      unsigned long run = 0;
      auto & f = fGFiles.emplace_back(std::unique_ptr<TFile>(TFile::Open(fname.c_str())));  // the file goes into the list here so the TTree we pull out of it never disappears
      if (f && !f->IsZombie())
      {
        tree = dynamic_cast<TTree *>(f->Get("gtree"));
        hdr = dynamic_cast<genie::NtpMCTreeHeader *>(f->Get("header"));
      }
      if (!tree || !hdr)
      {
        LOG.FATAL() << "Could not load TTree 'gtree' or associated header from supplied .ghep file: " << fname
                                          << "\n";
        abort();
      }

      run = hdr->runnu;
      if (run == 0)
      {
        // workaround for GENIE files where run number wasn't set
        std::regex pattern("(rock|nu)\\.(\\d+)");
        std::smatch matches;
        std::regex_search(fname, matches, pattern);
        if (matches.size() == 3)
          // this pattern from https://github.com/DUNE/2x2_sim/wiki/Production-changes-and-validation-findings#file-format-differences
          run = ((matches[1] == "rock" ? static_cast<int>(1e9) : 0) + std::stoull(matches[2]));
        else
        {
          LOG.ERROR() << "Got " << matches.size() << " pattern matches from this filename (expected: 2):\n";

          std::stringstream msg;
          msg << "   ";
          for (const auto & match : matches)
            msg << match << "\n   ";
          LOG.ERROR() << msg.str();

          msg.str("");
          msg << "Couldn't determine run number for events in GENIE file: '" << fname << "'\n";
          LOG.FATAL() << msg.str();
          throw std::runtime_error(msg.str());
        }
      }

      fGTrees[run] = tree;
      tree->SetBranchAddress("gmcrec", &fGEvt);

      LOG.INFO() << "Loaded TTree for run " << run << " from file: " << fname << "\n";
    }
  }

  // ------------------------------------------------------------
  const genie::NtpMCEventRecord * TruthMatcher::GTreeContainer::GEvt() const
  {
    return fGEvt;
  }

  // ------------------------------------------------------------
  void TruthMatcher::GTreeContainer::SelectEvent(unsigned long runNum, unsigned int evtNum)
  {
    auto it_tree = fGTrees.find(runNum);
    if (it_tree == fGTrees.end())
    {
      std::stringstream ss;
      ss << "Run number " << runNum << " was not found in this collection of .ghep files\n";
      LOG.FATAL() << ss.str();
      throw std::range_error(ss.str());
    }

    it_tree->second->GetEntry(evtNum);
  }

  // ------------------------------------------------------------
  void TruthMatcher::GTreeContainer::SetGEvtAddr(const genie::NtpMCEventRecord *evt)
  {
    fGEvt = evt;
  }
}

