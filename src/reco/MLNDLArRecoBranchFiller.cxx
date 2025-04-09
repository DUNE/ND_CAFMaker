
#include "MLNDLArRecoBranchFiller.h"

#include <limits>

// StandardRecord format
#include "duneanaobj/StandardRecord/StandardRecord.h"

// our headers
#include "DLP_h5_classes.h"
#include "Params.h"
#include "truth/FillTruth.h"

using namespace cafmaker::types::dlp;

// these overloads needed to make the ValidateOrCopy() templates functional over these types
std::ostream & operator<<(std::ostream& stream, InteractionMode mode)
{
  return stream << static_cast<std::underlying_type<InteractionMode>::type>(mode);
}

std::ostream & operator<<(std::ostream& stream, CurrentType curr)
{
  return stream << static_cast<std::underlying_type<CurrentType>::type>(curr);
}

namespace cafmaker
{
  caf::ScatteringMode DLP2CAF(cafmaker::types::dlp::InteractionMode mode)
  {
    using cafmaker::types::dlp::InteractionMode;

    switch(mode)
    {
      case InteractionMode::kQE:
        return caf::kQE;

      case InteractionMode::kDIS:
        return caf::kDIS;

      // for whatever reason they don't appear to have a catchall RES enumerator,
      // just all the various NUISANCE codes ... :-|
      case InteractionMode::kResCCNuBarDelta0PiMinus:        [[fallthrough]];
      case InteractionMode::kResCCNuBarDeltaMinusPiPlus:     [[fallthrough]];
      case InteractionMode::kResCCNuBarKaon0Lambda0:         [[fallthrough]];
      case InteractionMode::kResCCNuBarNeutronEta:           [[fallthrough]];
      case InteractionMode::kResCCNuBarNeutronPi0Pi0:        [[fallthrough]];
      case InteractionMode::kResCCNuBarNeutronPiMinus:       [[fallthrough]];
      case InteractionMode::kResCCNuBarNeutronPiPlusPiMinus: [[fallthrough]];
      case InteractionMode::kResCCNuBarNeutronRho0:          [[fallthrough]];
      case InteractionMode::kResCCNuBarNeutronRhoMinus:      [[fallthrough]];
      case InteractionMode::kResCCNuBarProtonPi0:            [[fallthrough]];
      case InteractionMode::kResCCNuBarProtonPi0Pi0:         [[fallthrough]];
      case InteractionMode::kResCCNuBarProtonPiMinus:        [[fallthrough]];
      case InteractionMode::kResCCNuBarSigma0Kaon0:          [[fallthrough]];
      case InteractionMode::kResCCNuBarSigmaMinusKaon0:      [[fallthrough]];
      case InteractionMode::kResCCNuDelta2PlusPiMinus:       [[fallthrough]];
      case InteractionMode::kResCCNuDeltaPlusPiPlus:         [[fallthrough]];
      case InteractionMode::kResCCNuKaonPlusLambda0:         [[fallthrough]];
      case InteractionMode::kResCCNuNeutronPi0:              [[fallthrough]];
      case InteractionMode::kResCCNuNeutronPiPlus:           [[fallthrough]];
      case InteractionMode::kResCCNuNeutronRhoPlus:          [[fallthrough]];
      case InteractionMode::kResCCNuProtonEta:               [[fallthrough]];
      case InteractionMode::kResCCNuProtonPi0Pi0:            [[fallthrough]];
      case InteractionMode::kResCCNuProtonPiPlus:            [[fallthrough]];
      case InteractionMode::kResCCNuProtonPiPlusPiMinus:     [[fallthrough]];
      case InteractionMode::kResCCNuProtonRhoPlus:           [[fallthrough]];
      case InteractionMode::kResCCNuSigmaPlusKaon0:          [[fallthrough]];
      case InteractionMode::kResCCNuSigmaPlusKaonPlus:       [[fallthrough]];
      case InteractionMode::kResNCNuBarNeutronPi0:           [[fallthrough]];
      case InteractionMode::kResNCNuBarNeutronPiMinus:       [[fallthrough]];
      case InteractionMode::kResNCNuBarProtonPi0:            [[fallthrough]];
      case InteractionMode::kResNCNuBarProtonPiPlus:         [[fallthrough]];
      case InteractionMode::kResNCNuNeutronPi0:              [[fallthrough]];
      case InteractionMode::kResNCNuNeutronPiMinus:          [[fallthrough]];
      case InteractionMode::kResNCNuProtonPi0:               [[fallthrough]];
      case InteractionMode::kResNCNuProtonPiPlus:
        return caf::kRes;

      case InteractionMode::kCoh:
        return caf::kCoh;

      case InteractionMode::kDiffractive:
        return caf::kDiffractive;

      case InteractionMode::kNuElectronElastic:
        return caf::kNuElectronElastic;

      case InteractionMode::kInverseMuDecay:
        return caf::kInvMuonDecay;

      case InteractionMode::kAMNuGamma:
        return caf::kAMNuGamma;

      case InteractionMode::kMEC:
        return caf::kMEC;

      case InteractionMode::kCohElastic:
        return caf::kCohElastic;

      case InteractionMode::kInverseBetaDecay:
        return caf::kInverseBetaDecay;

      case InteractionMode::kGlashowResonance:
        return caf::kGlashowResonance;

      case InteractionMode::kIMDAnnihilation:
        return caf::kIMDAnnihilation;

      case InteractionMode::kUnknownInteraction:
        return caf::kUnknownMode;

      default:
        std::cerr << "Unrecognized scattering mode: " << static_cast<int>(mode) << "\n";
        abort();
    }

  }


  // ------------------------------------------------------------------------------
  // todo: possibly build some mechanism for customizing the dataset names in the file here
  MLNDLArRecoBranchFiller::MLNDLArRecoBranchFiller(const std::string &h5filename)
    : IRecoBranchFiller("LArML"),
      fDSReader(h5filename,
                {{std::type_index(typeid(Particle)),                      "reco_particles"},
                 {std::type_index(typeid(Interaction)),                   "reco_interactions"},
                 {std::type_index(typeid(TrueParticle)),                  "truth_particles"},
                 {std::type_index(typeid(TrueInteraction)),               "truth_interactions"},
                 {std::type_index(typeid(Flash)),                         "flashes"},
                 {std::type_index(typeid(Event)),                         "events"},
                 {std::type_index(typeid(RunInfo)),                       "run_info"},
                 {std::type_index(typeid(cafmaker::types::dlp::Trigger)), "trigger"}}),  // needs to be disambiguated from CAFMaker's internal Trigger
      fTriggers(),
      fLastTriggerReqd(fTriggers.end())
  {
    // if we got this far, nothing bad happened trying to open the file or dataset
    SetConfigured(true);
  }

  // ------------------------------------------------------------------------------
  void
  MLNDLArRecoBranchFiller::_FillRecoBranches(const Trigger &trigger,
                                             caf::StandardRecord &sr,
                                             const cafmaker::Params &par,
                                             const TruthMatcher *truthMatcher) const

  {
    // figure out where in our list of triggers this event index is.
    // we should always be looking forwards, since we expect to be traversing in that direction
    auto it_start = (fLastTriggerReqd == fTriggers.end()) ? fTriggers.cbegin() : fLastTriggerReqd;
    auto itTrig = std::find(it_start, fTriggers.cend(), trigger);
    if (itTrig == fTriggers.end())
    {
      LOG.FATAL() << "Reco branch filler '" << GetName() << "' could not find trigger with evtID == " << trigger.evtID << "!  Abort.\n";
      abort();
    }
    std::size_t idx = std::distance(fTriggers.cbegin(), itTrig);

    LOG.VERBOSE() << "    Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID << ", internal evt idx = " << idx << ".\n";

    //Fill ND-LAr specific info in the meta branch
    H5DataView<cafmaker::types::dlp::RunInfo> run_info = fDSReader.GetProducts<cafmaker::types::dlp::RunInfo>(idx);
    sr.meta.lar2x2.enabled = true;
    for (const auto & runinf : run_info)
    {
      sr.meta.lar2x2.run = runinf.run;
      sr.meta.lar2x2.subrun = runinf.subrun;
      sr.meta.lar2x2.event = runinf.event;
    }
    
    sr.meta.lar2x2.readoutstart_s = trigger.triggerTime_s;
    sr.meta.lar2x2.readoutstart_ns = trigger.triggerTime_ns;

    H5DataView<cafmaker::types::dlp::Interaction> interactions = fDSReader.GetProducts<cafmaker::types::dlp::Interaction>(idx);
    H5DataView<cafmaker::types::dlp::TrueInteraction> trueInteractions = fDSReader.GetProducts<cafmaker::types::dlp::TrueInteraction>(idx);
    H5DataView<cafmaker::types::dlp::TrueParticle> trueParticles = fDSReader.GetProducts<cafmaker::types::dlp::TrueParticle>(idx);
    FillInteractions(interactions, trueInteractions, trueParticles, truthMatcher, sr);

    H5DataView<cafmaker::types::dlp::Particle> particles = fDSReader.GetProducts<cafmaker::types::dlp::Particle>(idx);
    FillParticles(particles, trueInteractions, trueParticles, truthMatcher, sr);

    FillTracks(particles, trueInteractions, trueParticles, truthMatcher, sr);
    FillShowers(particles, trueInteractions, trueParticles, truthMatcher, sr);
    H5DataView<cafmaker::types::dlp::Flash> flashes = fDSReader.GetProducts<cafmaker::types::dlp::Flash>(idx);
    FillFlashes(flashes, sr);

    // todo: now do some sanity checks:
    //       - compare the number of true particles in each dlp::TrueInteraction to the number discovered and filled in SRTrueInteraction
    //       - do the same with the reco particles
    //       - etc.


    // todo: figure out what to do with these
//    int64_t num_particles;
//    int64_t num_primaries;
//    std::array<int64_t, 6> particle_counts;
//    BufferView<int64_t> particle_ids;
//    std::array<int64_t, 6> primary_counts;
//    int64_t size;
//    char * topology;
//    BufferView<int64_t> truth_particle_counts;
//    BufferView<int64_t> truth_primary_counts;
//    char * truth_topology;
//    BufferView<double> truth_vertex;
//    char * units;
//    int64_t volume_id;

  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTrueInteraction(caf::SRTrueInteraction & srTrueInt,
                                                    const cafmaker::types::dlp::TrueInteraction & ptTrueInt /* pt = "pass-through" */) const
  {
    LOG.DEBUG() << "    now copying truth info from MLReco TrueInteraction to SRTrueInteraction...\n";

    const auto NaN = std::numeric_limits<float>::signaling_NaN();

    // vertices from ML-reco are adjusted to the edge of the sensitive detector volume
    // if they originate from outside it, so we can't use them
    //ValidateOrCopy(ptTrueInt.truth_vertex[0], srTrueInt.vtx.x, NaN, "SRTrueInteraction::vtx::x");
    //ValidateOrCopy(ptTrueInt.truth_vertex[1], srTrueInt.vtx.y, NaN, "SRTrueInteraction::vtx::y");
    //ValidateOrCopy(ptTrueInt.truth_vertex[2], srTrueInt.vtx.z, NaN, "SRTrueInteraction::vtx::z");

    const std::function<bool(const CurrentType &, const bool &)> nuCurrComp =
    [](const CurrentType & inCurr, const bool & outCurr)
    {
      return (outCurr && inCurr == cafmaker::types::dlp::CurrentType::kCC)
             || (!outCurr && inCurr == cafmaker::types::dlp::CurrentType::kNC);
    };
    const std::function<void(const CurrentType & inCurr, bool & outCurr)> nuCurrSet =
    [](const CurrentType & inCurr, bool & outCurr)
    {
      outCurr = inCurr == cafmaker::types::dlp::CurrentType::kCC;
    };

    // todo: these need us to propagate nu info through Supera.  WIP
    //ValidateOrCopy(ptTrueInt.nu_current_type, srTrueInt.iscc, false, nuCurrComp, nuCurrSet); //this is currently filled with -1 for iscc
//    ValidateOrCopy(ptTrueInt.nu_energy_init/1000., srTrueInt.E, NaN); //this is currently filled with many -inf
//    ValidateOrCopy(ptTrueInt.nu_interaction_mode, srTrueInt.mode, caf::ScatteringMode::kUnknownMode,
//                   [](const InteractionMode & inCurr, const caf::ScatteringMode & outCurr)
//                   {
//                     return DLP2CAF(inCurr) == outCurr;
//                   },
//                   [](const InteractionMode & inCurr, caf::ScatteringMode & outCurr)
//                   {
//                     outCurr = DLP2CAF(inCurr);
//                   });
    // InteractionType nu_interaction_type;    // this appears to be identical to nu_interaction_mode

    // int64_t image_id;      // ID of event passed to reco within the file.  use the event ID instead.
    // bool is_contained;     // If the whole event is contained.  we don't have a landing spot for this right now
    // bool is_neutrino;      // We really want the initiating PDG instead :-/
    // bool is_principal_match;          // for now at least we're going to focus on matching from the Reco end first
    // BufferView<int64_t> match_ids;        //   |
    // BufferView<float> match_overlaps;  //   |
    // uint8_t matched;                  //   v

    // int64_t nu_id;        // this is the index within the overlaid spill.  not really any more useful than just `id`

    // todo: figure out what to do with these
//    int64_t num_particles;
//    int64_t num_primaries;
//    std::array<int64_t, 6> particle_counts;
//    BufferView<int64_t> particle_ids;
//    std::array<int64_t, 6> primary_counts;
//    int64_t size;
//    char * topology;
//    BufferView<int64_t> truth_particle_counts;
//    BufferView<int64_t> truth_primary_counts;
//    char * truth_topology;
//    BufferView<double> truth_vertex;
//    char * units;
//    int64_t volume_id;

  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTrueParticle(caf::SRTrueParticle & srTruePart,
                                                 const cafmaker::types::dlp::TrueParticle & truePartPassthrough,
                                                 const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles) const
  {
    const auto NaN = std::numeric_limits<float>::signaling_NaN();
    ValidateOrCopy(truePartPassthrough.pdg_code, srTruePart.pdg, 0, "pdg_code");
    ValidateOrCopy(truePartPassthrough.track_id, srTruePart.G4ID, -1,"SRTrueParticle::track_id");

    ValidateOrCopy(truePartPassthrough.orig_interaction_id, srTruePart.interaction_id, -1L, "SRTrueParticle::interaction_id");

    const auto ancestorTypeComp = [](const char* inProc, const caf::TrueParticleID::PartType & outType)
    {
      // fixme: the process codes don't look like this
      if (strcmp(inProc, "primary") == 0)
        return outType == caf::TrueParticleID::kPrimary;
      else
        return outType == caf::TrueParticleID::kSecondary;
    };
    const auto ancestorTypeAssgn = [](const char* inProc, caf::TrueParticleID::PartType & outType)
    {
      // fixme: the process codes don't look like this
      if (strcmp(inProc, "primary") == 0)
        outType = caf::TrueParticleID::kPrimary;
      else
        outType = caf::TrueParticleID::kSecondary;
    };
    
    // todo: need to figure out how to translate "1::91" etc. to the enums...
//    ValidateOrCopy(truePartPassthrough.creation_process, srTruePart.start_process)
     ValidateOrCopy(truePartPassthrough.position[0], srTruePart.start_pos.x, NaN, "SRTrueParticle::start_pos.x");
     ValidateOrCopy(truePartPassthrough.position[1], srTruePart.start_pos.y, NaN, "SRTrueParticle::start_pos.y");
     ValidateOrCopy(truePartPassthrough.position[2], srTruePart.start_pos.z, NaN, "SRTrueParticle::start_pos.z");

     ValidateOrCopy(truePartPassthrough.end_position[0], srTruePart.end_pos.x, NaN, "SRTrueParticle::end_pos.x");
     ValidateOrCopy(truePartPassthrough.end_position[1], srTruePart.end_pos.y, NaN, "SRTrueParticle::end_pos.y");
     ValidateOrCopy(truePartPassthrough.end_position[2], srTruePart.end_pos.z, NaN, "SRTrueParticle::end_pos.z");

    ValidateOrCopy(truePartPassthrough.momentum[0]/1000., srTruePart.p.px, NaN, "SRTrueParticle::p.px");
    ValidateOrCopy(truePartPassthrough.momentum[1]/1000., srTruePart.p.py, NaN, "SRTrueParticle::p.py");
    ValidateOrCopy(truePartPassthrough.momentum[2]/1000., srTruePart.p.pz, NaN, "SRTrueParticle::p.pz");

    try
    {
      ValidateOrCopy(truePartPassthrough.energy_init / 1000., srTruePart.p.E, NaN, "SRTrueParticle::p.E");
    }
    catch (std::runtime_error & e)
    {
      auto diff = (truePartPassthrough.energy_init / 1000. - srTruePart.p.E);
      if (diff < 1) // < 1 MeV
      {
        LOG.WARNING() << "True particle energy (track id=" << srTruePart.G4ID << ", pdg=" << srTruePart.pdg << ", stored E=" << srTruePart.p.E << ")"
                      << " differs by " << diff << " MeV between stored (GENIE?) and ML-reco pass-through values";
      }
      else
        throw e;
    }


  }

  namespace
  {
    struct DLPIxnComp
    {
      long int ixnID;
      bool operator()(const cafmaker::types::dlp::TrueInteraction & ixn)
      {
        return ixn.id == ixnID;
      }
    };

    struct SRPartCmp
    {
      int trkid;
      bool operator()(const caf::SRTrueParticle & part) const
      {
        LOG_S("SRPartCmp").VERBOSE() << "       SRPartCmp::operator()():  looking for trk ID = " << trkid << ", this particle trkID = " << part.G4ID << "\n";
        return trkid == part.G4ID;
      }
    };
  }


  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillInteractions(const H5DataView<cafmaker::types::dlp::Interaction> &ixns,
                                                 const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueIxns,
                                                 const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                                                 const TruthMatcher * truthMatch,
                                                 caf::StandardRecord &sr) const
  {
    sr.common.ixn.dlp.reserve(ixns.size());
    sr.common.ixn.ndlp = ixns.size();

    sr.nd.lar.dlp.resize(ixns.size());
    sr.nd.lar.ndlp = ixns.size();
    
    LOG.DEBUG() << "Filling reco interactions...\n";
    int ixnidx = 0;
    for (const auto & ixn : ixns)
    {
      caf::SRInteraction interaction;
      interaction.id  = ixn.id;
      interaction.vtx  = caf::SRVector3D(ixn.vertex[0], ixn.vertex[1], ixn.vertex[2]);  // note: this branch suffers from "too many nested vectors" problem.  won't see vals in TBrowser unless using a FlatCAF
      LOG.VERBOSE() << " --> interaction id = "  << interaction.id << "\n";

      // if we *have* truth matches, we need to connect them now
      if (ixn.match_ids.size())
      {
        LOG.VERBOSE() << "  There are " << ixn.match_ids.size() << " matched true interactions:\n";
        for (std::size_t idx = 0; idx < ixn.match_ids.size(); idx++)
        {
          LOG.VERBOSE() << "  ** Match index " << idx << " --> truth ID " << ixn.match_ids[idx] << "\n";
          // here we need to search through the truth interactions and find the one with this ID (since it's no longer an index)
          static DLPIxnComp ixnCmp;
          ixnCmp.ixnID = ixn.match_ids[idx];
          auto itIxn = std::find_if(trueIxns.begin(), trueIxns.end(), ixnCmp);
          if (itIxn == trueIxns.end())
          {
            std::stringstream msg;
            msg << "Reco interaction claims to match to true interaction with ID " << ixnCmp.ixnID
                << ", but that interaction was not found in the list of true interactions\n";
            LOG.FATAL() << msg.str();
            throw std::out_of_range(msg.str());
          }
          cafmaker::types::dlp::TrueInteraction trueIxnPassThrough = *itIxn;

          LOG.VERBOSE() << "  Finding matched true interaction with ML-reco ID = " << trueIxnPassThrough.id
                        << " and interaction ID = " << trueIxnPassThrough.orig_id
                        << "\n";

          caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, trueIxnPassThrough.orig_id);

          LOG.VERBOSE() << "    --> resulting SRTrueInteraction has the following particles in it:\n";
          for (const caf::SRTrueParticle & part : srTrueInt.prim)
            LOG.VERBOSE() << "    (prim) id = " << part.G4ID << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";
          for (const caf::SRTrueParticle & part : srTrueInt.prefsi)
            LOG.VERBOSE() << "    (prefsi) id = " << part.G4ID << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";
          for (const caf::SRTrueParticle & part : srTrueInt.sec)
            LOG.VERBOSE() << "    (sec) id = " << part.G4ID  << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";

          // here we need to fill in any additional info
          // that GENIE didn't know about: e.g., secondary particles made by GEANT4
          FillTrueInteraction(srTrueInt, trueIxnPassThrough);

          // note that the interaction ID is GENIE's label for it, which may not be the same as the index in the vector
          std::size_t truthVecIdx = std::distance(sr.mc.nu.begin(),
                                                  std::find_if(sr.mc.nu.begin(),
                                                               sr.mc.nu.end(),
                                                               [&srTrueInt](const caf::SRTrueInteraction& nu)
                                                               {
                                                                 return nu.id == srTrueInt.id;
                                                               }));

          interaction.truth.push_back(truthVecIdx);
          interaction.truthOverlap.push_back(ixn.match_overlaps[idx]);

          LOG.VERBOSE() << "  ** end matched true interaction search for ML-reco ID " << trueIxnPassThrough.id << ".\n";
        }
      }

      sr.common.ixn.dlp.push_back(std::move(interaction));
      //Fill matched flash info
      caf::FlashMatch flashMatch;
      //flashMatch.id = ixn.flash_id;
      //flashMatch.time = ixn.flash_time;
      flashMatch.total_pe = ixn.flash_total_pe;
      flashMatch.hypothesis_pe = ixn.flash_hypo_pe;
      sr.nd.lar.dlp[ixnidx].flash.push_back(flashMatch);
      ixnidx++;
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillParticles(const H5DataView<cafmaker::types::dlp::Particle> &particles,
                                              const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueInxns,
                                              const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                                              const TruthMatcher * truthMatch,
                                              caf::StandardRecord &sr) const
  {
    LOG.DEBUG() << "Filling reco particles...\n";

    // note: used in the hack further below
    static SRPartCmp srPartCmp;

    //filling reco particles regardless of  type (track/shower)
    for (const auto & part : particles)
    {
      LOG.VERBOSE() << " --> reco particle id = "  << part.id << "\n";

      caf::SRRecoParticle reco_particle;
      if(part.is_primary) reco_particle.primary = true;
      reco_particle.start = caf::SRVector3D(part.start_point[0], part.start_point[1], part.start_point[2]);
      reco_particle.end = caf::SRVector3D(part.end_point[0], part.end_point[1], part.end_point[2]);
      reco_particle.contained = part.is_contained; // this is not just the vertex, but all energies are contained
      if(part.is_contained) reco_particle.tgtA = 40;
      reco_particle.pdg = part.pdg_code;
      reco_particle.p.x = part.momentum[0]/1000.;
      reco_particle.p.y = part.momentum[1]/1000.;
      reco_particle.p.z = part.momentum[2]/1000.;
      if(part.shape == types::dlp::Shape::kTrack)
      {
        if(part.is_contained)
        {
          reco_particle.E = part.csda_ke/1000.;
          reco_particle.E_method = caf::PartEMethod::kRange;
        }
        else
        {
      	  reco_particle.E = part.mcs_ke/1000.;
    	    reco_particle.E_method = caf::PartEMethod::kMCS;
        }
      }
      else
      {
        reco_particle.E = part.calo_ke/1000.;
        reco_particle.E_method = caf::PartEMethod::kCalorimetry;
      }

      if (part.match_ids.size())
      {
        for (std::size_t idx = 0; idx < part.match_ids.size(); idx++)
        {
          LOG.VERBOSE() << "   searching for matched true particle with ML reco index: " << part.match_ids[idx] << "\n";
          cafmaker::types::dlp::TrueParticle truePartPassThrough = trueParticles[part.match_ids[idx]];

          LOG.VERBOSE() << "      id = " << truePartPassThrough.id << "; "
                    << "track id = " << truePartPassThrough.track_id << "; "
                    << "interaction ID = " << truePartPassThrough.interaction_id << "; "
                    << "is primary = " << truePartPassThrough.is_primary << "; "
                    << "pdg = " << truePartPassThrough.pdg_code << "; "
                    << "energy = " << truePartPassThrough.energy_init
                    << "\n";

          // first ask for the right truth match from the matcher.
          // if we have GENIE info it'll come pre-filled with all its info & sub-particles
          static DLPIxnComp ixnCmp;
          ixnCmp.ixnID = truePartPassThrough.interaction_id;
          auto it_ixn = std::find_if(trueInxns.begin(), trueInxns.end(), ixnCmp);
          if (it_ixn == trueInxns.end())
          {
            std::stringstream ss;
            ss << "True particle ID " << truePartPassThrough.id << " claims to be associated with true interaction ID " << truePartPassThrough.interaction_id
               << " but no such interaction could be found!\n";
            LOG.FATAL() << ss.str();
            throw std::out_of_range(ss.str());
          }
          const cafmaker::types::dlp::TrueInteraction & trueIxn = *it_ixn;

          caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, trueIxn.orig_id, false);


          // we need this below because caf::TrueParticleID wants the *index* of the SRTrueInteraction
          int srTrueIntIdx = std::distance(sr.mc.nu.begin(),
                                           std::find_if(sr.mc.nu.begin(),
                                                        sr.mc.nu.end(),
                                                        [&srTrueInt](const caf::SRTrueInteraction& ixn) {return ixn.id == srTrueInt.id;}));

          bool is_primary = std::find_if(srTrueInt.prim.begin(), srTrueInt.prim.end(), 
                                   [&srTrueInt, &truePartPassThrough](const caf::SRTrueParticle& part) { return part.G4ID == truePartPassThrough.track_id; }) != srTrueInt.prim.end();
          srPartCmp.trkid = truePartPassThrough.track_id;
          caf::SRTrueParticle & srTruePart = is_primary ? truthMatch->GetTrueParticle(sr, srTrueInt, truePartPassThrough.track_id, srPartCmp, true, (!truthMatch->HaveGENIE()))
                                                        : truthMatch->GetTrueParticle(sr, srTrueInt, truePartPassThrough.track_id, srPartCmp, false, true);

          //  this will fill in any other fields that weren't copied from a GENIE record
          // (which also handles the case where this particle is a secondary)
          FillTrueParticle(srTruePart, truePartPassThrough, trueParticles);

          // the particle idx is within the GENIE vector, which may not be the same as the index in the vector here
          // first find the interaction that it goes with
          LOG.VERBOSE() << "      this particle is " << (is_primary ? "PRIMARY" : "SECONDARY") << "\n";
          std::vector<caf::SRTrueParticle> & collection = is_primary
                                                          ? srTrueInt.prim
                                                          : srTrueInt.sec;
          std::size_t truthVecIdx = std::distance(collection.begin(),
                                                  std::find_if(collection.begin(),
                                                               collection.end(),
                                                               srPartCmp));

          reco_particle.truth.push_back(caf::TrueParticleID{srTrueIntIdx,
                                                            is_primary ? caf::TrueParticleID::PartType::kPrimary
                                                                       :  caf::TrueParticleID::PartType::kSecondary,
                                                            static_cast<int>(truthVecIdx)});
          reco_particle.truthOverlap.push_back(part.match_overlaps[idx]);
        }
      }

      // note that interaction ID is not in general the same as the index within the sr.common.ixn.dlp vector
      // (some interaction IDs are filtered out as they're not beam triggers etc.)
      auto itIxn = std::find_if(sr.common.ixn.dlp.begin(), sr.common.ixn.dlp.end(),
                                [&part](const caf::SRInteraction & ixn){ return ixn.id == part.interaction_id; });
      if (itIxn == sr.common.ixn.dlp.end())
      {
        LOG.FATAL() << "Particle's interaction ID (" << part.interaction_id << ") does not match any in the DLP set!\n";
        abort();
      }
      sr.common.ixn.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].part.dlp.push_back(std::move(reco_particle));
      sr.common.ixn.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].part.ndlp++;

    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTracks(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                                           const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueInxns,
                                           const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                                           const TruthMatcher * truthMatch,
                                           caf::StandardRecord &sr) const
  {
    // note: used in the hack further below
    static SRPartCmp srPartCmp;

    for (const auto & part : particles)
    {
      // only choose 'particles' that correspond to Track type
      if (part.shape != types::dlp::Shape::kTrack)
        continue;


      caf::SRTrack track;
      // fill shower variables
      track.Evis = part.calo_ke/1000.;
      track.E = part.csda_ke/1000.; //range based energy
      track.start = caf::SRVector3D(part.start_point[0], part.start_point[1], part.start_point[2]);
      track.end = caf::SRVector3D(part.end_point[0], part.end_point[1], part.end_point[2]);
      track.dir = caf::SRVector3D(part.start_dir[0], part.start_dir[1], part.start_dir[2]);
      track.enddir = caf::SRVector3D(part.end_dir[0], part.end_dir[1], part.end_dir[2]);
      track.len_cm = sqrt(pow((part.start_point[0]-part.end_point[0]),2) + pow((part.start_point[1]-part.end_point[1]),2) + pow((part.start_point[2]-part.end_point[2]),2));
      if (part.match_ids.size())
      {
        for (std::size_t idx = 0; idx < part.match_ids.size(); idx++)
        {
          LOG.VERBOSE() << "   searching for matched true particle with ML reco index: " << part.match_ids[idx] << "\n";
          cafmaker::types::dlp::TrueParticle truePartPassThrough = trueParticles[part.match_ids[idx]];

          LOG.VERBOSE() << "      id = " << truePartPassThrough.id << "; "
                    << "track id = " << truePartPassThrough.track_id << "; "
                    << "interaction ID = " << truePartPassThrough.interaction_id << "; "
                    << "is primary = " << truePartPassThrough.is_primary << "; "
                    << "pdg = " << truePartPassThrough.pdg_code << "; "
                    << "energy = " << truePartPassThrough.energy_init
                    << "\n";

          // first ask for the right truth match from the matcher.
          // if we have GENIE info it'll come pre-filled with all its info & sub-particles
          static DLPIxnComp ixnCmp;
          ixnCmp.ixnID = truePartPassThrough.interaction_id;
          auto it_ixn = std::find_if(trueInxns.begin(), trueInxns.end(), ixnCmp);
          if (it_ixn == trueInxns.end())
          {
            std::stringstream ss;
            ss << "True particle ID " << truePartPassThrough.id << " claims to be associated with true interaction ID " << truePartPassThrough.interaction_id
               << " but no such interaction could be found!\n";
            LOG.FATAL() << ss.str();
            throw std::out_of_range(ss.str());
          }
          const cafmaker::types::dlp::TrueInteraction & trueIxn = *it_ixn;

          caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, trueIxn.orig_id, false);

          // we need this below because caf::TrueParticleID wants the *index* of the SRTrueInteraction
          int srTrueIntIdx = std::distance(sr.mc.nu.begin(),
                                           std::find_if(sr.mc.nu.begin(),
                                                        sr.mc.nu.end(),
                                                        [&srTrueInt](const caf::SRTrueInteraction& ixn) {return ixn.id == srTrueInt.id;}));
    	  
           
          bool is_primary = std::find_if(srTrueInt.prim.begin(), srTrueInt.prim.end(), 
                                   [&srTrueInt, &truePartPassThrough](const caf::SRTrueParticle& part) { return part.G4ID == truePartPassThrough.track_id; }) != srTrueInt.prim.end();
          srPartCmp.trkid = truePartPassThrough.track_id;

          // we want to make sure the particle is created, if it isn't there,
          // but we won't do anything further with it, so we throw the return value away
          if (is_primary)
            truthMatch->GetTrueParticle(sr, srTrueInt, truePartPassThrough.track_id, srPartCmp, true, (!truthMatch->HaveGENIE()));
          else
            truthMatch->GetTrueParticle(sr, srTrueInt, truePartPassThrough.track_id, srPartCmp, false, true);

          // the particle idx is within the GENIE vector, which may not be the same as the index in the vector here
          // first find the interaction that it goes with
          LOG.VERBOSE() << "      this particle is " << (is_primary ? "PRIMARY" : "SECONDARY") << "\n";
          std::vector<caf::SRTrueParticle> & collection = is_primary
                                                          ? srTrueInt.prim
                                                          : srTrueInt.sec;
          std::size_t truthVecIdx = std::distance(collection.begin(),
                                                  std::find_if(collection.begin(),
                                                               collection.end(),
                                                               srPartCmp));

          track.truth.push_back(caf::TrueParticleID{srTrueIntIdx,
                                                            is_primary ? caf::TrueParticleID::PartType::kPrimary
                                                                       :  caf::TrueParticleID::PartType::kSecondary,
                                                            static_cast<int>(truthVecIdx)});
          track.truthOverlap.push_back(part.match_overlaps[idx]);
        }
      }
      // note that interaction ID is not in general the same as the index within the sr.common.ixn.dlp vector
      // (some interaction IDs are filtered out as they're not beam triggers etc.)
      auto itIxn = std::find_if(sr.common.ixn.dlp.begin(), sr.common.ixn.dlp.end(),
                                [&part](const caf::SRInteraction & ixn){ return ixn.id == part.interaction_id; });
      if (itIxn == sr.common.ixn.dlp.end())
      {
        LOG.FATAL() << "Particle's interaction ID (" << part.interaction_id << ") does not match any in the DLP set!\n";
        abort();
      }
      sr.nd.lar.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].tracks.push_back(std::move(track));
      sr.nd.lar.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].ntracks++;
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillShowers(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                                            const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueInxns,
                                            const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                                            const TruthMatcher * truthMatch,
                                            caf::StandardRecord &sr) const
  {
    // note: used in the hack further below
    static SRPartCmp srPartCmp;

    for (const auto & part : particles)
    {
      if (part.shape != types::dlp::Shape::kShower)
        continue;

      caf::SRShower shower;
      // fill shower variables
      shower.Evis = part.calo_ke/1000.;
      shower.start = caf::SRVector3D(part.start_point[0], part.start_point[1], part.start_point[2]);
      shower.direction = caf::SRVector3D(part.start_dir[0], part.start_dir[1], part.start_dir[2]);
      if (part.match_ids.size())
      {
        for (std::size_t idx = 0; idx < part.match_ids.size(); idx++)
        {
          LOG.VERBOSE() << "   searching for matched true particle with ML reco index: " << part.match_ids[idx] << "\n";
          cafmaker::types::dlp::TrueParticle truePartPassThrough = trueParticles[part.match_ids[idx]];

          LOG.VERBOSE() << "      id = " << truePartPassThrough.id << "; "
                    << "track id = " << truePartPassThrough.track_id << "; "
                    << "interaction ID = " << truePartPassThrough.interaction_id << "; "
                    << "is primary = " << truePartPassThrough.is_primary << "; "
                    << "pdg = " << truePartPassThrough.pdg_code << "; "
                    << "energy = " << truePartPassThrough.energy_init
                    << "\n";

          // first ask for the right truth match from the matcher.
          // if we have GENIE info it'll come pre-filled with all its info & sub-particles
          static DLPIxnComp ixnCmp;
          ixnCmp.ixnID = truePartPassThrough.interaction_id;
          auto it_ixn = std::find_if(trueInxns.begin(), trueInxns.end(), ixnCmp);
          if (it_ixn == trueInxns.end())
          {
            std::stringstream ss;
            ss << "True particle ID " << truePartPassThrough.id << " claims to be associated with true interaction ID " << truePartPassThrough.interaction_id
               << " but no such interaction could be found!\n";
            LOG.FATAL() << ss.str();
            throw std::out_of_range(ss.str());
          }
          const cafmaker::types::dlp::TrueInteraction & trueIxn = *it_ixn;

          caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, trueIxn.orig_id, false);

          // we need this below because caf::TrueParticleID wants the *index* of the SRTrueInteraction
          int srTrueIntIdx = std::distance(sr.mc.nu.begin(),
                                           std::find_if(sr.mc.nu.begin(),
                                                        sr.mc.nu.end(),
                                                        [&srTrueInt](const caf::SRTrueInteraction& ixn) {return ixn.id == srTrueInt.id;}));

    	    bool is_primary = std::find_if(srTrueInt.prim.begin(), srTrueInt.prim.end(), 
                                   [&srTrueInt, &truePartPassThrough](const caf::SRTrueParticle& part) { return part.G4ID == truePartPassThrough.track_id; }) != srTrueInt.prim.end();
          srPartCmp.trkid = truePartPassThrough.track_id;
          // we don't actually need the return value here for anything,
          // but we do want the TruthMatcher to *create* a new particle when that's appropriate
          is_primary ? truthMatch->GetTrueParticle(sr, srTrueInt, truePartPassThrough.track_id, srPartCmp, true, (!truthMatch->HaveGENIE()))
                     : truthMatch->GetTrueParticle(sr, srTrueInt, truePartPassThrough.track_id, srPartCmp, false, true);


          // the particle idx is within the GENIE vector, which may not be the same as the index in the vector here
          // first find the interaction that it goes with
          LOG.VERBOSE() << "      this particle is " << (is_primary ? "PRIMARY" : "SECONDARY") << "\n";
          std::vector<caf::SRTrueParticle> & collection = is_primary
                                                          ? srTrueInt.prim
                                                          : srTrueInt.sec;
          std::size_t truthVecIdx = std::distance(collection.begin(),
                                                  std::find_if(collection.begin(),
                                                               collection.end(),
                                                               srPartCmp));

          shower.truth.push_back(caf::TrueParticleID{srTrueIntIdx,
                                                            is_primary ? caf::TrueParticleID::PartType::kPrimary
                                                                       :  caf::TrueParticleID::PartType::kSecondary,
                                                            static_cast<int>(truthVecIdx)});
          shower.truthOverlap.push_back(part.match_overlaps[idx]);
        }
      }
      // note that interaction ID is not in general the same as the index within the sr.common.ixn.dlp vector
      // (some interaction IDs are filtered out as they're not beam triggers etc.)
      auto itIxn = std::find_if(sr.common.ixn.dlp.begin(), sr.common.ixn.dlp.end(),
                                [&part](const caf::SRInteraction & ixn){ return ixn.id == part.interaction_id; });
      if (itIxn == sr.common.ixn.dlp.end())
      {
        LOG.FATAL() << "Particle's interaction ID (" << part.interaction_id << ") does not match any in the DLP set!\n";
        abort();
      }
      sr.nd.lar.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].showers.push_back(std::move(shower));
      sr.nd.lar.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].nshowers++;

    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillFlashes(const H5DataView<cafmaker::types::dlp::Flash> & flashes,
                                            caf::StandardRecord &sr) const
  {

    for (const auto & flash : flashes)
    {

      caf::SROpticalFlash opflash;
      // fill flash variables for all flashes

      opflash.id = flash.id;
      //opflash.tpc_id = flash.tpc; //TODO
      opflash.time = flash.time;
      opflash.time_width = flash.time_width;
      opflash.total_pe = flash.total_pe;

      sr.nd.lar.flashes.push_back(std::move(opflash));
      sr.nd.lar.nflashes++;

    }
    
  }
  // ------------------------------------------------------------------------------
  std::deque<Trigger> MLNDLArRecoBranchFiller::GetTriggers(int triggerType) const
  {
    if (fTriggers.empty())
    {
      auto triggersIn = fDSReader.GetProducts<cafmaker::types::dlp::Trigger>(-1); // get ALL the Trigger products
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "' from " << triggersIn.size() << " ND-LAr RunInfo products:\n";
      fTriggers.reserve(triggersIn.size());
      for (const cafmaker::types::dlp::Trigger &trigger: triggersIn)
      {
        const int placeholderTriggerType = 0;
        // fixme: this check needs to be fixed when we have trigger type info
        if (triggerType >= 0 && triggerType != placeholderTriggerType)
        {
          LOG.VERBOSE() << "    skipping trigger ID=" << trigger.id << "\n";
          continue;
        }

        fTriggers.emplace_back();
        Trigger & trig = fTriggers.back();
        trig.evtID = trigger.id;

        trig.triggerType = trigger.type;
        trig.triggerTime_s = trigger.time_s;
        trig.triggerTime_ns = trigger.time_ns;

        LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
                      << ", triggerType=" << trig.triggerType
                      << ", triggerTime_s=" << trig.triggerTime_s
                      << ", triggerTime_ns=" << trig.triggerTime_ns
                      << "\n";
      }
      fLastTriggerReqd = fTriggers.end();  // since we just modified the list, any iterators have been invalidated
    }

    std::deque<Trigger> triggers;
    for (const Trigger & trigger : fTriggers)
    {
      if (triggerType < 0 || triggerType == fTriggers.back().triggerType)
        triggers.push_back(trigger);
    }

    return triggers;
  }

} // namespace cafmaker
