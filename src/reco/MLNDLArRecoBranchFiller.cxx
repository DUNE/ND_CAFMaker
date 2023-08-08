#include "MLNDLArRecoBranchFiller.h"

#include <limits>

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "DLP_h5_classes.h"
#include "Params.h"
#include "truth/FillTruth.h"

using namespace cafmaker::types::dlp;

// these overloads needed to make the ValidateOrCopy() templates functional over these types
std::ostream & operator<<(std::ostream& stream, NuInteractionMode mode)
{
  return stream << static_cast<std::underlying_type<NuInteractionMode>::type>(mode);
}

std::ostream & operator<<(std::ostream& stream, NuCurrentType curr)
{
  return stream << static_cast<std::underlying_type<NuCurrentType>::type>(curr);
}

namespace cafmaker
{
  caf::ScatteringMode DLP2CAF(cafmaker::types::dlp::NuInteractionMode mode)
  {
    using cafmaker::types::dlp::NuInteractionMode;

    switch(mode)
    {
      case NuInteractionMode::kQE:
        return caf::kQE;

      case NuInteractionMode::kDIS:
        return caf::kDIS;

      // for whatever reason they don't appear to have a catchall RES enumerator,
      // just all the various NUISANCE codes ... :-|
      case NuInteractionMode::kResCCNuBarDelta0PiMinus:        [[fallthrough]];
      case NuInteractionMode::kResCCNuBarDeltaMinusPiPlus:     [[fallthrough]];
      case NuInteractionMode::kResCCNuBarKaon0Lambda0:         [[fallthrough]];
      case NuInteractionMode::kResCCNuBarNeutronEta:           [[fallthrough]];
      case NuInteractionMode::kResCCNuBarNeutronPi0Pi0:        [[fallthrough]];
      case NuInteractionMode::kResCCNuBarNeutronPiMinus:       [[fallthrough]];
      case NuInteractionMode::kResCCNuBarNeutronPiPlusPiMinus: [[fallthrough]];
      case NuInteractionMode::kResCCNuBarNeutronRho0:          [[fallthrough]];
      case NuInteractionMode::kResCCNuBarNeutronRhoMinus:      [[fallthrough]];
      case NuInteractionMode::kResCCNuBarProtonPi0:            [[fallthrough]];
      case NuInteractionMode::kResCCNuBarProtonPi0Pi0:         [[fallthrough]];
      case NuInteractionMode::kResCCNuBarProtonPiMinus:        [[fallthrough]];
      case NuInteractionMode::kResCCNuBarSigma0Kaon0:          [[fallthrough]];
      case NuInteractionMode::kResCCNuBarSigmaMinusKaon0:      [[fallthrough]];
      case NuInteractionMode::kResCCNuDelta2PlusPiMinus:       [[fallthrough]];
      case NuInteractionMode::kResCCNuDeltaPlusPiPlus:         [[fallthrough]];
      case NuInteractionMode::kResCCNuKaonPlusLambda0:         [[fallthrough]];
      case NuInteractionMode::kResCCNuNeutronPi0:              [[fallthrough]];
      case NuInteractionMode::kResCCNuNeutronPiPlus:           [[fallthrough]];
      case NuInteractionMode::kResCCNuNeutronRhoPlus:          [[fallthrough]];
      case NuInteractionMode::kResCCNuProtonEta:               [[fallthrough]];
      case NuInteractionMode::kResCCNuProtonPi0Pi0:            [[fallthrough]];
      case NuInteractionMode::kResCCNuProtonPiPlus:            [[fallthrough]];
      case NuInteractionMode::kResCCNuProtonPiPlusPiMinus:     [[fallthrough]];
      case NuInteractionMode::kResCCNuProtonRhoPlus:           [[fallthrough]];
      case NuInteractionMode::kResCCNuSigmaPlusKaon0:          [[fallthrough]];
      case NuInteractionMode::kResCCNuSigmaPlusKaonPlus:       [[fallthrough]];
      case NuInteractionMode::kResNCNuBarNeutronPi0:           [[fallthrough]];
      case NuInteractionMode::kResNCNuBarNeutronPiMinus:       [[fallthrough]];
      case NuInteractionMode::kResNCNuBarProtonPi0:            [[fallthrough]];
      case NuInteractionMode::kResNCNuBarProtonPiPlus:         [[fallthrough]];
      case NuInteractionMode::kResNCNuNeutronPi0:              [[fallthrough]];
      case NuInteractionMode::kResNCNuNeutronPiMinus:          [[fallthrough]];
      case NuInteractionMode::kResNCNuProtonPi0:               [[fallthrough]];
      case NuInteractionMode::kResNCNuProtonPiPlus:
        return caf::kRes;

      case NuInteractionMode::kCoh:
        return caf::kCoh;

      case NuInteractionMode::kDiffractive:
        return caf::kDiffractive;

      case NuInteractionMode::kNuElectronElastic:
        return caf::kNuElectronElastic;

      case NuInteractionMode::kInverseMuDecay:
        return caf::kInvMuonDecay;

      case NuInteractionMode::kAMNuGamma:
        return caf::kAMNuGamma;

      case NuInteractionMode::kMEC:
        return caf::kMEC;

      case NuInteractionMode::kCohElastic:
        return caf::kCohElastic;

      case NuInteractionMode::kInverseBetaDecay:
        return caf::kInverseBetaDecay;

      case NuInteractionMode::kGlashowResonance:
        return caf::kGlashowResonance;

      case NuInteractionMode::kIMDAnnihilation:
        return caf::kIMDAnnihilation;

      case NuInteractionMode::kUnknownInteraction:
        return caf::kUnknownMode;

      default:
        std::cerr << "Unrecognized scattering mode: " << static_cast<int>(mode) << "\n";
        abort();
    }

  }


  // ------------------------------------------------------------------------------
  // todo: possibly build some mechanism for customizing the dataset names in the file here
  MLNDLArRecoBranchFiller::MLNDLArRecoBranchFiller(const std::string &h5filename)
    : fDSReader(h5filename,
                {{std::type_index(typeid(Particle)),         "particles"},
                 {std::type_index(typeid(Interaction)),      "interactions"},
                 {std::type_index(typeid(TrueParticle)),     "truth_particles"},
                 {std::type_index(typeid(TrueInteraction)),  "truth_interactions"},
                 {std::type_index(typeid(Event)),            "events"},
                 {std::type_index(typeid(RunInfo)), "run_info"}}),
      fTriggers(),
      fLastTriggerReqd(fTriggers.end())
  {
    // if we got this far, nothing bad happened trying to open the file or dataset
    SetConfigured(true);
    name = "LArML";
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
      std::cerr << "Reco branch filler '" << GetName() << "' could not find trigger with evtID == " << trigger.evtID << "!  Abort.\n";
      abort();
    }
    std::size_t idx = std::distance(fTriggers.cbegin(), itTrig);


    std::cout << "    Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID << ", internal evt idx = " << idx << ".\n";
    H5DataView<cafmaker::types::dlp::Interaction> interactions = fDSReader.GetProducts<cafmaker::types::dlp::Interaction>(idx);
    H5DataView<cafmaker::types::dlp::TrueInteraction> trueInteractions = fDSReader.GetProducts<cafmaker::types::dlp::TrueInteraction>(idx);
    FillInteractions(interactions, trueInteractions, truthMatcher, sr);

    H5DataView<cafmaker::types::dlp::Particle> particles = fDSReader.GetProducts<cafmaker::types::dlp::Particle>(idx);
    H5DataView<cafmaker::types::dlp::TrueParticle> trueParticles = fDSReader.GetProducts<cafmaker::types::dlp::TrueParticle>(idx);
    FillParticles(particles, trueParticles, truthMatcher, sr);

    FillTracks(particles, sr);
    FillShowers(particles, sr);

    // todo: now do some sanity checks:
    //       - compare the number of true particles in each dlp::TrueInteraction to the number discovered and filled in SRTrueInteraction
    //       - do the same with the reco particles
    //       - etc.

    //Fill ND-LAr specificinfo in the meta branch
    sr.meta.nd_lar.enabled = true;
    sr.meta.nd_lar.run = par().runInfo().run();
    sr.meta.nd_lar.subrun = par().runInfo().subrun();
    sr.meta.nd_lar.event = evtIdx;
  
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTrueInteraction(caf::SRTrueInteraction & srTrueInt,
                                                    const cafmaker::types::dlp::TrueInteraction & ptTrueInt /* pt = "pass-through" */) const
  {
    const auto NaN = std::numeric_limits<float>::signaling_NaN();

    ValidateOrCopy(ptTrueInt.id, srTrueInt.id, -1);

    ValidateOrCopy(ptTrueInt.vertex[0], srTrueInt.vtx.x, NaN);
    ValidateOrCopy(ptTrueInt.vertex[1], srTrueInt.vtx.y, NaN);
    ValidateOrCopy(ptTrueInt.vertex[2], srTrueInt.vtx.z, NaN);

    const std::function<bool(const NuCurrentType &, const bool &)> nuCurrComp =
      [](const NuCurrentType & inCurr, const bool & outCurr)
      {
        return (outCurr && inCurr == cafmaker::types::dlp::NuCurrentType::kCC)
               || (!outCurr && inCurr == cafmaker::types::dlp::NuCurrentType::kNC);
      };
    const std::function<void(const NuCurrentType & inCurr, bool & outCurr)> nuCurrSet =
      [](const NuCurrentType & inCurr, bool & outCurr)
      {
        outCurr = inCurr == cafmaker::types::dlp::NuCurrentType::kCC;
      };

    ValidateOrCopy(ptTrueInt.nu_current_type, srTrueInt.iscc, false,
                   nuCurrComp, nuCurrSet);

    ValidateOrCopy(ptTrueInt.nu_energy_init, srTrueInt.E, NaN);

    // int64_t image_id;      // ID of event passed to reco within the file.  use the event ID instead.
    // bool is_contained;     // If the whole event is contained.  we don't have a landing spot for this right now
    // bool is_neutrino;      // We really want the initiating PDG instead :-\
    // bool is_principal_match;          // for now at least we're going to focus on matching from the Reco end first
    // BufferView<int64_t> match;        //   |
    // BufferView<float> match_overlap;  //   |
    // uint8_t matched;                  //   v

    // int64_t nu_id;        // this is the index within the overlaid spill.  not really any more useful than just `id`

    ValidateOrCopy(ptTrueInt.nu_interaction_mode, srTrueInt.mode, caf::ScatteringMode::kUnknownMode,
                   [](const NuInteractionMode & inCurr, const caf::ScatteringMode & outCurr)
                   {
                     return DLP2CAF(inCurr) == outCurr;
                   },
                   [](const NuInteractionMode & inCurr, caf::ScatteringMode & outCurr)
                   {
                     outCurr = DLP2CAF(inCurr);
                   });

    // NuInteractionType nu_interaction_type;    // this appears to be identical to nu_interaction_mode


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
                                                 const cafmaker::types::dlp::TrueParticle & truePartPassthrough) const
  {
    const auto NaN = std::numeric_limits<float>::signaling_NaN();

    ValidateOrCopy(truePartPassthrough.interaction_id, srTruePart.interaction_id, -1);
    ValidateOrCopy(truePartPassthrough.ancestor_track_id, srTruePart.ancestor_id.ixn, -1);

    const auto ancestorTypeComp = [](const char* inProc, const caf::TrueParticleID::PartType & outType)
                                  {
                                    if (strcmp(inProc, "primary") == 0)
                                      return outType == caf::TrueParticleID::kPrimary;
                                    else
                                      return outType == caf::TrueParticleID::kSecondary;
                                  };
    const auto ancestorTypeAssgn = [](const char* inProc, caf::TrueParticleID::PartType & outType)
                                   {
                                     if (strcmp(inProc, "primary") == 0)
                                       outType = caf::TrueParticleID::kPrimary;
                                     else
                                       outType = caf::TrueParticleID::kSecondary;
                                   };
    ValidateOrCopy(truePartPassthrough.ancestor_creation_process, srTruePart.ancestor_id.type, caf::TrueParticleID::kUnknown,
                   ancestorTypeComp, ancestorTypeAssgn);

    // todo: this is incorrect; the track_id (what we have) won't be the same as the index of the ancestor SRParticle (what we want).
    //       to fix this I think we need access to the SRTrueInteraction for this particle too, so we can dig around in its particle vectors
    ValidateOrCopy(truePartPassthrough.ancestor_track_id, srTruePart.ancestor_id.part, -1);

    ValidateOrCopy(truePartPassthrough.parent_track_id, srTruePart.parent, -1);

    // todo: need to figure out how to translate "1::91" etc. to the enums...
//    ValidateOrCopy(truePartPassthrough.creation_process, srTruePart.start_process)

    ValidateOrCopy(truePartPassthrough.start_point[0], srTruePart.start_pos.x, NaN);
    ValidateOrCopy(truePartPassthrough.start_point[1], srTruePart.start_pos.y, NaN);
    ValidateOrCopy(truePartPassthrough.start_point[2], srTruePart.start_pos.z, NaN);

    ValidateOrCopy(truePartPassthrough.end_point[0], srTruePart.end_pos.x, NaN);
    ValidateOrCopy(truePartPassthrough.end_point[1], srTruePart.end_pos.y, NaN);
    ValidateOrCopy(truePartPassthrough.end_point[2], srTruePart.end_pos.z, NaN);

    ValidateOrCopy(truePartPassthrough.energy_init, srTruePart.p.E, NaN);
    ValidateOrCopy(truePartPassthrough.momentum[0], srTruePart.p.px, NaN);
    ValidateOrCopy(truePartPassthrough.momentum[1], srTruePart.p.py, NaN);
    ValidateOrCopy(truePartPassthrough.momentum[2], srTruePart.p.pz, NaN);


  }
  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillInteractions(const H5DataView<cafmaker::types::dlp::Interaction> &Inxns,
                                                 const H5DataView<cafmaker::types::dlp::TrueInteraction> &trueInxns,
                                                 const TruthMatcher * truthMatch,
                                                 caf::StandardRecord &sr) const
  {
    sr.common.ixn.dlp.resize(Inxns.size());
    sr.common.ixn.ndlp = Inxns.size();

    //filling interactions
    for (const auto & inx : Inxns)
    {
      caf::SRInteraction interaction;
      interaction.id  = inx.id;
      interaction.vtx  = caf::SRVector3D(inx.vertex[0], inx.vertex[1], inx.vertex[2]);  // note: this branch suffers from "too many nested vectors" problem.  won't see vals in TBrowser unless using a FlatCAF

      // if we *have* truth matches, we need to connect them now
      if (inx.matched)
      {
        for (std::size_t idx = 0; idx < inx.match.size(); idx++)
        {
          cafmaker::types::dlp::TrueInteraction trueIxnPassThrough = trueInxns[inx.match[idx]];

          // first ask for the right truth match from the matcher.
          // if we have GENIE info it'll come pre-filled with all its info & sub-particles
          caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, trueIxnPassThrough.id);

          // here we need to fill in any additional info
          // that GENIE didn't know about: e.g., secondary particles made by GEANT4
          FillTrueInteraction(srTrueInt, trueIxnPassThrough);

          // note that the interaction ID is GENIE's label for it, which may not be the same as the index in the vector
          std::size_t truthVecIdx = std::distance(sr.mc.nu.begin(),
                                                  std::find_if(sr.mc.nu.begin(),
                                                               sr.mc.nu.end(),
                                                               [&trueIxnPassThrough](const caf::SRTrueInteraction& nu)
                                                               {
                                                                 return nu.id == trueIxnPassThrough.id;
                                                               }));

          interaction.truth.push_back(truthVecIdx);
          interaction.truthOverlap.push_back(inx.match_overlap[idx]);
        }
      }

      sr.common.ixn.dlp.push_back(std::move(interaction));
      sr.common.ixn.ndlp++;
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillParticles(const H5DataView<cafmaker::types::dlp::Particle> &particles,
                                              const H5DataView<cafmaker::types::dlp::TrueParticle> &trueParticles,
                                              const TruthMatcher * truthMatch,
                                              caf::StandardRecord &sr) const
  {
    //filling reco particles regardless of semantic type (track/shower)
    for (const auto & part : particles)
    {
      caf::SRRecoParticle reco_particle;
      if(part.is_primary) reco_particle.primary  = true;
      reco_particle.start = caf::SRVector3D(part.start_point[0], part.start_point[1], part.start_point[2]);
      reco_particle.end = caf::SRVector3D(part.end_point[0], part.end_point[1], part.end_point[2]); 
      reco_particle.E = part.depositions_sum;
      //To do: momentum mcs is currently filled with just -1
/*      reco_particle.p.x = part.momentum_mcs[0];
      reco_particle.p.y = part.momentum_mcs[1];
      reco_particle.p.z = part.momentum_mcs[2];
  */

      if (part.matched)
      {
        for (std::size_t idx = 0; idx < part.match.size(); idx++)
        {
          cafmaker::types::dlp::TrueParticle truePartPassThrough = trueParticles[part.match[idx]];

          // find the true particle this reco particle goes with.
          // if we had GENIE info and it was a primary, it should already be filled in.
          caf::SRTrueParticle & srTruePart = truthMatch->GetTrueParticle(sr,
                                                                         truePartPassThrough.interaction_id,
                                                                         truePartPassThrough.id,
                                                                         truePartPassThrough.is_primary);

          // however this will fill in any other fields that weren't copied from a GENIE record
          // (which also handles the case where this particle is a secondary)
          FillTrueParticle(srTruePart, truePartPassThrough);

          // the particle idx is within the GENIE vector, which may not be the same as the index in the vector here
          // first find the interaction that it goes with
          std::vector<caf::SRTrueParticle> & collection = (truePartPassThrough.is_primary)
                                                          ? sr.mc.nu[srTruePart.interaction_id].prim
                                                          : sr.mc.nu[srTruePart.interaction_id].sec;
          std::size_t truthVecIdx = std::distance(collection.begin(),
                                                  std::find_if(collection.begin(),
                                                               collection.end(),
                                                               [&truePartPassThrough](const caf::SRTrueParticle& part)
                                                               {
                                                                 return part.G4ID == truePartPassThrough.id;
                                                               }));

          reco_particle.truth.push_back(caf::TrueParticleID{srTruePart.interaction_id,
                                                            (truePartPassThrough.is_primary) ? caf::TrueParticleID::PartType::kPrimary :  caf::TrueParticleID::PartType::kSecondary,
                                                            static_cast<int>(truthVecIdx)});
          reco_particle.truthOverlap.push_back(truePartPassThrough.match_overlap[idx]);
        }
      }

      // note that interaction ID is not in general the same as the index within the sr.common.ixn.dlp vector
      // (some interaction IDs are filtered out as they're not beam triggers etc.)
      auto itIxn = std::find_if(sr.common.ixn.dlp.begin(), sr.common.ixn.dlp.end(),
                                [&part](const caf::SRInteraction & ixn){ return ixn.id == part.interaction_id; });
      if (itIxn == sr.common.ixn.dlp.end())
      {
        std::cerr << "ERROR: Particle's interaction ID (" << part.interaction_id << ") does not match any in the DLP set!\n";
        abort();
      }
      sr.common.ixn.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].part.dlp.push_back(std::move(reco_particle));
  
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillTracks(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                                           caf::StandardRecord &sr) const
  {
    sr.nd.lar.dlp.resize(sr.common.ixn.dlp.size());

    for (const auto & part : particles)
    {
      // only choose 'particles' that correspond to Track semantic type
      if (part.semantic_type != types::dlp::SemanticType::kTrack)
        continue;


      caf::SRTrack track;
      // fill shower variables
      track.Evis = part.depositions_sum;
      track.start = caf::SRVector3D(part.start_point[0], part.start_point[1], part.start_point[2]);
      track.end = caf::SRVector3D(part.end_point[0], part.end_point[1], part.end_point[2]); 
      track.dir = caf::SRVector3D(part.start_dir[0], part.start_dir[1], part.start_dir[2]);
      track.enddir = caf::SRVector3D(part.end_dir[0], part.end_dir[1], part.end_dir[2]);
      track.len_cm = sqrt(pow((part.start_point[0]-part.end_point[0]),2) + pow((part.start_point[1]-part.end_point[1]),2) + pow((part.start_point[2]-part.end_point[2]),2));
      track.truth.ixn = part.interaction_id;
      if(part.is_primary)track.truth.type = caf::TrueParticleID::kPrimary;
      track.truth.part = part.id;

      // note that interaction ID is not in general the same as the index within the sr.common.ixn.dlp vector
      // (some interaction IDs are filtered out as they're not beam triggers etc.)
      auto itIxn = std::find_if(sr.common.ixn.dlp.begin(), sr.common.ixn.dlp.end(),
                                [&part](const caf::SRInteraction & ixn){ return ixn.id == part.interaction_id; });
      if (itIxn == sr.common.ixn.dlp.end())
      {
        std::cerr << "ERROR: Particle's interaction ID (" << part.interaction_id << ") does not match any in the DLP set!\n";
        abort();
      }
      sr.nd.lar.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].tracks.push_back(std::move(track)); 
    }
  }

  // ------------------------------------------------------------------------------
  void MLNDLArRecoBranchFiller::FillShowers(const H5DataView<cafmaker::types::dlp::Particle> & particles,
                                            caf::StandardRecord &sr) const
  { 
    for (const auto & part : particles)
    {
      if (part.semantic_type != types::dlp::SemanticType::kShower)
        continue;

      caf::SRShower shower;
      // fill shower variables
      shower.Evis = part.depositions_sum;
      shower.start = caf::SRVector3D(part.start_point[0], part.start_point[1], part.start_point[2]);
      shower.direction = caf::SRVector3D(part.start_dir[0], part.start_dir[1], part.start_dir[2]);
      shower.truth.ixn = part.interaction_id;
      if(part.is_primary)shower.truth.type = caf::TrueParticleID::kPrimary;
      shower.truth.part = part.id;

      // note that interaction ID is not in general the same as the index within the sr.common.ixn.dlp vector
      // (some interaction IDs are filtered out as they're not beam triggers etc.)
      auto itIxn = std::find_if(sr.common.ixn.dlp.begin(), sr.common.ixn.dlp.end(),
                                [&part](const caf::SRInteraction & ixn){ return ixn.id == part.interaction_id; });
      if (itIxn == sr.common.ixn.dlp.end())
      {
        std::cerr << "ERROR: Particle's interaction ID (" << part.interaction_id << ") does not match any in the DLP set!\n";
        abort();
      }
      sr.nd.lar.dlp[std::distance(sr.common.ixn.dlp.begin(), itIxn)].showers.push_back(std::move(shower)); 
     
    }
  }

  // ------------------------------------------------------------------------------
  std::deque<Trigger> MLNDLArRecoBranchFiller::GetTriggers(int triggerType) const
  {
    auto runInfos = fDSReader.GetProducts<cafmaker::types::dlp::RunInfo>(-1); // get ALL the RunInfo products

    std::deque<Trigger> triggers;
    if (fTriggers.empty())
    {
      std::cout << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "' from " << runInfos.size() << " ND-LAr RunInfo products:\n";
      fTriggers.reserve(runInfos.size());
      for (const cafmaker::types::dlp::RunInfo &runInfo: runInfos)
      {
        const int placeholderTriggerType = 0;
        // todo: this check needs to be fixed in when we have trigger type info
        if (triggerType >= 0 && triggerType != placeholderTriggerType)
        {
          std::cout << "    skipping runinfo with event=" << runInfo.event << "\n";
          continue;
        }

        fTriggers.emplace_back();
        Trigger & trig = fTriggers.back();
        trig.evtID = runInfo.event;

        // todo: these are placeholder values until we can propagate enough info through the reco files
        trig.triggerType = 0;
        trig.triggerTime_s = runInfo.event;
        trig.triggerTime_ns = 0.;

        triggers.push_back(trig);

        std::cout << "  added trigger:  evtID=" << trig.evtID
                  << ", triggerType=" << trig.triggerType
                  << ", triggerTime_s=" << trig.triggerTime_s
                  << ", triggerTime_ns=" << trig.triggerTime_ns
                  << "\n";
      }
      fLastTriggerReqd = fTriggers.end();  // since we just modified the list, any iterators have been invalidated
    }
    else
    {
      for (const Trigger & trigger : fTriggers)
      {
        if (triggerType < 0 || triggerType == fTriggers.back().triggerType)
          triggers.push_back(trigger);
      }
    }

    return triggers;
  }

} // namespace cafmaker
