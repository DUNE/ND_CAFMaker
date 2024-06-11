/// Fill NDGAr reco branches using NDGAr reco data
///
/// \author  F. Martinez Lopez
/// \date    Oct. 2022
///

#include "NDGArRecoBranchFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include "Params.h"
#include "truth/FillTruth.h"

namespace cafmaker
{
  
  NDGArRecoBranchFiller::NDGArRecoBranchFiller(const std::string &NDGArRecoFilename)
    : IRecoBranchFiller("NDGAr"),
      fTriggers(),
      fLastTriggerReqd(fTriggers.end())
  {  
    fNDGArRecoFile = new TFile(NDGArRecoFilename.c_str(), "READ");
    if(!fNDGArRecoFile->IsZombie()){
      SetConfigured(true);

      // fTree = fNDGArRecoFile->Get<TTree>("anatree/GArAnaTree");
      NDGArRecoTree = dynamic_cast<TTree*>(fNDGArRecoFile->Get("anatree/GArAnaTree"));

      //NDGArRecoTree->SetBranchAddress("Run",    &fRun);
      //NDGArRecoTree->SetBranchAddress("SubRun", &fSubRun);
      NDGArRecoTree->SetBranchAddress("Event",  &fEvent);

      NDGArRecoTree->SetBranchAddress("GPartIdx", &fGPartIdx);

      NDGArRecoTree->SetBranchAddress("NType", &fMCNuPDG);
      NDGArRecoTree->SetBranchAddress("MCNuPx", &fMCNuPX);
      NDGArRecoTree->SetBranchAddress("MCNuPy", &fMCNuPY);
      NDGArRecoTree->SetBranchAddress("MCNuPz", &fMCNuPZ);

      NDGArRecoTree->SetBranchAddress("MCVertX", &fMCVertX);
      NDGArRecoTree->SetBranchAddress("MCVertY", &fMCVertY);
      NDGArRecoTree->SetBranchAddress("MCVertZ", &fMCVertZ);

      // MC Particles info
      NDGArRecoTree->SetBranchAddress("MCPTrkID", &fMCTrkID);
      NDGArRecoTree->SetBranchAddress("PDG", &fMCPDG);
      NDGArRecoTree->SetBranchAddress("MotherIndex", &fMCMotherIndex);
      NDGArRecoTree->SetBranchAddress("PDGMother", &fMCPDGMother);

      NDGArRecoTree->SetBranchAddress("MCPStartX", &fMCStartX);
      NDGArRecoTree->SetBranchAddress("MCPStartY", &fMCStartY);
      NDGArRecoTree->SetBranchAddress("MCPStartZ", &fMCStartZ);
      NDGArRecoTree->SetBranchAddress("MCPStartPX", &fMCStartPX);
      NDGArRecoTree->SetBranchAddress("MCPStartPY", &fMCStartPY);
      NDGArRecoTree->SetBranchAddress("MCPStartPZ", &fMCStartPZ);

      //NDGArRecoTree->SetBranchAddress("MCPEndX", &fMCEndX);
      //NDGArRecoTree->SetBranchAddress("MCPEndY", &fMCEndY);
      //NDGArRecoTree->SetBranchAddress("MCPEndZ", &fMCEndZ);
      //NDGArRecoTree->SetBranchAddress("MCPEndPX", &fMCEndPX);
      //NDGArRecoTree->SetBranchAddress("MCPEndPY", &fMCEndPY);
      //NDGArRecoTree->SetBranchAddress("MCPEndPZ", &fMCEndPZ);

      // Track-related info
      NDGArRecoTree->SetBranchAddress("TrackStartX", &fTrackStartX);
      NDGArRecoTree->SetBranchAddress("TrackStartY", &fTrackStartY);
      NDGArRecoTree->SetBranchAddress("TrackStartZ", &fTrackStartZ);
      NDGArRecoTree->SetBranchAddress("TrackEndX", &fTrackEndX);
      NDGArRecoTree->SetBranchAddress("TrackEndY", &fTrackEndY);
      NDGArRecoTree->SetBranchAddress("TrackEndZ", &fTrackEndZ);

      NDGArRecoTree->SetBranchAddress("TrackStartPX", &fTrackStartPX);
      NDGArRecoTree->SetBranchAddress("TrackStartPY", &fTrackStartPY);
      NDGArRecoTree->SetBranchAddress("TrackStartPZ", &fTrackStartPZ);
      NDGArRecoTree->SetBranchAddress("TrackEndPX", &fTrackEndPX);
      NDGArRecoTree->SetBranchAddress("TrackEndPY", &fTrackEndPY);
      NDGArRecoTree->SetBranchAddress("TrackEndPZ", &fTrackEndPZ);

      NDGArRecoTree->SetBranchAddress("TrackLenF", &fTrackLenF);
      NDGArRecoTree->SetBranchAddress("TrackLenB", &fTrackLenB);
      //NDGArRecoTree->SetBranchAddress("TrackPF", &fTrackPF);
      //NDGArRecoTree->SetBranchAddress("TrackPB", &fTrackPB);
      NDGArRecoTree->SetBranchAddress("TrackAvgIonF", &fTrackAvgIonF);
      NDGArRecoTree->SetBranchAddress("TrackAvgIonB", &fTrackAvgIonB);

      NDGArRecoTree->SetBranchAddress("TrackIDNumber", &fTrackIDNumber);
      NDGArRecoTree->SetBranchAddress("NTPCClustersOnTrack", &fTrackNClusters);

      NDGArRecoTree->SetBranchAddress("TrackPIDF", &fTrackPIDF);
      NDGArRecoTree->SetBranchAddress("TrackPIDProbF", &fTrackPIDProbF);
      NDGArRecoTree->SetBranchAddress("TrackPIDB", &fTrackPIDB);
      NDGArRecoTree->SetBranchAddress("TrackPIDProbB", &fTrackPIDProbB);

      NDGArRecoTree->SetBranchAddress("TrackMCindex", &fTrackMCindex);
      NDGArRecoTree->SetBranchAddress("TrackMCfrac", &fTrackMCfrac);

      // ECAL-related info
      NDGArRecoTree->SetBranchAddress("ClusterX", &fECALClusterX);
      NDGArRecoTree->SetBranchAddress("ClusterY", &fECALClusterY);
      NDGArRecoTree->SetBranchAddress("ClusterZ", &fECALClusterZ);

      NDGArRecoTree->SetBranchAddress("ClusterIDNumber", &fECALClusterIDNumber);

      NDGArRecoTree->SetBranchAddress("ClusterEnergy", &fECALClusterEnergy);
      NDGArRecoTree->SetBranchAddress("ClusterNhits",  &fECALClusterNhits);

      NDGArRecoTree->SetBranchAddress("ClusterMCindex", &fECALClusterMCindex);
      NDGArRecoTree->SetBranchAddress("ClusterMCfrac",  &fECALClusterMCfrac);

      // ECAL-track associations
      NDGArRecoTree->SetBranchAddress("ECALAssn_ClusIDNumber",  &fECALAssn_ClusterID);
      NDGArRecoTree->SetBranchAddress("ECALAssn_TrackIDNumber", &fECALAssn_TrackID);

      // MuID-related info
      NDGArRecoTree->SetBranchAddress("ClusterX_MuID", &fMuIDClusterX);
      NDGArRecoTree->SetBranchAddress("ClusterY_MuID", &fMuIDClusterY);
      NDGArRecoTree->SetBranchAddress("ClusterZ_MuID", &fMuIDClusterZ);

      NDGArRecoTree->SetBranchAddress("ClusterIDNumber_MuID", &fMuIDClusterIDNumber);

      NDGArRecoTree->SetBranchAddress("ClusterEnergy_MuID", &fMuIDClusterEnergy);
      NDGArRecoTree->SetBranchAddress("ClusterNhits_MuID",  &fMuIDClusterNhits);

      NDGArRecoTree->SetBranchAddress("ClusterMCindex_MuID", &fMuIDClusterMCindex);
      NDGArRecoTree->SetBranchAddress("ClusterMCfrac_MuID",  &fMuIDClusterMCfrac);

      // MuID-track associations
      NDGArRecoTree->SetBranchAddress("MuIDAssn_ClusIDNumber",  &fMuIDAssn_ClusterID);
      NDGArRecoTree->SetBranchAddress("MuIDAssn_TrackIDNumber", &fMuIDAssn_TrackID);

      // Reco particle info
      NDGArRecoTree->SetBranchAddress("RecoParticleIDNumber",    &fRecoParticleID);

      NDGArRecoTree->SetBranchAddress("RecoParticleMomentum",    &fRecoParticleMomentum);

      NDGArRecoTree->SetBranchAddress("RecoParticleTotalCaloEnergy",    &fRecoParticleTotalCaloEnergy);
      NDGArRecoTree->SetBranchAddress("RecoParticleMeanCaloEnergy",     &fRecoParticleMeanCaloEnergy);

      NDGArRecoTree->SetBranchAddress("RecoParticleTotalECALEnergy",  &fRecoParticleTotalECALEnergy);
      NDGArRecoTree->SetBranchAddress("RecoParticleNHitsECAL",        &fRecoParticleNHitsECAL);
      NDGArRecoTree->SetBranchAddress("RecoParticleTotalMuIDEnergy",  &fRecoParticleTotalMuIDEnergy);
      NDGArRecoTree->SetBranchAddress("RecoParticleNHitsMuID",        &fRecoParticleNHitsMuID);

      NDGArRecoTree->SetBranchAddress("RecoParticleECALToFTime",    &fRecoParticleToFTime);
      NDGArRecoTree->SetBranchAddress("RecoParticleECALToFBeta",    &fRecoParticleToFBeta);

      NDGArRecoTree->SetBranchAddress("RecoParticleCharge",      &fRecoParticleCharge);

      NDGArRecoTree->SetBranchAddress("RecoParticleStartX",      &fRecoParticleStartX);
      NDGArRecoTree->SetBranchAddress("RecoParticleStartY",      &fRecoParticleStartY);
      NDGArRecoTree->SetBranchAddress("RecoParticleStartZ",      &fRecoParticleStartZ);
      NDGArRecoTree->SetBranchAddress("RecoParticleEndX",        &fRecoParticleEndX);
      NDGArRecoTree->SetBranchAddress("RecoParticleEndY",        &fRecoParticleEndY);
      NDGArRecoTree->SetBranchAddress("RecoParticleEndZ",        &fRecoParticleEndZ);
      NDGArRecoTree->SetBranchAddress("RecoParticleDirectionX",  &fRecoParticleDirectionX);
      NDGArRecoTree->SetBranchAddress("RecoParticleDirectionY",  &fRecoParticleDirectionY);
      NDGArRecoTree->SetBranchAddress("RecoParticleDirectionZ",  &fRecoParticleDirectionZ);

      NDGArRecoTree->SetBranchAddress("RecoParticleMuonScore",        &fRecoParticleMuonScore);
      NDGArRecoTree->SetBranchAddress("RecoParticleProtonCaloScore",  &fRecoParticleProtonCaloScore);
      NDGArRecoTree->SetBranchAddress("RecoParticleProtonToFScore",   &fRecoParticleProtonToFScore);

      NDGArRecoTree->SetBranchAddress("RecoParticlePID",     &fRecoParticlePID);
      NDGArRecoTree->SetBranchAddress("RecoParticleEnergy",  &fRecoParticleEnergy);

      NDGArRecoTree->SetBranchAddress("RecoParticleMCindex", &fRecoParticleMCindex);
      NDGArRecoTree->SetBranchAddress("RecoParticleMCfrac",  &fRecoParticleMCfrac);

      // Reco particle associations
      NDGArRecoTree->SetBranchAddress("RecoParticleTrackAssn_RecoPIDNumber",     &fParticleTrackAssn_ParticleID);
      NDGArRecoTree->SetBranchAddress("RecoParticleTrackAssn_TrackIDNumber",     &fParticleTrackAssn_TrackID);

      NDGArRecoTree->SetBranchAddress("RecoParticleECalAssn_RecoPIDNumber",      &fParticleECALAssn_ParticleID);
      NDGArRecoTree->SetBranchAddress("RecoParticleECalAssn_ClusterIDNumber",    &fParticleECALAssn_ClusterID);

      NDGArRecoTree->SetBranchAddress("RecoParticleMuIDAssn_RecoPIDNumber",      &fParticleMuIDAssn_ParticleID);
      NDGArRecoTree->SetBranchAddress("RecoParticleMuIDAssn_ClusterIDNumber",    &fParticleMuIDAssn_ClusterID);

      // Reco interaction info
      NDGArRecoTree->SetBranchAddress("RecoNuEnergy",  &fRecoNuEnergy);

    } else {
      NDGArRecoTree = NULL;
      std::cerr << "Did not find input ND-GAr reco file you provided: " << NDGArRecoFilename << std::endl;
      std::cerr << "Are you sure it exists?" << std::endl;
      throw;
    }
  }

  namespace
  {
    // F says: slightly modified version of the SRPartCmp from
    //         MLNDLArRecoBranchFiller, to account for the fact
    //         that we didn't store the particle energy but the
    //         3-momentum
    struct SRPartCmp
    {
      float p;
      int trkid;
      bool operator()(const caf::SRTrueParticle & part) const
      {
        LOG_S("SRPartCmp").VERBOSE() << "       SRPartCmp::operator()():  looking for p = " << p << "; this particle p = " << part.p.Mag() << ","
                                     << "trk ID = " << trkid << ", this particle trkID = " << part.G4ID << "\n";
        return part.p.Mag() == p && (trkid < 0 || part.G4ID < 0 || trkid == part.G4ID);
      }
    };
  }

  // ------------------------------------------------------------------------------
  void
  NDGArRecoBranchFiller::_FillRecoBranches(const Trigger &trigger,
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

    NDGArRecoTree->GetEntry(idx);

    //Fill ND-GAr specific info in the meta branch
    sr.meta.nd_gar.enabled = true;
    //sr.meta.nd_gar.run = fRun;
    //sr.meta.nd_gar.subrun = fSubRun;
    LOG.VERBOSE() << "    GArSoft event number: " << fEvent << ".\n";
    sr.meta.nd_gar.event = fEvent;

    FillInteractions(truthMatcher, sr);

    FillTracks(truthMatcher, sr);
    FillECalClusters(truthMatcher, sr);
    FillMuIDClusters(truthMatcher, sr);
    FillAssociations(sr);

    FillParticles(truthMatcher, sr);
    LOG.VERBOSE() << "    Done with event: " << fEvent << ".\n";

    //legacy variables
    //using forward variables only!

    // int pi_pl_mult = 0;
    // int pi_min_mult = 0;
    //
    // for (size_t i=0; i< fTrackStartX->size(); ++i){
    //    sr.nd.gar.pdg.push_back(fTrackPIDCheatedF->at(i));
    //    sr.nd.gar.trkLen.push_back(fTrackLenF->at(i));
    //    sr.nd.gar.ptrue.push_back(fTrackPF->at(i));       //filling with "reco" momentum for the moment
    //    if (fTrackPIDCheatedF->at(i) == 211){
    //      ++pi_pl_mult;
    //    } else if (fTrackPIDCheatedF->at(i) == -211){
    //      ++pi_min_mult;
    //    }
    // }
    // sr.nd.gar.gastpc_pi_pl_mult = pi_pl_mult;
    // sr.nd.gar.gastpc_pi_min_mult = pi_min_mult;

  }

  // ------------------------------------------------------------------------------
  void NDGArRecoBranchFiller::FillTrueInteraction(caf::SRTrueInteraction & srTrueInt) const
  {

    const auto NaN = std::numeric_limits<float>::signaling_NaN();

    srTrueInt.vtx.x = fMCVertX->at(0);
    srTrueInt.vtx.y = fMCVertY->at(0);
    srTrueInt.vtx.z = fMCVertZ->at(0);

    LOG.DEBUG() << "    true nu vertex: (" << srTrueInt.vtx.x << ", " << srTrueInt.vtx.y << ", " << srTrueInt.vtx.z <<")\n.";

  }

  // ------------------------------------------------------------------------------
  void NDGArRecoBranchFiller::FillInteractions(const TruthMatcher * truthMatch,
                                               caf::StandardRecord &sr) const
  {
    // F says: our GArSoft samples contain one neutrino interaction
    //         per event/trigger (I think we can easily change that?)
    //         so we only need one interaction object in the common
    //         reco branch so far
    sr.common.ixn.gsft.reserve(1);
    sr.common.ixn.ngsft = 1;

    // F says: again, right now we have one interaction per event/trigger
    //         so we simply need to create an SRInteraction and assign a
    //         dummy id of 1 to it

    caf::SRInteraction interaction;
    interaction.id  = 1;

    // F says: no reco interaction vertex for now, don't fill
    /* interaction.vtx  = caf::SRVector3D(ixn.vertex[0], ixn.vertex[1], ixn.vertex[2]); */

    interaction.Enu.lep_calo = fRecoNuEnergy->at(0); // only one true interaction per trigger, only one reco energy in vector...

    // F says: we store all the GENIE PartIdx in our GArSoft tree,
    //         so we can simply take the first one 
    caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, 10000000000+fEvent-1, true);

    LOG.VERBOSE() << "    --> resulting SRTrueInteraction has the following particles in it:\n";
    for (const caf::SRTrueParticle & part : srTrueInt.prim)
      LOG.VERBOSE() << "    (prim) pdg = " << part.pdg << ", energy = " << part.p.E << ", G4ID = " << part.G4ID << "\n";
    for (const caf::SRTrueParticle & part : srTrueInt.prefsi)
      LOG.VERBOSE() << "    (prefsi) pdg = " << part.pdg << ", energy = " << part.p.E << ", G4ID = " << part.G4ID << "\n";
    for (const caf::SRTrueParticle & part : srTrueInt.sec)
      LOG.VERBOSE() << "    (sec) pdg = " << part.pdg << ", energy = " << part.p.E << ", G4ID = " << part.G4ID << "\n";

    LOG.VERBOSE() << "    is CC = " << srTrueInt.iscc << "\n";

    // here we need to fill in any additional info
    // that GENIE didn't know about: e.g., secondary particles made by GEANT4

    // F says: for some reason this also fills the vertex position (?)
    FillTrueInteraction(srTrueInt);

    // note that the interaction ID is GENIE's label for it, which may not be the same as the index in the vector
    std::size_t truthVecIdx = std::distance(sr.mc.nu.begin(),
                                            std::find_if(sr.mc.nu.begin(),
                                                          sr.mc.nu.end(),
                                                          [&srTrueInt](const caf::SRTrueInteraction& nu)
                                                          {
                                                            return nu.id == srTrueInt.id;
                                                          }));

    interaction.truth.push_back(truthVecIdx);

    LOG.VERBOSE() << "  ** end matched true interaction search" << ".\n";

    sr.common.ixn.gsft.push_back(std::move(interaction));

  }

  // ------------------------------------------------------------------------------
  void NDGArRecoBranchFiller::FillParticles(const TruthMatcher * truthMatch,
                                            caf::StandardRecord &sr) const
  {
    LOG.DEBUG() << "Filling reco particles...\n";

    // note: used in the hack further below
    //static SRPartCmp srPartCmp;

    size_t n_particles = fTrackStartX->size();
    LOG.VERBOSE() << "    GArSoft number of reco particles: " << n_particles << ".\n";
    //sr.nd.gar.ixn[0].nparticles = n_particles;

    for (size_t iParticle=0; iParticle<n_particles; iParticle++){
        LOG.VERBOSE() << "        Filling particle " << iParticle << ".\n";

      // F says: we create both a SRRecoParticle (which will go in sr.common.ixn.gsft.part)
      //         and a SRGArParticle (which will go in sr.nd.gar.ixn.particles) at the same
      //         time. They both contain almost the same information, but SRGArParticle has
      //         all the variables required to recompute the PID.
      caf::SRRecoParticle reco_particle;

      // if(part.is_primary) reco_particle.primary = true;
      //reco_particle.start = caf::SRVector3D(part.start_point[0], part.start_point[1], part.start_point[2]);
      //reco_particle.end = caf::SRVector3D(part.end_point[0], part.end_point[1], part.end_point[2]);
      //reco_particle.contained = part.is_contained; // this is not just the vertex, but all energies are contained

      reco_particle.tgtA = 40;

      reco_particle.pdg = fRecoParticlePID->at(iParticle);

      float momentum = fRecoParticleMomentum->at(iParticle);
      reco_particle.p.x = momentum*fRecoParticleDirectionX->at(iParticle);
      reco_particle.p.y = momentum*fRecoParticleDirectionY->at(iParticle);
      reco_particle.p.z = momentum*fRecoParticleDirectionZ->at(iParticle);

      reco_particle.start = caf::SRVector3D(fRecoParticleStartX->at(iParticle),
                                            fRecoParticleStartY->at(iParticle),
                                            fRecoParticleStartZ->at(iParticle));

      reco_particle.end   = caf::SRVector3D(fRecoParticleEndX->at(iParticle),
                                            fRecoParticleEndY->at(iParticle),
                                            fRecoParticleEndZ->at(iParticle));
      
      reco_particle.E = fRecoParticleEnergy->at(iParticle);
      reco_particle.E_method = caf::PartEMethod::kCurvature; // only method available for GAr atm

      reco_particle.score.gsft_pid.muon_score        = fRecoParticleMuonScore->at(iParticle);
      reco_particle.score.gsft_pid.proton_dEdx_score = fRecoParticleProtonCaloScore->at(iParticle);
      reco_particle.score.gsft_pid.proton_tof_score  = fRecoParticleProtonToFScore->at(iParticle);

      caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, 10000000000+fEvent-1, false);
      // we need this below because caf::TrueParticleID wants the *index* of the SRTrueInteraction
      int srTrueIntIdx = std::distance(sr.mc.nu.begin(),
                                       std::find_if(sr.mc.nu.begin(),
                                                    sr.mc.nu.end(),
                                                    [&srTrueInt](const caf::SRTrueInteraction& ixn) {return ixn.id == srTrueInt.id;}));

      int iMCParticleRecoP = fRecoParticleMCindex->at(iParticle);
      if(iMCParticleRecoP >= 0){

        caf::SRVector3D mc_p(fMCStartPX->at(iMCParticleRecoP), fMCStartPY->at(iMCParticleRecoP), fMCStartPZ->at(iMCParticleRecoP));


        // find the true particle this reco particle goes with.
        // if we had GENIE info and it was a primary, it should already be filled in.
        // we use the comparison version because the G4ID from the pass-through
        // counts up monotonically from 0 across the whole FILE,
        // whereas the GENIE events start over at every interaction.
        // moreover, the cafmaker::types::gsft::TrueParticle::is_primary flag
        // is currently broken (upstream info from Supera is screwed up)
        // so we need to try both collections :(
        static SRPartCmp srPartCmp;
        srPartCmp.p = mc_p.Mag();
        srPartCmp.trkid = fMCTrkID->at(iMCParticleRecoP);

        bool isPrim = false;
        caf::SRTrueParticle * srTruePart = nullptr;
        try
        {
          // F says: our GENIE particles won't have a valid G4ID
          //         therefore the comparsion will only look at 
          //         the module of the 3-momentum
          srTruePart = &truthMatch->GetTrueParticle(sr, srTrueInt, srPartCmp, true, false);
          isPrim = true;
        }
        catch ( std::runtime_error& err )
        {
          // guess if it wasn't a primary, it must be a secondary :(
          srTruePart = &truthMatch->GetTrueParticle(sr, srTrueInt, srPartCmp, false, true);
        }

        // however this will fill in any other fields that weren't copied from a GENIE record
        // (which also handles the case where this particle is a secondary)
        FillTruth(*srTruePart, iMCParticleRecoP);

        // the particle idx is within the GENIE vector, which may not be the same as the index in the vector here
        // first find the interaction that it goes with
        LOG.VERBOSE() << "      this particle is " << (isPrim ? "PRIMARY" : "SECONDARY") << "\n";
        std::vector<caf::SRTrueParticle> & collection = (isPrim)
                                                        ? srTrueInt.prim
                                                        : srTrueInt.sec;
        std::size_t truthVecIdx = std::distance(collection.begin(),
                                                std::find_if(collection.begin(),
                                                              collection.end(),
                                                              srPartCmp));

        LOG.VERBOSE() << "      index of SRParticle in the SRInteraction " << truthVecIdx << "\n";

        // F says: currently our reco only gives truth matching
        //         to one MCParticle, but we can change that if
        //         we need to...
        reco_particle.truth.push_back(caf::TrueParticleID{srTrueIntIdx,
                                                          (isPrim) ? caf::TrueParticleID::PartType::kPrimary :  caf::TrueParticleID::PartType::kSecondary,
                                                          static_cast<int>(truthVecIdx)});

        reco_particle.truthOverlap.push_back(fRecoParticleMCfrac->at(iParticle));

      }

      sr.common.ixn.gsft[0].part.gsft.push_back(std::move(reco_particle));
      sr.common.ixn.gsft[0].part.ngsft++;

    }
  }

  // ------------------------------------------------------------------------------
  void NDGArRecoBranchFiller::FillTracks(const TruthMatcher * truthMatch,
                                         caf::StandardRecord &sr) const
  {
    // F says: sr.nd.gar.ixn should be resize using the size of
    //         sr.common.ixn.gsft (why? isn't sr.common.ixn.ngsft
    //         the same?), but because it's empty for now let's use 
    //         sr.common.ixn.ngsft
    sr.nd.gar.ixn.resize(sr.common.ixn.ngsft);
    sr.nd.gar.nixn = sr.common.ixn.ngsft;

    size_t n_tracks = fTrackStartX->size();
    LOG.VERBOSE() << "    GArSoft number of reco tracks: " << n_tracks << ".\n";
    sr.nd.gar.ixn[0].ntracks = n_tracks;

    //size_t pid_counter = 0;
    caf::SRGArTrack track;
    for (size_t iTrack=0; iTrack<n_tracks; iTrack++){
      LOG.VERBOSE() << "        Filling track " << iTrack << ".\n";
      
      caf::SRVector3D start(fTrackStartX->at(iTrack), fTrackStartY->at(iTrack), fTrackStartZ->at(iTrack));
      track.start = start;
      
      caf::SRVector3D end(fTrackEndX->at(iTrack), fTrackEndY->at(iTrack), fTrackEndZ->at(iTrack));
      track.end = end;
      
      caf::SRVector3D dir(fTrackStartPX->at(iTrack), fTrackStartPY->at(iTrack), fTrackStartPZ->at(iTrack));
      track.dir = dir.Unit();

      caf::SRVector3D enddir(fTrackEndPX->at(iTrack), fTrackEndPY->at(iTrack), fTrackEndPZ->at(iTrack));
      track.enddir = enddir.Unit();

      track.len_cm_fwd = fTrackLenF->at(iTrack);
      track.len_cm_bak = fTrackLenB->at(iTrack);
      //track.p_fwd = fTrackPF->at(iTrack);
      //track.p_bak = fTrackPB->at(iTrack);
      track.dQdx_fwd = fTrackAvgIonF->at(iTrack);
      track.dQdx_bak = fTrackAvgIonB->at(iTrack);

      // F says: Names are a mess, soz...
      track.dEdx_fwd  = fRecoParticleMeanCaloEnergy->at(iTrack);
      //track.dEdx_bak = 
      track.dEdx_total = fRecoParticleTotalCaloEnergy->at(iTrack);

      track.garsoft_trk_id = fTrackIDNumber->at(iTrack);
      track.clusters_in_track = fTrackNClusters->at(iTrack);

      // Currently not working in reco, will uncomment when solved
      /* LOG.VERBOSE() << "        PID counter: " << pid_counter << ".\n";
      for (size_t iPID=6*pid_counter; iPID<6*(pid_counter+1); ++iPID){
        LOG.VERBOSE() << "            iPID in loop: " << iPID << ".\n";
        track.pid_fwd.push_back(fTrackPIDF->at(iPID));
        track.pid_prob_fwd.push_back(fTrackPIDProbF->at(iPID));
        track.pid_bak.push_back(fTrackPIDB->at(iPID));
        track.pid_prob_bak.push_back(fTrackPIDProbB->at(iPID));
      }
      ++pid_counter; */

      caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, 10000000000+fEvent-1, false);
      // we need this below because caf::TrueParticleID wants the *index* of the SRTrueInteraction
      int srTrueIntIdx = std::distance(sr.mc.nu.begin(),
                                        std::find_if(sr.mc.nu.begin(),
                                                    sr.mc.nu.end(),
                                                    [&srTrueInt](const caf::SRTrueInteraction& ixn) {return ixn.id == srTrueInt.id;}));

      int iMCParticleTrack = fTrackMCindex->at(iTrack);
      if(iMCParticleTrack >= 0){

        caf::SRVector3D mc_p(fMCStartPX->at(iMCParticleTrack), fMCStartPY->at(iMCParticleTrack), fMCStartPZ->at(iMCParticleTrack));


        // find the true particle this reco particle goes with.
        // if we had GENIE info and it was a primary, it should already be filled in.
        // we use the comparison version because the G4ID from the pass-through
        // counts up monotonically from 0 across the whole FILE,
        // whereas the GENIE events start over at every interaction.
        // moreover, the cafmaker::types::gsft::TrueParticle::is_primary flag
        // is currently broken (upstream info from Supera is screwed up)
        // so we need to try both collections :(
        static SRPartCmp srPartCmp;
        srPartCmp.p = mc_p.Mag();
        srPartCmp.trkid = fMCTrkID->at(iMCParticleTrack);

        bool isPrim = false;
        caf::SRTrueParticle * srTruePart = nullptr;
        try
        {
          // F says: our GENIE particles won't have a valid G4ID
          //         therefore the comparsion will only look at 
          //         the module of the 3-momentum
          srTruePart = &truthMatch->GetTrueParticle(sr, srTrueInt, srPartCmp, true, false);
          isPrim = true;
        }
        catch ( std::runtime_error& err )
        {
          // guess if it wasn't a primary, it must be a secondary :(
          srTruePart = &truthMatch->GetTrueParticle(sr, srTrueInt, srPartCmp, false, true);
        }

        // however this will fill in any other fields that weren't copied from a GENIE record
        // (which also handles the case where this particle is a secondary)
        FillTruth(*srTruePart, iMCParticleTrack);

        // the particle idx is within the GENIE vector, which may not be the same as the index in the vector here
        // first find the interaction that it goes with
        LOG.VERBOSE() << "      this particle is " << (isPrim ? "PRIMARY" : "SECONDARY") << "\n";
        std::vector<caf::SRTrueParticle> & collection = (isPrim)
                                                        ? srTrueInt.prim
                                                        : srTrueInt.sec;
        std::size_t truthVecIdx = std::distance(collection.begin(),
                                                std::find_if(collection.begin(),
                                                              collection.end(),
                                                              srPartCmp));

        LOG.VERBOSE() << "      index of SRParticle in the SRInteraction " << truthVecIdx << "\n";

        // F says: currently our reco only gives truth matching
        //         to one MCParticle, but we can change that if
        //         we need to...
        track.truth.push_back(caf::TrueParticleID{srTrueIntIdx,
                                          (isPrim) ? caf::TrueParticleID::PartType::kPrimary :  caf::TrueParticleID::PartType::kSecondary,
                                          static_cast<int>(truthVecIdx)});
        
        track.truthOverlap.push_back(fTrackMCfrac->at(iTrack));

      }

      sr.nd.gar.ixn[0].tracks.push_back(track);

    }
  }

  // ------------------------------------------------------------------------------
  void NDGArRecoBranchFiller::FillTruth(caf::SRTrueParticle & srTruePart, size_t iMCParticleTrack) const
  {

    srTruePart.G4ID = fMCTrkID->at(iMCParticleTrack);

    srTruePart.start_pos.x = fMCStartX->at(iMCParticleTrack);
    srTruePart.start_pos.y = fMCStartY->at(iMCParticleTrack);
    srTruePart.start_pos.z = fMCStartZ->at(iMCParticleTrack);

    srTruePart.p.px = fMCStartPX->at(iMCParticleTrack);
    srTruePart.p.py = fMCStartPY->at(iMCParticleTrack);
    srTruePart.p.pz = fMCStartPZ->at(iMCParticleTrack);

    srTruePart.parent = fMCMotherIndex->at(iMCParticleTrack);

  }

  // ------------------------------------------------------------------------------
  void NDGArRecoBranchFiller::FillECalClusters(const TruthMatcher * truthMatch,
                                               caf::StandardRecord &sr) const
  {
    // F says: sr.nd.gar.ixn already has the right size, so
    //         no need to worry about that again here

    size_t n_clusters = fECALClusterX->size();
    LOG.VERBOSE() << "    GArSoft number of ECal clusters: " << n_clusters << ".\n";
    sr.nd.gar.ixn[0].necalclusters = n_clusters;

    size_t n_assns = fECALAssn_ClusterID->size();
    caf::SRGArCalo cluster;
    cluster.cluster_type = caf::GArCaloType::kECalCluster;
    for (size_t iECAL=0; iECAL<n_clusters; iECAL++){
        LOG.VERBOSE() << "        Filling cluster " << iECAL << ".\n";

        caf::SRVector3D position(fECALClusterX->at(iECAL), fECALClusterY->at(iECAL), fECALClusterZ->at(iECAL));
        cluster.position = position;

        cluster.E = fECALClusterEnergy->at(iECAL);
        cluster.hits_in_cluster = fECALClusterNhits->at(iECAL);

        cluster.garsoft_calo_id = fECALClusterIDNumber->at(iECAL);

        for (size_t iAssn=0; iAssn<n_assns; ++iAssn){
          LOG.VERBOSE() << "            iAssn in loop: " << iAssn << ".\n";
          if (cluster.garsoft_calo_id == fECALAssn_ClusterID->at(iAssn)){
            cluster.garsoft_trk_assn = fECALAssn_TrackID->at(iAssn);
          }
        }

      caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, 10000000000+fEvent-1, false);
      // we need this below because caf::TrueParticleID wants the *index* of the SRTrueInteraction
      int srTrueIntIdx = std::distance(sr.mc.nu.begin(),
                                        std::find_if(sr.mc.nu.begin(),
                                                    sr.mc.nu.end(),
                                                    [&srTrueInt](const caf::SRTrueInteraction& ixn) {return ixn.id == srTrueInt.id;}));

      int iMCParticleECAL = fECALClusterMCindex->at(iECAL);
      if(iMCParticleECAL >= 0){

        caf::SRVector3D mc_p(fMCStartPX->at(iMCParticleECAL), fMCStartPY->at(iMCParticleECAL), fMCStartPZ->at(iMCParticleECAL));


        // find the true particle this reco particle goes with.
        // if we had GENIE info and it was a primary, it should already be filled in.
        // we use the comparison version because the G4ID from the pass-through
        // counts up monotonically from 0 across the whole FILE,
        // whereas the GENIE events start over at every interaction.
        // moreover, the cafmaker::types::gsft::TrueParticle::is_primary flag
        // is currently broken (upstream info from Supera is screwed up)
        // so we need to try both collections :(
        static SRPartCmp srPartCmp;
        srPartCmp.p = mc_p.Mag();
        srPartCmp.trkid = fMCTrkID->at(iMCParticleECAL);

        bool isPrim = false;
        caf::SRTrueParticle * srTruePart = nullptr;
        try
        {
          // F says: our GENIE particles won't have a valid G4ID
          //         therefore the comparsion will only look at 
          //         the module of the 3-momentum
          srTruePart = &truthMatch->GetTrueParticle(sr, srTrueInt, srPartCmp, true, false);
          isPrim = true;
        }
        catch ( std::runtime_error& err )
        {
          // guess if it wasn't a primary, it must be a secondary :(
          srTruePart = &truthMatch->GetTrueParticle(sr, srTrueInt, srPartCmp, false, true);
        }

        // however this will fill in any other fields that weren't copied from a GENIE record
        // (which also handles the case where this particle is a secondary)
        FillTruth(*srTruePart, iMCParticleECAL);

        // the particle idx is within the GENIE vector, which may not be the same as the index in the vector here
        // first find the interaction that it goes with
        LOG.VERBOSE() << "      this particle is " << (isPrim ? "PRIMARY" : "SECONDARY") << "\n";
        std::vector<caf::SRTrueParticle> & collection = (isPrim)
                                                        ? srTrueInt.prim
                                                        : srTrueInt.sec;
        std::size_t truthVecIdx = std::distance(collection.begin(),
                                                std::find_if(collection.begin(),
                                                              collection.end(),
                                                              srPartCmp));

        LOG.VERBOSE() << "      index of SRParticle in the SRInteraction " << truthVecIdx << "\n";

        // F says: currently our reco only gives truth matching
        //         to one MCParticle, but we can change that if
        //         we need to...
        cluster.truth.push_back(caf::TrueParticleID{srTrueIntIdx,
                                          (isPrim) ? caf::TrueParticleID::PartType::kPrimary :  caf::TrueParticleID::PartType::kSecondary,
                                          static_cast<int>(truthVecIdx)});
        
        cluster.truthOverlap.push_back(fECALClusterMCfrac->at(iECAL));

      }

        sr.nd.gar.ixn[0].ecalclusters.push_back(cluster);

    }

  }

  // ------------------------------------------------------------------------------
  void NDGArRecoBranchFiller::FillMuIDClusters(const TruthMatcher * truthMatch,
                                               caf::StandardRecord &sr) const
  {
    // F says: sr.nd.gar.ixn already has the right size, so
    //         no need to worry about that again here

    size_t n_clusters = fMuIDClusterX->size();
    LOG.VERBOSE() << "    GArSoft number of MuID clusters: " << n_clusters << ".\n";
    sr.nd.gar.ixn[0].nmuidclusters = n_clusters;

    size_t n_assns = fMuIDAssn_ClusterID->size();
    caf::SRGArCalo cluster;
    cluster.cluster_type = caf::GArCaloType::kMuIDCluster;
    for (size_t iMuID=0; iMuID<n_clusters; iMuID++){
        LOG.VERBOSE() << "        Filling cluster " << iMuID << ".\n";

        caf::SRVector3D position(fMuIDClusterX->at(iMuID), fMuIDClusterY->at(iMuID), fMuIDClusterZ->at(iMuID));
        cluster.position = position;

        cluster.E = fMuIDClusterEnergy->at(iMuID);
        cluster.hits_in_cluster = fMuIDClusterNhits->at(iMuID);

        cluster.garsoft_calo_id = fMuIDClusterIDNumber->at(iMuID);

        for (size_t iAssn=0; iAssn<n_assns; ++iAssn){
          LOG.VERBOSE() << "            iAssn in loop: " << iAssn << ".\n";
          if (cluster.garsoft_calo_id == fMuIDAssn_ClusterID->at(iAssn)){
            cluster.garsoft_trk_assn = fMuIDAssn_TrackID->at(iAssn);
          }
        }

      caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, 10000000000+fEvent-1, false);
      // we need this below because caf::TrueParticleID wants the *index* of the SRTrueInteraction
      int srTrueIntIdx = std::distance(sr.mc.nu.begin(),
                                        std::find_if(sr.mc.nu.begin(),
                                                    sr.mc.nu.end(),
                                                    [&srTrueInt](const caf::SRTrueInteraction& ixn) {return ixn.id == srTrueInt.id;}));

      int iMCParticleMuID = fMuIDClusterMCindex->at(iMuID);
      if(iMCParticleMuID >= 0){

        caf::SRVector3D mc_p(fMCStartPX->at(iMCParticleMuID), fMCStartPY->at(iMCParticleMuID), fMCStartPZ->at(iMCParticleMuID));


        // find the true particle this reco particle goes with.
        // if we had GENIE info and it was a primary, it should already be filled in.
        // we use the comparison version because the G4ID from the pass-through
        // counts up monotonically from 0 across the whole FILE,
        // whereas the GENIE events start over at every interaction.
        // moreover, the cafmaker::types::gsft::TrueParticle::is_primary flag
        // is currently broken (upstream info from Supera is screwed up)
        // so we need to try both collections :(
        static SRPartCmp srPartCmp;
        srPartCmp.p = mc_p.Mag();
        srPartCmp.trkid = fMCTrkID->at(iMCParticleMuID);

        bool isPrim = false;
        caf::SRTrueParticle * srTruePart = nullptr;
        try
        {
          // F says: our GENIE particles won't have a valid G4ID
          //         therefore the comparsion will only look at 
          //         the module of the 3-momentum
          srTruePart = &truthMatch->GetTrueParticle(sr, srTrueInt, srPartCmp, true, false);
          isPrim = true;
        }
        catch ( std::runtime_error& err )
        {
          // guess if it wasn't a primary, it must be a secondary :(
          srTruePart = &truthMatch->GetTrueParticle(sr, srTrueInt, srPartCmp, false, true);
        }

        // however this will fill in any other fields that weren't copied from a GENIE record
        // (which also handles the case where this particle is a secondary)
        FillTruth(*srTruePart, iMCParticleMuID);

        // the particle idx is within the GENIE vector, which may not be the same as the index in the vector here
        // first find the interaction that it goes with
        LOG.VERBOSE() << "      this particle is " << (isPrim ? "PRIMARY" : "SECONDARY") << "\n";
        std::vector<caf::SRTrueParticle> & collection = (isPrim)
                                                        ? srTrueInt.prim
                                                        : srTrueInt.sec;
        std::size_t truthVecIdx = std::distance(collection.begin(),
                                                std::find_if(collection.begin(),
                                                              collection.end(),
                                                              srPartCmp));

        LOG.VERBOSE() << "      index of SRParticle in the SRInteraction " << truthVecIdx << "\n";

        // F says: currently our reco only gives truth matching
        //         to one MCParticle, but we can change that if
        //         we need to...
        cluster.truth.push_back(caf::TrueParticleID{srTrueIntIdx,
                                          (isPrim) ? caf::TrueParticleID::PartType::kPrimary :  caf::TrueParticleID::PartType::kSecondary,
                                          static_cast<int>(truthVecIdx)});
        
        cluster.truthOverlap.push_back(fMuIDClusterMCfrac->at(iMuID));

      }

        sr.nd.gar.ixn[0].muidclusters.push_back(cluster);

    }

  }

  // ------------------------------------------------------------------------------
  void NDGArRecoBranchFiller::FillAssociations(caf::StandardRecord &sr) const
  {
    // F says: sr.nd.gar.ixn already has the right size, so
    //         no need to worry about that again here

    size_t n_associations = fRecoParticleID->size();
    LOG.VERBOSE() << "    GArSoft number of associations: " << n_associations << ".\n";
    sr.nd.gar.ixn[0].nassociations = n_associations;

    // Get the number of associations in the event
    size_t n_assns_tracks = fParticleTrackAssn_ParticleID->size();
    size_t n_assns_ecal   = fParticleECALAssn_ParticleID->size();
    size_t n_assns_muid   = fParticleMuIDAssn_ParticleID->size();

    caf::SRGArAssociation association;

    for (size_t iAssociation=0; iAssociation<n_associations; iAssociation++){

      association.garsoft_assn_id = fRecoParticleID->at(iAssociation);

      association.ECAL_total_energy = fRecoParticleTotalECALEnergy->at(iAssociation);
      association.ECAL_n_hits       = fRecoParticleNHitsECAL->at(iAssociation);
      association.MuID_total_energy = fRecoParticleTotalMuIDEnergy->at(iAssociation);
      association.MuID_n_hits       = fRecoParticleNHitsMuID->at(iAssociation);

      association.ToF_time = fRecoParticleToFTime->at(iAssociation);
      association.ToF_beta = fRecoParticleToFBeta->at(iAssociation);

      association.charge = fRecoParticleCharge->at(iAssociation);

      // Filling associations to other reco objects
      for (size_t iAssnTrack=0; iAssnTrack<n_assns_tracks; ++iAssnTrack){
        LOG.VERBOSE() << "            iAssnTrack in loop: " << iAssnTrack
                      << "            ParticleID:         " << fParticleTrackAssn_ParticleID->at(iAssnTrack)
                      << "            TrackID:            " << fParticleTrackAssn_TrackID->at(iAssnTrack) << ".\n";
        if (association.garsoft_assn_id == fParticleTrackAssn_ParticleID->at(iAssnTrack)){
          LOG.VERBOSE() << "                association found!" << ".\n";
          association.garsoft_trk_assn = fParticleTrackAssn_TrackID->at(iAssnTrack);
        }
      }

      for (size_t iAssnECAL=0; iAssnECAL<n_assns_ecal; ++iAssnECAL){
        LOG.VERBOSE() << "            iAssnECAL in loop:  " << iAssnECAL
                      << "            ParticleID:         " << fParticleECALAssn_ParticleID->at(iAssnECAL)
                      << "            ECalID:             " << fParticleECALAssn_ClusterID->at(iAssnECAL) << ".\n";
        if (association.garsoft_assn_id == fParticleECALAssn_ParticleID->at(iAssnECAL)){
          LOG.VERBOSE() << "                association found!" << ".\n";
          association.garsoft_ecal_assns.push_back(fParticleECALAssn_ClusterID->at(iAssnECAL));
        }
      }

      for (size_t iAssnMuID=0; iAssnMuID<n_assns_muid; ++iAssnMuID){
        LOG.VERBOSE() << "            iAssnMuID in loop:  " << iAssnMuID
                      << "            ParticleID:         " << fParticleMuIDAssn_ParticleID->at(iAssnMuID)
                      << "            MuIDID:             " << fParticleMuIDAssn_ClusterID->at(iAssnMuID) << ".\n";
        if (association.garsoft_assn_id == fParticleMuIDAssn_ParticleID->at(iAssnMuID)){
          LOG.VERBOSE() << "                association found!" << ".\n";
          association.garsoft_muid_assns.push_back(fParticleMuIDAssn_ClusterID->at(iAssnMuID));
        }
      }

      sr.nd.gar.ixn[0].associations.push_back(association);

    }

  }

  // ------------------------------------------------------------------------------
  std::deque<Trigger> NDGArRecoBranchFiller::GetTriggers(int triggerType) const
  {

    size_t numEntries = NDGArRecoTree->GetEntries();
    std::deque<Trigger> triggers;
    if (fTriggers.empty())
    {
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "\n";
      fTriggers.reserve(numEntries);
      for (size_t i = 0; i < numEntries; i++)
      {
        const int placeholderTriggerType = 0;
        // fixme: this check needs to be fixed when we have trigger type info
        if (triggerType >= 0 && triggerType != placeholderTriggerType)
        {
          LOG.VERBOSE() << "    skipping this event" << "\n";
          continue;
        }

        NDGArRecoTree->GetEntry(i);

        fTriggers.emplace_back();
        Trigger & trig = fTriggers.back();
        trig.evtID = fEvent;

        // todo: these are placeholder values until we can propagate enough info through the reco files
        trig.triggerType = 0;
        trig.triggerTime_s = fEvent;
        trig.triggerTime_ns = 0.;

        triggers.push_back(trig);

        LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
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
}
