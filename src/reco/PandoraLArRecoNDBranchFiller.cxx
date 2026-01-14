#include "PandoraLArRecoNDBranchFiller.h"

#include "Params.h"
#include <cmath>
#include <limits>

namespace cafmaker
{

  PandoraLArRecoNDBranchFiller::PandoraLArRecoNDBranchFiller(const std::string &pandoraLArRecoNDFilename,
                                                             const float LArDensity)
      : IRecoBranchFiller("PandoraLArRecoND"),
        m_Triggers(),
        m_LastTriggerReqd(m_Triggers.end()),
        m_LArDensity(LArDensity)
  {
    // Open Pandora LArRecoND hierarchy analysis ROOT file
    m_LArRecoNDFile.reset(TFile::Open(pandoraLArRecoNDFilename.c_str(), "read"));

    if (m_LArRecoNDFile && m_LArRecoNDFile->IsOpen())
    {

      LOG.VERBOSE() << " Using PandoraLArRecoND file " << pandoraLArRecoNDFilename << "\n";

      // Input tree
      m_LArRecoNDTree.reset(dynamic_cast<TTree *>(m_LArRecoNDFile->Get("LArRecoND")));
      if (!m_LArRecoNDTree)
      {
        LOG.FATAL() << "Did not find LArRecoND tree in input file " << pandoraLArRecoNDFilename << "\n";
        throw;
      }

      // Set branch addresses
      m_LArRecoNDTree->SetBranchAddress("event", &m_eventId);
      m_LArRecoNDTree->SetBranchAddress("run", &m_run);
      m_LArRecoNDTree->SetBranchAddress("subRun", &m_subRun);
      m_LArRecoNDTree->SetBranchAddress("unixTime", &m_unixTime);
      m_LArRecoNDTree->SetBranchAddress("unixTimeUsec", &m_unixTimeUsec);
      m_LArRecoNDTree->SetBranchAddress("startTime", &m_startTime);
      m_LArRecoNDTree->SetBranchAddress("triggers", &m_triggerType);
      m_LArRecoNDTree->SetBranchAddress("isShower", &m_isShowerVect);
      m_LArRecoNDTree->SetBranchAddress("sliceId", &m_sliceIdVect);
      m_LArRecoNDTree->SetBranchAddress("startX", &m_startXVect);
      m_LArRecoNDTree->SetBranchAddress("startY", &m_startYVect);
      m_LArRecoNDTree->SetBranchAddress("startZ", &m_startZVect);
      m_LArRecoNDTree->SetBranchAddress("endX", &m_endXVect);
      m_LArRecoNDTree->SetBranchAddress("endY", &m_endYVect);
      m_LArRecoNDTree->SetBranchAddress("endZ", &m_endZVect);
      m_LArRecoNDTree->SetBranchAddress("dirX", &m_dirXVect);
      m_LArRecoNDTree->SetBranchAddress("dirY", &m_dirYVect);
      m_LArRecoNDTree->SetBranchAddress("dirZ", &m_dirZVect);
      m_LArRecoNDTree->SetBranchAddress("energy", &m_energyVect);
      m_LArRecoNDTree->SetBranchAddress("n3DHits", &m_n3DHitsVect);
      m_LArRecoNDTree->SetBranchAddress("mcNuId", &m_mcNuIdVect);
      m_LArRecoNDTree->SetBranchAddress("mcLocalId", &m_mcLocalIdVect);
      m_LArRecoNDTree->SetBranchAddress("isPrimary", &m_isPrimaryVect);
      m_LArRecoNDTree->SetBranchAddress("completeness", &m_completenessVect);
      m_LArRecoNDTree->SetBranchAddress("nuVtxX", &m_nuVtxXVect);
      m_LArRecoNDTree->SetBranchAddress("nuVtxY", &m_nuVtxYVect);
      m_LArRecoNDTree->SetBranchAddress("nuVtxZ", &m_nuVtxZVect);
      m_LArRecoNDTree->SetBranchAddress("isRecoPrimary", &m_isRecoPrimaryVect);
      m_LArRecoNDTree->SetBranchAddress("recoPDG", &m_recoPDGVect);
      // TRACK VARIABLES (PANDORA OUTERFACE)
      m_LArRecoNDTree->SetBranchAddress("trackScore", &m_trackScoreVect);
      m_LArRecoNDTree->SetBranchAddress("trkfitPID_Mu", &m_trkfitPID_Mu);
      m_LArRecoNDTree->SetBranchAddress("trkfitPID_Pro", &m_trkfitPID_Pro);
      m_LArRecoNDTree->SetBranchAddress("trkfitPID_NDF", &m_trkfitPID_NDF);
      m_LArRecoNDTree->SetBranchAddress("trkfitContained", &m_trkfitContained);
      m_LArRecoNDTree->SetBranchAddress("trkfitWallDist", &m_trkfitWallDist);
      m_LArRecoNDTree->SetBranchAddress("trkfitLength", &m_trkfitLength);
      m_LArRecoNDTree->SetBranchAddress("trkfitKEFromLengthMuon", &m_trkfitKEFromLengthMuon);
      m_LArRecoNDTree->SetBranchAddress("trkfitKEFromLengthProton", &m_trkfitKEFromLengthProton);
      m_LArRecoNDTree->SetBranchAddress("trkfitPFromLengthMuon", &m_trkfitPFromLengthMuon);
      m_LArRecoNDTree->SetBranchAddress("trkfitPFromLengthProton", &m_trkfitPFromLengthProton);
      m_LArRecoNDTree->SetBranchAddress("trkfitStartX", &m_trkfitStartX);
      m_LArRecoNDTree->SetBranchAddress("trkfitStartY", &m_trkfitStartY);
      m_LArRecoNDTree->SetBranchAddress("trkfitStartZ", &m_trkfitStartZ);
      m_LArRecoNDTree->SetBranchAddress("trkfitEndX", &m_trkfitEndX);
      m_LArRecoNDTree->SetBranchAddress("trkfitEndY", &m_trkfitEndY);
      m_LArRecoNDTree->SetBranchAddress("trkfitEndZ", &m_trkfitEndZ);
      m_LArRecoNDTree->SetBranchAddress("trkfitStartDirX", &m_trkfitStartDirX);
      m_LArRecoNDTree->SetBranchAddress("trkfitStartDirY", &m_trkfitStartDirY);
      m_LArRecoNDTree->SetBranchAddress("trkfitStartDirZ", &m_trkfitStartDirZ);
      m_LArRecoNDTree->SetBranchAddress("trkfitEndDirX", &m_trkfitEndDirX);
      m_LArRecoNDTree->SetBranchAddress("trkfitEndDirY", &m_trkfitEndDirY);
      m_LArRecoNDTree->SetBranchAddress("trkfitEndDirZ", &m_trkfitEndDirZ);
      // track fit calo
      m_LArRecoNDTree->SetBranchAddress("trkfitTrackCaloE", &m_trkfitTrackCaloE);
      m_LArRecoNDTree->SetBranchAddress("trkfitVisE", &m_trkfitVisE);
      m_LArRecoNDTree->SetBranchAddress("trkfitSliceId", &m_trkfitSliceId);
      m_LArRecoNDTree->SetBranchAddress("trkfitPfoId", &m_trkfitPfoId);
      m_LArRecoNDTree->SetBranchAddress("trkfitX", &m_trkfitX);
      m_LArRecoNDTree->SetBranchAddress("trkfitY", &m_trkfitY);
      m_LArRecoNDTree->SetBranchAddress("trkfitZ", &m_trkfitZ);
      m_LArRecoNDTree->SetBranchAddress("trkfitQ", &m_trkfitQ);
      m_LArRecoNDTree->SetBranchAddress("trkfitRR", &m_trkfitRR);
      m_LArRecoNDTree->SetBranchAddress("trkfitdx", &m_trkfitdx);
      m_LArRecoNDTree->SetBranchAddress("trkfitdQdx", &m_trkfitdQdx);
      m_LArRecoNDTree->SetBranchAddress("trkfitdEdx", &m_trkfitdEdx);
      // SHOWER VARIABLES (PANDORA OUTERFACE)
      m_LArRecoNDTree->SetBranchAddress("shwrfitLength", &m_shwrfitLength);
      m_LArRecoNDTree->SetBranchAddress("shwrfitCentroidX", &m_shwrfitCentroidX);
      m_LArRecoNDTree->SetBranchAddress("shwrfitCentroidY", &m_shwrfitCentroidY);
      m_LArRecoNDTree->SetBranchAddress("shwrfitCentroidZ", &m_shwrfitCentroidZ);
      m_LArRecoNDTree->SetBranchAddress("shwrfitStartX", &m_shwrfitStartX);
      m_LArRecoNDTree->SetBranchAddress("shwrfitStartY", &m_shwrfitStartY);
      m_LArRecoNDTree->SetBranchAddress("shwrfitStartZ", &m_shwrfitStartZ);
      m_LArRecoNDTree->SetBranchAddress("shwrfitDirX", &m_shwrfitDirX);
      m_LArRecoNDTree->SetBranchAddress("shwrfitDirY", &m_shwrfitDirY);
      m_LArRecoNDTree->SetBranchAddress("shwrfitDirZ", &m_shwrfitDirZ);
      m_LArRecoNDTree->SetBranchAddress("shwrSliceId", &m_shwrSliceId);
      m_LArRecoNDTree->SetBranchAddress("shwrClusterId", &m_shwrClusterId);
      m_LArRecoNDTree->SetBranchAddress("shwrdEdx", &m_shwrdEdx);
      m_LArRecoNDTree->SetBranchAddress("shwrEnergy", &m_shwrEnergy);
      
      // We have setup the input tree
      SetConfigured(true);
    }
  }

  // Copy all of the Pandora LArRecoND info to the PandoraLArRecoND branch of the StandardRecord object
  void PandoraLArRecoNDBranchFiller::_FillRecoBranches(const Trigger &trigger,
                                                       caf::StandardRecord &sr,
                                                       const cafmaker::Params &par,
                                                       const TruthMatcher *truthMatch) const
  {
    // Figure out where in our list of triggers this event index is.
    // We should always be looking forwards, since we expect to be traversing in that direction
    auto it_start = (m_LastTriggerReqd == m_Triggers.end()) ? m_Triggers.cbegin() : m_LastTriggerReqd;
    auto itTrig = std::find(it_start, m_Triggers.cend(), trigger);
    if (itTrig == m_Triggers.end())
    {
      LOG.FATAL() << " Reco branch filler '" << GetName() << "' could not find trigger with evtID == "
                  << trigger.evtID << "!  Abort.\n";
      abort();
    }
    std::size_t idx = std::distance(m_Triggers.cbegin(), itTrig);

    LOG.VERBOSE() << " Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID
                  << ", internal evt idx = " << idx << ".\n";

    // Get the event entry
    m_LArRecoNDTree->GetEntry(fEntryMap[idx]);

    // Set the event and run numbers
    sr.meta.nd_lar.enabled = true;
    sr.meta.nd_lar.event = m_eventId;
    sr.meta.nd_lar.run = m_run;
    sr.meta.nd_lar.subrun = m_subRun;

    // Number of PFO clusters (tracks or showers) = size of the ntuple vectors for this event
    const int nClusters = (m_sliceIdVect != nullptr) ? m_sliceIdVect->size() : 0;

    // Find the unique neutrino slices and their vertices (1 slice = 1 neutrino PFO).
    // Start to build the neutrino interaction (vertex) objects as well
    std::vector<int> uniqueSliceIDs;
    std::vector<caf::SRInteraction> nuInteractions;

    int oldSliceId{-1};
    for (int i = 0; i < nClusters; i++)
    {
      const int sliceId = (m_sliceIdVect != nullptr) ? (*m_sliceIdVect)[i] : -1;
      if (sliceId != oldSliceId)
      {
        uniqueSliceIDs.emplace_back(sliceId);

        // Create standard record interaction object
        caf::SRInteraction interaction;

        // Id
        interaction.id = sliceId;

        // Reconstructed neutrino interaction vertex
        const float vtxX = (m_nuVtxXVect != nullptr) ? (*m_nuVtxXVect)[i] : 0.0;
        const float vtxY = (m_nuVtxYVect != nullptr) ? (*m_nuVtxYVect)[i] : 0.0;
        const float vtxZ = (m_nuVtxZVect != nullptr) ? (*m_nuVtxZVect)[i] : 0.0;
        interaction.vtx = caf::SRVector3D(vtxX, vtxY, vtxZ);

        // Initialise total neutrino energy to zero
        interaction.Enu.calo = 0.0;

        // Add to interaction vector
        nuInteractions.emplace_back(interaction);

        // Keep track of the current sliceId
        oldSliceId = sliceId;
      }
    }

    // Set the size for the Pandora NDLAr standard record for this trigger/event
    const int nNeutrinos = nuInteractions.size();
    sr.nd.lar.pandora.resize(nNeutrinos);
    sr.nd.lar.npandora = sr.nd.lar.pandora.size();

    // Check if input file has Pandora Outerface fields. If false, fill
    // Standard Record classes with default filler
    const bool inputHasOuterfaceBranches = HasOuterfaceBranches();
    LOG.VERBOSE() << "Input file has Pandora Outerface branches? " << inputHasOuterfaceBranches << "\n";

    if (inputHasOuterfaceBranches)
    {
      // Fill track and shower info. Both use the same clusters (PFOs), and no
      // distinction is made (yet) to identify which are tracks or showers.
      // This also fills in the interaction common variables
      FillRecoParticles(sr, nClusters, uniqueSliceIDs, nuInteractions, truthMatch);
    }
    else
    {
      // TODO : delete this when pandora outerfaces has been merged into LArRecoND main
      FillRecoParticlesDefault(sr, nClusters, uniqueSliceIDs, nuInteractions, truthMatch);
    }


    // Store the common neutrino interactions
    sr.common.ixn.pandora.reserve(nNeutrinos);
    sr.common.ixn.npandora = nNeutrinos;
    for (auto interaction : nuInteractions)
      sr.common.ixn.pandora.emplace_back(std::move(interaction));
  }

  bool PandoraLArRecoNDBranchFiller::HasOuterfaceBranches() const
  { // should come up with something more robust but ok for now? (to be deleted in the future anyway)
    if (m_trkfitPID_Mu == nullptr &&
        m_trkfitPID_Pro == nullptr &&
        m_trkfitPID_NDF == nullptr &&
        m_shwrfitLength == nullptr &&
        m_shwrfitStartX == nullptr) return false;

    return true;
  }

  // ------------------------------------------------------------------------------
  bool PandoraLArRecoNDBranchFiller::FillTrack(const int i, caf::SRRecoParticle& recoParticle) const
  {  // Assign the pdg for the fitting hypotesis that has the lowest chi2
     // as version 0, the track is chosen between muon and proton
     // TODO consider also pion and kaons
     
     const int assigned_pdg = ((*m_trkfitPID_Mu)[i] < (*m_trkfitPID_Pro)[i]) ? m_muonPDG : m_protonPDG;
     recoParticle.pdg = assigned_pdg;

     const float trkfitStartX = (*m_trkfitStartX)[i];
     const float trkfitStartY = (*m_trkfitStartY)[i];
     const float trkfitStartZ = (*m_trkfitStartZ)[i];
     const float trkfitEndX = (*m_trkfitEndX)[i];
     const float trkfitEndY = (*m_trkfitEndY)[i];
     const float trkfitEndZ = (*m_trkfitEndZ)[i];
     const float trkfitStartDirX = (*m_trkfitStartDirX)[i];
     const float trkfitStartDirY = (*m_trkfitStartDirY)[i];
     const float trkfitStartDirZ = (*m_trkfitStartDirZ)[i];
     const float trkfitStartDirMag = sqrt(trkfitStartDirX * trkfitStartDirX + trkfitStartDirY * trkfitStartDirY + trkfitStartDirZ * trkfitStartDirZ); 
     
     float p_mod = 1.; 
     
     recoParticle.E_method = caf::PartEMethod::kRange; 
     
     if(m_muonPDG == assigned_pdg)
     {
       recoParticle.E = (*m_trkfitKEFromLengthMuon)[i] + m_mMuon; // [GeV]
       p_mod = (*m_trkfitPFromLengthMuon)[i];
       recoParticle.score = (*m_trkfitPID_Mu)[i];
     }
     else if (m_protonPDG == assigned_pdg)
     {
       recoParticle.E = (*m_trkfitKEFromLengthProton)[i] + m_mProton; // [GeV]
       p_mod = (*m_trkfitPFromLengthProton)[i];
       recoParticle.score = (*m_trkfitPID_Pro)[i];
     }
     else  
     {
       // TODO: to be filled with other particle options once available
     } 
    
     caf::SRVector3D start{trkfitStartX, trkfitStartY, trkfitStartZ};
     caf::SRVector3D end{trkfitEndX, trkfitEndY, trkfitEndZ};
     caf::SRVector3D p{trkfitStartDirX / trkfitStartDirMag * p_mod, trkfitStartDirY / trkfitStartDirMag* p_mod, trkfitStartDirZ / trkfitStartDirMag* p_mod};

     recoParticle.start = start;
     recoParticle.end = end;
     recoParticle.p = p;
     recoParticle.walldist = (*m_trkfitWallDist)[i];
     recoParticle.origRecoObjType = caf::RecoObjType::kTrack; 

     return true;

 }
 // ------------------------------------------------------------------------------
  bool PandoraLArRecoNDBranchFiller::FillShower(const int i, caf::SRRecoParticle& recoParticle) const
  {
     const float shwrfitStartX = (*m_shwrfitStartX)[i];
     const float shwrfitStartY = (*m_shwrfitStartY)[i];
     const float shwrfitStartZ = (*m_shwrfitStartZ)[i];
     const float vertexX = (m_nuVtxXVect != nullptr) ? (*m_nuVtxXVect)[i] : 0.;
     const float vertexY = (m_nuVtxYVect != nullptr) ? (*m_nuVtxYVect)[i] : 0.;
     const float vertexZ = (m_nuVtxZVect != nullptr) ? (*m_nuVtxZVect)[i] : 0.;
     const float conversionGap = (m_nuVtxXVect != nullptr) ?  sqrt(std::pow((shwrfitStartX - vertexX),2) + std::pow((shwrfitStartY - vertexY),2) + std::pow((shwrfitStartZ - vertexZ),2)) : -999;
     const float dEdx = (m_shwrdEdx != nullptr) ? (*m_shwrdEdx)[i] : -999.;
      
     if ((*m_trackScoreVect)[i] >= m_TrackShowerCut) // track w/ failed trackfit: we still want to fill the recoParticle with the shower info (at least we have something)
     {
       recoParticle.pdg = m_antiprotonPDG; // still a track --> assign as default antiproton
       recoParticle.origRecoObjType = caf::RecoObjType::kTrack;
       recoParticle.E = (*m_shwrEnergy)[i]; // TODO what rest mass to assign?
     }
     else
     {
       recoParticle.origRecoObjType = caf::RecoObjType::kShower; 
       if ((conversionGap > m_ConversionGapCut) && (dEdx > m_dEdxShowerCut) ) // TODO (or even delete? ok for v0). Some considerations here: the final decision on gamma vs e- should be left to the analyzer, based on the information propagated to the CAF output. We shouldn't discard any hyphotesis beforehand. Future development will also depend on how we restructure the Standard Reco track,shower,particle classes.Last but not leas, if for any reason vertex is null, m_mElectron is assigned.
       {
         recoParticle.pdg = m_gammaPDG;
         recoParticle.E = (*m_shwrEnergy)[i];
       }
       else
       {
         recoParticle.pdg = m_electronPDG;
         recoParticle.E = (*m_shwrEnergy)[i] + m_mElectron * 1e-3; // TODO: check if default units are MeV for this branch.
       }
     }
        
     float shwrfitDirX = (*m_shwrfitDirX)[i];
     float shwrfitDirY = (*m_shwrfitDirY)[i];
     float shwrfitDirZ = (*m_shwrfitDirZ)[i];
     float shwrLength = (*m_shwrfitLength)[i];
     float shwrfitEndX = shwrfitStartX + shwrfitDirX * shwrLength;
     float shwrfitEndY = shwrfitStartY + shwrfitDirY * shwrLength;
     float shwrfitEndZ = shwrfitStartZ + shwrfitDirZ * shwrLength;
      
     caf::SRVector3D start{shwrfitStartX, shwrfitStartY, shwrfitStartZ};
     caf::SRVector3D end{shwrfitEndX, shwrfitEndY, shwrfitEndZ};
     caf::SRVector3D dir{shwrfitDirX, shwrfitDirY, shwrfitDirZ};
     
     recoParticle.E_method = caf::PartEMethod::kRange;
     recoParticle.start = start;    
     recoParticle.end = end;    
     recoParticle.p = dir * recoParticle.E; 
     recoParticle.contained = false; // TODO read from dedicated branch
     recoParticle.walldist = -999.; // TODO read from dedicated branch

     return true;
  }

// ------------------------------------------------------------------------------

 void PandoraLArRecoNDBranchFiller::FillTruthInfo(const unsigned i, const TruthMatcher *truthMatch, caf::StandardRecord &sr, caf::TrueParticleID& truePartID) const
 {
   const long mcNuId = (m_mcNuIdVect != nullptr) ? (*m_mcNuIdVect)[i] : 0;
   const long mcId = (m_mcLocalIdVect != nullptr) ? (*m_mcLocalIdVect)[i] : 0;
   const int isPrimary = (m_isPrimaryVect != nullptr) ? (*m_isPrimaryVect)[i] : -1;

   if (isPrimary == 1)
   {
     truePartID.type = caf::TrueParticleID::kPrimary;
     truePartID.part = mcId;
   }
   else if (isPrimary == -1)
   {
     truePartID.type = caf::TrueParticleID::kUnknown;
   }
   else
   {
     truePartID.type = caf::TrueParticleID::kSecondary;
   }

   if (mcNuId != 0)
   {
     caf::SRTrueInteraction &srTrueInt = truthMatch->GetTrueInteraction(sr, mcNuId);
     const auto predicate = [&srTrueInt](const caf::SRTrueInteraction &ixn)
     { return ixn.id == srTrueInt.id; };
     // Get the truth interaction index
     const int srTrueIntIdx = std::distance(sr.mc.nu.begin(), std::find_if(sr.mc.nu.begin(), sr.mc.nu.end(), predicate));
     truePartID.ixn = srTrueIntIdx;

      // If the particle is not a primary, we might want to create a new particle if it wasn't created originally
      if (isPrimary != 1)
      {
        const auto pred = [&mcId](const caf::SRTrueParticle &part)
        { return part.G4ID == mcId; };
        truePartID.part = std::distance(srTrueInt.sec.begin(),
                                    std::find_if(srTrueInt.sec.begin(), srTrueInt.sec.end(), pred));
      }
   }
 }

// ------------------------------------------------------------------------------
 void PandoraLArRecoNDBranchFiller::FillRecoParticlesDefault(caf::StandardRecord &sr, const int nClusters,
                                               const std::vector<int> &uniqueSliceIDs,
                                                std::vector<caf::SRInteraction> &nuInteractions,
                                                const TruthMatcher *truthMatch) const
 {
   LOG.DEBUG() << " Using default Pandora Reco ND branch filler \n";
   LOG.VERBOSE() << " Pandora LArRecoND FillTracks using " << nClusters << " PFO clusters\n";

    const caf::TrueParticleID nullTrueID;
    // Direction of the longest track
    float maxTrackLength{0.0};
    caf::SRVector3D longestTrackDir; 

    for (int i = 0; i < nClusters; i++)
    {
      // Check that the PFO is a track and not a shower
      const int isShower = (m_isShowerVect != nullptr) ? (*m_isShowerVect)[i] : 0;

      // Starting position (vertex or first hit location)
      const float startX = (m_startXVect != nullptr) ? (*m_startXVect)[i] : 0.0;
      const float startY = (m_startYVect != nullptr) ? (*m_startYVect)[i] : 0.0;
      const float startZ = (m_startZVect != nullptr) ? (*m_startZVect)[i] : 0.0;
      const caf::SRVector3D start{startX, startY, startZ};

      // End position
      const float endX = (m_endXVect != nullptr) ? (*m_endXVect)[i] : 0.0;
      const float endY = (m_endYVect != nullptr) ? (*m_endYVect)[i] : 0.0;
      const float endZ = (m_endZVect != nullptr) ? (*m_endZVect)[i] : 0.0;
      const caf::SRVector3D end{endX, endY, endZ};

      // Principal axis direction
      const float dirX = (m_dirXVect != nullptr) ? (*m_dirXVect)[i] : 0.0;
      const float dirY = (m_dirYVect != nullptr) ? (*m_dirYVect)[i] : 0.0;
      const float dirZ = (m_dirZVect != nullptr) ? (*m_dirZVect)[i] : 0.0;
      const caf::SRVector3D dir{dirX, dirY, dirZ};

      // Energy (GeV)
      const float energy = (m_energyVect != nullptr) ? (*m_energyVect)[i] : 0.0;

      // Total number of 3D hits in the cluster
      const int n3DHits = (m_n3DHitsVect != nullptr) ? (*m_n3DHitsVect)[i] : 0;

      // Slice id of the PFO cluster
      const int sliceId = (m_sliceIdVect != nullptr) ? (*m_sliceIdVect)[i] : 0;

      const int isRecoPrimary = (m_isRecoPrimaryVect != nullptr) ? (*m_isRecoPrimaryVect)[i] : 0;

      const float completeness = (m_completenessVect != nullptr) ? (*m_completenessVect)[i] : 0.0;
      std::vector<float> truthOverlap;
      truthOverlap.emplace_back(completeness);

      // Neutrino index number 0 to N-1 for N neutrinos
      const int nuIndex = std::distance(uniqueSliceIDs.begin(), std::find(uniqueSliceIDs.begin(),
                                                                          uniqueSliceIDs.end(), sliceId));
      std::vector<caf::TrueParticleID> truePartIDVect;
      caf::TrueParticleID truePartID;
      FillTruthInfo(i, truthMatch, sr, truePartID);
      truePartIDVect.emplace_back(truePartID);

      if (isShower) // fill SR shower
      {
        caf::SRShower shower;

        shower.start = start;
        shower.direction = dir;
        shower.Evis = energy;
        shower.truth = truePartIDVect;
        shower.truthOverlap = truthOverlap;
        
        sr.nd.lar.pandora[nuIndex].showers.emplace_back(std::move(shower));
        sr.nd.lar.pandora[nuIndex].nshowers++;
      }
      else // fill SR track
      {
        caf::SRTrack track;

        track.start = start;
        track.end = end;
        track.dir = dir;
        track.enddir = dir;
        track.Evis = energy;
        track.E = energy;
        track.qual = n3DHits * 1.0;
        track.truth = truePartIDVect;
        track.truthOverlap = truthOverlap;

        // Cluster length from start and end points
        const float dX = endX - startX;
        const float dY = endY - startY;
        const float dZ = endZ - startZ;
        track.len_cm = sqrt(dX * dX + dY * dY + dZ * dZ);
        track.len_gcm2 = track.len_cm * m_LArDensity;

        sr.nd.lar.pandora[nuIndex].tracks.emplace_back(std::move(track));
        sr.nd.lar.pandora[nuIndex].ntracks++;
      }
      
      // now fill the SR reco particle
      caf::SRRecoParticle recoParticle;

      recoParticle.E = energy;
      recoParticle.E_method = caf::PartEMethod::kCalorimetry;
      recoParticle.p = energy * dir;
      recoParticle.start = start;
      recoParticle.end = end;
      recoParticle.truth = truePartIDVect;
      recoParticle.truthOverlap = truthOverlap;
      recoParticle.primary = (isRecoPrimary == 1) ? true : false;
      recoParticle.pdg = (m_recoPDGVect != nullptr) ? (*m_recoPDGVect)[i] : 0;
      
      // Add particle to the interaction
      caf::SRInteraction &interaction = nuInteractions[nuIndex];
      interaction.part.pandora.emplace_back(std::move(recoParticle));
      interaction.part.npandora++;

      // Add track truth info
      const caf::TrueParticleID nullTrueID;
      const caf::TrueParticleID trackTrueID = (truePartIDVect.size() > 0) ? truePartIDVect[0] : nullTrueID;
      const float trackOverlap = (truthOverlap.size() > 0) ? truthOverlap[0] : 0.0;
      const int trackIxn = trackTrueID.ixn;

      interaction.truth.emplace_back(trackIxn);
      interaction.truthOverlap.emplace_back(trackOverlap);

      // Update total interaction neutrino energy
      interaction.Enu.calo += energy;

    }
 }
 // ------------------------------------------------------------------------------
 void PandoraLArRecoNDBranchFiller::FillRecoParticles(caf::StandardRecord &sr, const int nClusters,
                                               const std::vector<int> &uniqueSliceIDs,
                                                std::vector<caf::SRInteraction> &nuInteractions,
                                                const TruthMatcher *truthMatch) const
  {
    // Create tracks for each PFO (cluster) in the event
    LOG.VERBOSE() << " Pandora LArRecoND FillTracks using " << nClusters << " PFO clusters\n";
 
    const caf::TrueParticleID nullTrueID;

    // Value and direction of the longest track, value and direction of the most energetic shower
    
    float longestTrack{0.0};
    float maxShowerE{0.0};
    caf::SRVector3D longestTrackDir; // TODO fill it in the loop below
    caf::SRVector3D maxShowerEDir;

    for (int i = 0; i < nClusters; i++)
    {
      // Use truth matching info from Pandora's LArContent hierarchy tools.
      // For LArRecoND MC SpacePoints, we use the unique file_traj_id's
      // for the MCParticles and the neutrino ID is the unique vertex_id.
      // We also store the local traj_id's for the CAF truth matching tools
      
      // Neutrino ID = vertex_id
      const long mcNuId = (m_mcNuIdVect != nullptr) ? (*m_mcNuIdVect)[i] : 0;
      // MC particle ID = traj_id
      const long mcId = (m_mcLocalIdVect != nullptr) ? (*m_mcLocalIdVect)[i] : 0;
      const int isPrimary = (m_isPrimaryVect != nullptr) ? (*m_isPrimaryVect)[i] : -1;

      caf::TrueParticleID truePartID;

      if (isPrimary == 1)
      {
        truePartID.type = caf::TrueParticleID::kPrimary;
        truePartID.part = mcId;
      }
      else if (isPrimary == -1)
      {
        truePartID.type = caf::TrueParticleID::kUnknown;
      }
      else
      {
        truePartID.type = caf::TrueParticleID::kSecondary;
      }

      if (mcNuId != 0)
      {
        // Get the true interaction in the stack
        caf::SRTrueInteraction &srTrueInt = truthMatch->GetTrueInteraction(sr, mcNuId);
        const auto predicate = [&srTrueInt](const caf::SRTrueInteraction &ixn)
        { return ixn.id == srTrueInt.id; };
        // Get the truth interaction index
        const int srTrueIntIdx = std::distance(sr.mc.nu.begin(), std::find_if(sr.mc.nu.begin(), sr.mc.nu.end(), predicate));
        truePartID.ixn = srTrueIntIdx;

        // If the particle is not a primary, we might want to create a new particle if it wasn't created originally
        if (isPrimary != 1)
        {
          const auto pred = [&mcId](const caf::SRTrueParticle &part)
          { return part.G4ID == mcId; };
          truePartID.part = std::distance(srTrueInt.sec.begin(),
                                          std::find_if(srTrueInt.sec.begin(), srTrueInt.sec.end(), pred));
        }
      }

      // Just store the best MC match
      std::vector<caf::TrueParticleID> truePartIDVect;
      truePartIDVect.emplace_back(truePartID);
      
      // Fraction of true MC hits that are captured by the reconstructed cluster
      const float completeness = (m_completenessVect != nullptr) ? (*m_completenessVect)[i] : 0.0;
      std::vector<float> truthOverlap;
      truthOverlap.emplace_back(completeness);

      // Add track to the record
      // Slice id of the PFO cluster
      const int sliceId = (m_sliceIdVect != nullptr) ? (*m_sliceIdVect)[i] : 0;
      // Neutrino index number 0 to N-1 for N neutrinos
      const int nuIndex = std::distance(uniqueSliceIDs.begin(), std::find(uniqueSliceIDs.begin(),
                                                                          uniqueSliceIDs.end(), sliceId));

      // Store interaction info
      caf::SRInteraction &interaction = nuInteractions[nuIndex];

      // Track reco particle
      caf::SRRecoParticle recoParticle;
            
      // Truth info
      recoParticle.truth = truePartIDVect;
      recoParticle.truthOverlap = truthOverlap;
      const caf::TrueParticleID trackTrueID = (truePartIDVect.size() > 0) ? truePartIDVect[0] : nullTrueID;
      const int trackIxn = trackTrueID.ixn;
      const float trackOverlap = (truthOverlap.size() > 0) ? truthOverlap[0] : 0.0;

      // Is reco primary & reco PDG hypothesis
      const int isRecoPrimary = (m_isRecoPrimaryVect != nullptr) ? (*m_isRecoPrimaryVect)[i] : 0;
      recoParticle.primary = (isRecoPrimary == 1) ? true : false;

      if(m_trackScoreVect != nullptr)
      {
        recoParticle.contained = (m_trkfitContained != nullptr) ? (*m_trkfitContained)[i] : false;
        
        bool isTrackFitFailed = ((*m_trkfitLength)[i] < std::numeric_limits<float>::epsilon());
        
        if((*m_trackScoreVect)[i] >= m_TrackShowerCut && !isTrackFitFailed) // reco particle is a track and trackfit was succesful
        { 
          FillTrack(i, recoParticle);
          
          // Create standard record track
          caf::SRTrack track;
          track.E = recoParticle.E;
          track.Evis = (m_trkfitVisE != nullptr) ? (*m_trkfitVisE)[i] : -999.;

          // Total number of 3D hits in the cluster
          const int n3DHits = (m_n3DHitsVect != nullptr) ? (*m_n3DHitsVect)[i] : 0;
          track.qual = n3DHits * 1.0;
          track.start = recoParticle.start;
          track.end = recoParticle.end;
          track.len_cm = (*m_trkfitLength)[i];
          track.len_gcm2 = track.len_cm * m_LArDensity;
          track.truth = truePartIDVect; 
          track.truthOverlap = truthOverlap;
          track.dir = recoParticle.p.Unit(); // p vector build from trkfitStartDir
          caf::SRVector3D enddir = {(*m_trkfitEndDirX)[i],(*m_trkfitEndDirY)[i],(*m_trkfitEndDirZ)[i]};
          track.enddir = enddir;
          
          if (track.len_cm > longestTrack)
          {
            longestTrack = track.len_cm;
            longestTrackDir = track.dir;
          }

          // Initialise total neutrino energy to zero
          interaction.Enu.calo = track.Evis;

          sr.nd.lar.pandora[nuIndex].tracks.emplace_back(std::move(track));
          sr.nd.lar.pandora[nuIndex].ntracks++;

        }
        else // reco particle is shower or track with failed trackfit
        {
          FillShower(i, recoParticle);

          // create standard record shower and fill it 
          caf::SRShower shower;
          shower.Evis = (*m_shwrEnergy)[i]; 
          shower.start = recoParticle.start;
          shower.direction = recoParticle.p.Unit(); // p build from shwrfitDir
          shower.truth = truePartIDVect;
          shower.truthOverlap = truthOverlap;

          if (shower.Evis > maxShowerE)
          {
            maxShowerE = shower.Evis;
            maxShowerEDir = shower.direction;
          }
          
          // Update total interaction neutrino energy
          interaction.Enu.calo += shower.Evis;

          sr.nd.lar.pandora[nuIndex].showers.emplace_back(std::move(shower));
          sr.nd.lar.pandora[nuIndex].nshowers++;

        }

        LOG.DEBUG() << "trackScore "      << (*m_trackScoreVect)[i] 
                    << ", is track? "     << ((*m_trackScoreVect)[i] > 0.5)
                    << ", is shower? "    << ((*m_trackScoreVect)[i] < 0.5)
                    << ", NDF "           << (*m_trkfitPID_NDF)[i]
                    << ", muon score "    << (*m_trkfitPID_Mu)[i] 
                    << ", proton score "  << (*m_trkfitPID_Pro)[i] 
                    << ", assigned pdg "  << recoParticle.pdg 
                    << ", contained ?  "  << recoParticle.contained
                    << ", E  "            << recoParticle.E
                    << ", E_method  "     << recoParticle.E_method
                    << ", start = ("      << recoParticle.start.X()
                    << ", "               << recoParticle.start.Y()
                    << ", "               << recoParticle.start.Z()
                    << ")"
                    << ", end = ("        << recoParticle.end.X()
                    << ", "               << recoParticle.end.Y()
                    << ", "               << recoParticle.end.Z()
                    << ")"
                    << ", momentum = ("   << recoParticle.p.X()
                    << ", "               << recoParticle.p.Y()
                    << ", "               << recoParticle.p.Z()
                    << ")"
                    << "\n";
      }
      else // trackfit failed: neither track nor shower
      {
        LOG.DEBUG() << "trackfit failed for particle number : " << i << ", assigning pdg 0\n";
        // TODO : fill with hierarchy analysis info
        recoParticle.pdg = 0;
      }

      // Add particle to the interaction
      interaction.part.pandora.emplace_back(std::move(recoParticle));
      interaction.part.npandora++;

      // Add track truth info
      interaction.truth.emplace_back(trackIxn);
      interaction.truthOverlap.emplace_back(trackOverlap);

      // Update interaction longest track direction
      interaction.dir.lngtrk = longestTrackDir;
      interaction.dir.heshw = maxShowerEDir;

    }
  }

  // ------------------------------------------------------------------------------
  bool PandoraLArRecoNDBranchFiller::IsBeamTrigger(int triggerType) const
  {
    if (triggerType == 5 || triggerType == 1) // io_group = 5 for beam trigger and 1 for MC in Flow
    {
      return true;
    }
    return false;
  }

  // ------------------------------------------------------------------------------
  std::deque<Trigger> PandoraLArRecoNDBranchFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    int iTrigger = 0;

    if (m_Triggers.empty())
    {
      const int nEvents = m_LArRecoNDTree->GetEntries();
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName()
                  << "' from " << nEvents << " Pandora LArRecoND tree entries:\n";

      m_Triggers.reserve(nEvents);
      for (int entry = 0; entry < nEvents; entry++)
      {
        m_LArRecoNDTree->GetEntry(entry);

        if ((triggerType >= 0 && m_triggerType != triggerType) || (beamOnly && !IsBeamTrigger(m_triggerType))) // skip if not the right type
        {
          LOG.VERBOSE() << "    skipping trigger ID=" << m_triggerType << "\n";
          continue;
        }

        fEntryMap[iTrigger] = entry;
        iTrigger += 1;

        m_Triggers.emplace_back();
        Trigger &trig = m_Triggers.back();
        // Event number
        trig.evtID = m_eventId;

        trig.triggerType = m_triggerType;
        // unix_ts trigger time (seconds)
        trig.triggerTime_s = m_unixTime;
        // unix_time_usec ticks (microseconds) converted to nanoseconds
        trig.triggerTime_ns = m_unixTimeUsec * 1000;

        LOG.VERBOSE() << "  added trigger: evtID = " << trig.evtID
                      << ", triggerType = " << trig.triggerType
                      << ", triggerTime_s = " << trig.triggerTime_s
                      << ", triggerTime_ns = " << trig.triggerTime_ns
                      << "\n";
      }
      // Since we just modified the list, any iterators have been invalidated
      m_LastTriggerReqd = m_Triggers.end();
    }

    std::deque<Trigger> triggers;
    for (const Trigger &trigger : m_Triggers)
    {
      if (triggerType < 0 || triggerType == m_Triggers.back().triggerType)
        triggers.push_back(trigger);
    }

    return triggers;
  }

} // end namespace
