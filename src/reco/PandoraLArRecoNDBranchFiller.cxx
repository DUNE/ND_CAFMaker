#include "PandoraLArRecoNDBranchFiller.h"

namespace cafmaker
{

  PandoraLArRecoNDBranchFiller::~PandoraLArRecoNDBranchFiller()
  {
      if (m_LArRecoNDFile && m_LArRecoNDFile->IsOpen())
      {
	  delete m_LArRecoNDTree; m_LArRecoNDTree = nullptr;
      }
      delete m_LArRecoNDFile; m_LArRecoNDFile = nullptr;      
  }

  PandoraLArRecoNDBranchFiller::PandoraLArRecoNDBranchFiller(const std::string &pandoraLArRecoNDFilename)
  : IRecoBranchFiller("PandoraLArRecoND"),
    m_Triggers(),
    m_LastTriggerReqd(m_Triggers.end())
  {
      // Open Pandora LArRecoND hierarchy analysis ROOT file
      m_LArRecoNDFile = TFile::Open(pandoraLArRecoNDFilename.c_str(), "read");
      if (m_LArRecoNDFile && m_LArRecoNDFile->IsOpen()) {

	  LOG.VERBOSE() << " Using PandoraLArRecoND file " << pandoraLArRecoNDFilename << "\n";
	  
	  // Input tree
	  m_LArRecoNDTree = dynamic_cast<TTree*>(m_LArRecoNDFile->Get("LArRecoND"));
	  if (!m_LArRecoNDTree) {
	      LOG.FATAL() << "Did not find LArRecoND tree in input file " << pandoraLArRecoNDFilename << "\n";
	      throw;
	  }

	  // Set branch addresses
	  m_LArRecoNDTree->SetBranchAddress("event", &m_eventId);
	  m_LArRecoNDTree->SetBranchAddress("run", &m_run);
	  m_LArRecoNDTree->SetBranchAddress("subRun", &m_subRun);
	  m_LArRecoNDTree->SetBranchAddress("startTime", &m_startTime);
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
	  m_LArRecoNDTree->SetBranchAddress("isPrimary", &m_isPrimaryVect);
	  m_LArRecoNDTree->SetBranchAddress("mcId", &m_mcIdVect);
	  m_LArRecoNDTree->SetBranchAddress("completeness", &m_completenessVect);
	  
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
    m_LArRecoNDTree->GetEntry(idx);

    // Set the event and run numbers
    sr.meta.nd_lar.enabled = true;
    sr.meta.nd_lar.event = m_eventId;
    sr.meta.nd_lar.run = m_run;
    sr.meta.nd_lar.subrun = m_subRun;

    // Set the size for the Pandora NDLAr standard record for this trigger/event
    const int nClusters = (m_sliceIdVect != nullptr) ? m_sliceIdVect->size() : 0;
    sr.nd.lar.pandora.resize(nClusters);
    sr.nd.lar.npandora = sr.nd.lar.pandora.size();

    // Fill track and shower info. Both use the same clusters (PFOs), and no
    // distinction is made (yet) to identify which are tracks or showers
    FillTracks(sr, nClusters);
    FillShowers(sr, nClusters);
  }

  // ------------------------------------------------------------------------------
    void PandoraLArRecoNDBranchFiller::FillTracks(caf::StandardRecord& sr, const int nClusters) const
  {
    // Create tracks for each PFO (cluster) in the event
    LOG.VERBOSE() << " Pandora LArRecoND FillTracks using " << nClusters <<" PFO clusters\n";
    
    for (int i = 0; i < nClusters; i++)
    {
	// Check that the PFO is a track and not a shower
	const int isShower = (*m_isShowerVect)[i];
	if (isShower == 1)
	    continue;

	// Slice id of the PFO cluster
	const int sliceId = (*m_sliceIdVect)[i];

	// Create standard record track
	caf::SRTrack track;
	// Starting position (vertex or first hit location)
	const float startX = (m_startXVect != nullptr) ? (*m_startXVect)[i] : 0.0;
	const float startY = (m_startYVect != nullptr) ? (*m_startYVect)[i] : 0.0;
	const float startZ = (m_startZVect != nullptr) ? (*m_startZVect)[i] : 0.0;
	track.start = caf::SRVector3D(startX, startY, startZ);

	// End position
	const float endX = (m_endXVect != nullptr) ? (*m_endXVect)[i] : 0.0;
	const float endY = (m_endYVect != nullptr) ? (*m_endYVect)[i] : 0.0;
	const float endZ = (m_endZVect != nullptr) ? (*m_endZVect)[i] : 0.0;
	track.end = caf::SRVector3D(endX, endY, endZ);

	// Principal axis direction
	const float dirX = (m_dirXVect != nullptr) ? (*m_dirXVect)[i] : 0.0;
	const float dirY = (m_dirYVect != nullptr) ? (*m_dirYVect)[i] : 0.0;
	const float dirZ = (m_dirZVect != nullptr) ? (*m_dirZVect)[i] : 0.0;
	track.dir = caf::SRVector3D(dirX, dirY, dirZ);
	track.enddir = caf::SRVector3D(dirX, dirY, dirZ);

	// Energy (GeV)
	const float energy = (m_energyVect != nullptr) ? (*m_energyVect)[i] : 0.0;
	track.Evis = energy;
	track.E = energy;

	// Total number of 3D hits in the cluster
	const int n3DHits = (m_n3DHitsVect != nullptr) ? (*m_n3DHitsVect)[i] : 0;
	track.qual = n3DHits*1.0;

	// Cluster length from start and end points
	const float dX = endX - startX;
	const float dY = endY - startY;
	const float dZ = endZ - startZ;
	track.len_cm = sqrt(dX*dX + dY*dY + dZ*dZ);

	// Cluster length multiplied by LAr density (g/cm2)
	track.len_gcm2 = track.len_cm*m_LArRho;

	// Use truth matching info from Pandora's LArContent hierarchy tools.
	// For LArRecoND MC SpacePoints, we offset the MCId's to make them all unique:
	// nuId = orig_nuId + 10^8, so orig_nuId = nuId - 10^8
	// mcId = orig_mcId + nuIndex*10^6, where nuIndex = 0 to N-1 neutrinos
	// nuIndex = int(mcId/10^6), so orig_mcId = mcId - nuIndex*10^6
	// Use the original MCId values
	const int mcNuId = (m_mcNuIdVect != nullptr) ? (*m_mcNuIdVect)[i] : 0;
	const int origMCNuId = (mcNuId > m_nuIdOffset) ? mcNuId - m_nuIdOffset : mcNuId;
	const int isPrimary = (m_isPrimaryVect != nullptr) ? (*m_isPrimaryVect)[i] : -1;
	const int mcId = (m_mcIdVect != nullptr) ? (*m_mcIdVect)[i] : 0;
	const int nuIndex = int(mcId/m_maxMCId);
	const int origMCId = mcId - nuIndex*m_maxMCId;
	
	caf::TrueParticleID trueID;
	trueID.ixn = origMCNuId;
	if (isPrimary == 1) {
	    trueID.type = caf::TrueParticleID::kPrimary;
	} else if (isPrimary == -1) {
	    trueID.type = caf::TrueParticleID::kUnknown;
	} else {
	    trueID.type = caf::TrueParticleID::kSecondary;
	}
	trueID.part = origMCId;

	// Just store the best MC match
	std::vector<caf::TrueParticleID> trueIDVect;
	trueIDVect.emplace_back(trueID);
	track.truth = trueIDVect;

	// Fraction of true MC hits that are captured by the reconstructed cluster
	const float completeness = (m_completenessVect != nullptr) ? (*m_completenessVect)[i] : 0.0;
	std::vector<float> truthOverlap;
	truthOverlap.emplace_back(completeness);
	track.truthOverlap = truthOverlap;
	
	// Add track to the record
	sr.nd.lar.pandora[sliceId].tracks.emplace_back(std::move(track));
	sr.nd.lar.pandora[sliceId].ntracks++;
    }
  }
    
  // ------------------------------------------------------------------------------
    void PandoraLArRecoNDBranchFiller::FillShowers(caf::StandardRecord& sr, const int nClusters) const
  {
    // Create showers for each PFO (cluster) in the event
    LOG.VERBOSE() << " Pandora LArRecoND FillShowers using " << nClusters <<" PFO clusters\n";

    for (int i = 0; i < nClusters; i++)
    {
	// Check that the PFO is a shower
	const int isShower = (*m_isShowerVect)[i];
	if (isShower == 0)
	    continue;

	// Slice id of the PFO cluster
	const int sliceId = (*m_sliceIdVect)[i];

	// Create standard record shower
	caf::SRShower shower;
	// Starting position
	const float startX = (m_startXVect != nullptr) ? (*m_startXVect)[i] : 0.0;
	const float startY = (m_startYVect != nullptr) ? (*m_startYVect)[i] : 0.0;
	const float startZ = (m_startZVect != nullptr) ? (*m_startZVect)[i] : 0.0;
	shower.start = caf::SRVector3D(startX, startY, startZ);

	// Principal axis direction
	const float dirX = (m_dirXVect != nullptr) ? (*m_dirXVect)[i] : 0.0;
	const float dirY = (m_dirYVect != nullptr) ? (*m_dirYVect)[i] : 0.0;
	const float dirZ = (m_dirZVect != nullptr) ? (*m_dirZVect)[i] : 0.0;
	shower.direction = caf::SRVector3D(dirX, dirY, dirZ);

	// Energy (GeV)
	const float energy = (m_energyVect != nullptr) ? (*m_energyVect)[i] : 0.0;
	shower.Evis = energy;

	// Use truth matching info from Pandora's LArContent hierarchy tools.
	// For LArRecoND MC SpacePoints, we offset the MCId's to make them all unique:
	// nuId = orig_nuId + 10^8, so orig_nuId = nuId - 10^8
	// mcId = orig_mcId + nuIndex*10^6, where nuIndex = 0 to N-1 neutrinos
	// nuIndex = int(mcId/10^6), so orig_mcId = mcId - nuIndex*10^6
	// Use the original MCId values
	const int mcNuId = (m_mcNuIdVect != nullptr) ? (*m_mcNuIdVect)[i] : 0;
	const int origMCNuId = (mcNuId > m_nuIdOffset) ? mcNuId - m_nuIdOffset : mcNuId;
	const int isPrimary = (m_isPrimaryVect != nullptr) ? (*m_isPrimaryVect)[i] : -1;
	const int mcId = (m_mcIdVect != nullptr) ? (*m_mcIdVect)[i] : 0;
	const int nuIndex = int(mcId/m_maxMCId);
	const int origMCId = mcId - nuIndex*m_maxMCId;

	caf::TrueParticleID trueID;
	trueID.ixn = origMCNuId;
	if (isPrimary == 1) {
	    trueID.type = caf::TrueParticleID::kPrimary;
	} else if (isPrimary == -1) {
	    trueID.type = caf::TrueParticleID::kUnknown;
	} else {
	    trueID.type = caf::TrueParticleID::kSecondary;
	}
	trueID.part = origMCId;

	// Just store the best MC match
	std::vector<caf::TrueParticleID> trueIDVect;
	trueIDVect.emplace_back(trueID);
	shower.truth = trueIDVect;

	// Fraction of true MC hits that are captured by the reconstructed cluster
	const float completeness = (m_completenessVect != nullptr) ? (*m_completenessVect)[i] : 0.0;
	std::vector<float> truthOverlap;
	truthOverlap.emplace_back(completeness);
	shower.truthOverlap = truthOverlap;

	// Add shower to the record
	sr.nd.lar.pandora[sliceId].showers.emplace_back(std::move(shower));
	sr.nd.lar.pandora[sliceId].nshowers++;
    }
  }

  // ------------------------------------------------------------------------------
  std::deque<Trigger> PandoraLArRecoNDBranchFiller::GetTriggers(int triggerType) const
  {
    if (m_Triggers.empty())
    {
	const int nEvents = m_LArRecoNDTree->GetEntries();
	LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName()
		    << "' from " << nEvents << " Pandora LArRecoND tree entries:\n";
	
	m_Triggers.reserve(nEvents);
	for (int entry = 0; entry < nEvents; entry++)
	{
	    m_LArRecoNDTree->GetEntry(entry);
	    
	    m_Triggers.emplace_back();
	    Trigger &trig = m_Triggers.back();
	    // Event number
	    trig.evtID = m_eventId;
	    // Pandora SpacePoint (SP) H5Flow-to-ROOT format doesn't store trigger type, so just select "all"
	    trig.triggerType = -1;
	    // Pandora SP format uses timestamp ticks from /charge/events/ts_start
	    trig.triggerTime_s = m_startTime;
	    // Use the same (integer) time for now
	    trig.triggerTime_ns = m_startTime;

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
