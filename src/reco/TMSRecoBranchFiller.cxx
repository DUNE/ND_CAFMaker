#include "TMSRecoBranchFiller.h"
#include "truth/FillTruth.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

namespace
{
  bool TMSRecoTimingEnabled()
  {
    static const bool enabled = []()
    {
      const char * env = std::getenv("ND_CAFMAKER_TIMING");
      if (!env) return false;
      const std::string value(env);
      return !value.empty() && value != "0" && value != "false" && value != "FALSE";
    }();
    return enabled;
  }

  class ScopedTiming
  {
    public:
      ScopedTiming(bool enabled, cafmaker::TMSRecoBranchFiller::TimingSummary & summary)
        : fEnabled(enabled), fSummary(summary)
      {
        if (fEnabled)
          fStart = std::chrono::steady_clock::now();
      }

      ~ScopedTiming()
      {
        if (!fEnabled) return;
        ++fSummary.calls;
        fSummary.totalMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - fStart).count();
      }

    private:
      bool fEnabled = false;
      cafmaker::TMSRecoBranchFiller::TimingSummary & fSummary;
      std::chrono::steady_clock::time_point fStart;
  };

  void PrintTimingSummaryLine(const char * label, const cafmaker::TMSRecoBranchFiller::TimingSummary & summary)
  {
    std::cerr << "  " << label << ": total=" << summary.totalMs << " ms"
              << ", calls=" << summary.calls;
    if (summary.calls > 0)
      std::cerr << ", avg=" << (summary.totalMs / summary.calls) << " ms/call";
    std::cerr << "\n";
  }
}

/*
 * Liam O'Sullivan <liam.osullivan@uni-mainz.de>  -  Oct 2024
 * Put together mostly from the previous example and the MINERvA RecoBranchFiller
 */

namespace cafmaker
{

  TMSRecoBranchFiller::TMSRecoBranchFiller(const std::string &tmsRecoFilename)
    : IRecoBranchFiller("TMS"),
    fTriggers(),
    fLastTriggerReqd(fTriggers.end())
  {
    fTiming.enabled = TMSRecoTimingEnabled();

    fTMSRecoFile = new TFile(tmsRecoFilename.c_str(), "READ");
    name = std::string("TMS");

    if (!fTMSRecoFile->IsZombie()) {
      SetConfigured(true);

      // Save pointer to input tree
      TMSRecoTree = dynamic_cast<TTree*>(fTMSRecoFile->Get("Reco_Tree"));
      TMSLCTree = dynamic_cast<TTree*>(fTMSRecoFile->Get("Line_Candidates"));
      if (!TMSRecoTree) {
        std::cerr << "Did not find TMS reco tree Reco_Tree in input file " << tmsRecoFilename << std::endl;
        std::cerr << "Are you sure this is a TMS reco file?" << std::endl;
        throw;
      }

      if (!TMSLCTree) {
        std::cerr << "Did not find TMS reco tree Line_Candidates in input file " << tmsRecoFilename << std::endl;
        std::cerr << "Are you sure this is a TMS reco file?" << std::endl;
        throw;
      }
      // Save pointer to truth tree
      TMSTrueTree = dynamic_cast<TTree*>(fTMSRecoFile->Get("Truth_Info"));
      if (!TMSTrueTree) {
        std::cerr << "Did not find TMS true tree Truth_Info in input file " << tmsRecoFilename << std::endl;
        std::cerr << "Are you sure this is a TMS reco file?" << std::endl;
        throw;
      }
      // Save pointer to truth spill tree
      TMSTrueSpill = dynamic_cast<TTree*>(fTMSRecoFile->Get("Truth_Spill"));
      if (!TMSTrueSpill) {
        std::cerr << "Did not find TMS true tree Truth_Spill in input file " << tmsRecoFilename << std::endl;
        std::cerr << "Are you sure this is a TMS reco file?" << std::endl;
        throw;
      }

      TMSRecoTree->SetBranchAddress("EventNo",               &_EventNo);
      TMSRecoTree->SetBranchAddress("SliceNo",               &_SliceNo);
      TMSRecoTree->SetBranchAddress("SpillNo",               &_SpillNo);
      TMSRecoTree->SetBranchAddress("nTracks",               &_nTracks);
      TMSRecoTree->SetBranchAddress("RunNo",                 &_RunNo);
      TMSRecoTree->SetBranchAddress("nHits",                 _nHitsInTrack);
      TMSRecoTree->SetBranchAddress("Length_3D",             _TrackLength);
      TMSRecoTree->SetBranchAddress("Length",                _TrackArealDensity);
      TMSRecoTree->SetBranchAddress("Momentum",              _TrackMomentum);
      TMSRecoTree->SetBranchAddress("Charge",                _TrackCharge);
      TMSRecoTree->SetBranchAddress("EnergyRange",           _TrackTotalEnergy);
      TMSRecoTree->SetBranchAddress("EnergyDeposit",         _TrackEnergyDeposit);
      //TMSRecoTree->SetBranchAddress("Time",                  _TrackTime); // TODO: Uncomment for prod >= n4p2
      //TMSRecoTree->SetBranchAddress("Occupancy",             _Occupancy); // TODO: Uncomment when Occupancy filled by TMS

      TMSRecoTree->SetBranchAddress("TrackHitPos",           _TrackRecoHitPos);
      TMSRecoTree->SetBranchAddress("StartPos",              _TrackStartPos);
      TMSRecoTree->SetBranchAddress("KalmanPos",             _TrackHitPos);
      TMSRecoTree->SetBranchAddress("EndPos",                _TrackEndPos);
      TMSRecoTree->SetBranchAddress("StartDirection",        _TrackStartDirection);
      TMSRecoTree->SetBranchAddress("EndDirection",          _TrackEndDirection);
      //TMSRecoTree->SetBranchAddress("TimeSliceStartTime",    &_TimeSliceStartTime);
      TMSLCTree->SetBranchAddress("TMSStartTime",            _TMSStartTime); // Temporary for prod n4p1, time avialable in Reco_Tree for future prods

      // Truth_Info is only used to map reco tracks to truth-particle indices.
      TMSTrueTree->SetBranchAddress("RecoTrackPrimaryParticleIndex",   _RecoTruePartId);
      TMSTrueTree->SetBranchAddress("RecoTrackSecondaryParticleIndex", _RecoTruePartIdSec);

      TMSTrueSpill->SetBranchAddress("SpillNo",                        &_TruthSpillSpillNo);
      TMSTrueSpill->SetBranchAddress("RunNo",                          &_TruthSpillRunNo);
      TMSTrueSpill->SetBranchAddress("nTrueParticles",                 &_TruthSpillNTrueParticles);
      TMSTrueSpill->SetBranchAddress("VertexID",                       _TruthSpillParticleVertexID);
      TMSTrueSpill->SetBranchAddress("Parent",                         _TruthSpillParent);
      TMSTrueSpill->SetBranchAddress("TrackId",                        _TruthSpillTrackID);
      TMSTrueSpill->SetBranchAddress("BirthPosition",                  _TruthSpillBirthPosition);
      TMSTrueSpill->SetBranchAddress("TrueVtxN",                       &_TruthSpillTrueVtxN);
      TMSTrueSpill->SetBranchAddress("TrueVtxID",                      _TruthSpillTrueVtxID);
      TMSTrueSpill->SetBranchAddress("TrueVtxX",                       _TruthSpillTrueVtxX);
      TMSTrueSpill->SetBranchAddress("TrueVtxY",                       _TruthSpillTrueVtxY);
      TMSTrueSpill->SetBranchAddress("TrueVtxZ",                       _TruthSpillTrueVtxZ);

    } else {
      fTMSRecoFile = NULL;
      TMSRecoTree  = NULL;
      TMSLCTree = NULL;
      std::cerr << "The TMS reco file you provided: " << tmsRecoFilename 
                << " appears to be a Zombie 🧟" << std::endl;
      throw;
    }
  }


  TMSRecoBranchFiller::~TMSRecoBranchFiller() {
    if (fTiming.enabled)
    {
      std::cerr << "[ND_CAFMaker timing] TMSRecoBranchFiller summary\n";
      PrintTimingSummaryLine("LoadTruthSpillEntry", fTiming.loadTruthSpillEntry);
      PrintTimingSummaryLine("BuildTruthSpillEntryMap", fTiming.buildTruthSpillEntryMap);
      PrintTimingSummaryLine("ResolveTrueInteractionIDFromVertexIndex", fTiming.resolveTrueInteractionIDFromVertexIndex);
      PrintTimingSummaryLine("ResolveRecoTrackTruthParticleIndex", fTiming.resolveRecoTrackTruthParticleIndex);
      PrintTimingSummaryLine("FindTruthSpillParticleIndex", fTiming.findTruthSpillParticleIndex);
      PrintTimingSummaryLine("ResolvePrimaryTruthParticleIndex", fTiming.resolvePrimaryTruthParticleIndex);
      PrintTimingSummaryLine("ResolveRecoTrackInteractionID", fTiming.resolveRecoTrackInteractionID);
      PrintTimingSummaryLine("GetTrueInteraction", fTiming.getTrueInteraction);
      PrintTimingSummaryLine("GetTrueParticle", fTiming.getTrueParticle);
      PrintTimingSummaryLine("FindSRTrueInteractionIndex", fTiming.findSRTrueInteractionIndex);
      PrintTimingSummaryLine("FindSRTrueParticleIndex", fTiming.findSRTrueParticleIndex);
      PrintTimingSummaryLine("_FillRecoBranches", fTiming.fillRecoBranches);
      PrintTimingSummaryLine("FillInteractions", fTiming.fillInteractions);
      PrintTimingSummaryLine("GetTriggers", fTiming.getTriggers);
      std::cerr << "  BuildTruthSpillEntryMap: entries_scanned=" << fTiming.buildTruthSpillEntryMapEntries
                << ", map_size=" << fTruthSpillEntryBySpillNo.size() << "\n";
      std::cerr << "  LoadTruthSpillEntry: misses=" << fTiming.loadTruthSpillEntryMisses << "\n";
      std::cerr << "  ResolvePrimaryTruthParticleIndex: parent_steps=" << fTiming.resolvePrimaryParentSteps
                << ", loop_fallbacks=" << fTiming.resolvePrimaryLoopFallbacks
                << ", missing_parent_fallbacks=" << fTiming.resolvePrimaryMissingParentFallbacks;
      if (fTiming.resolvePrimaryTruthParticleIndex.calls > 0)
        std::cerr << ", avg_steps/call="
                  << (static_cast<double>(fTiming.resolvePrimaryParentSteps) / fTiming.resolvePrimaryTruthParticleIndex.calls);
      std::cerr << "\n";
      std::cerr << "  FindTruthSpillParticleIndex: entries_scanned=" << fTiming.findTruthSpillParticleIndexEntriesScanned
                << ", hits=" << fTiming.findTruthSpillParticleIndexHits;
      if (fTiming.findTruthSpillParticleIndex.calls > 0)
        std::cerr << ", avg_entries/call="
                  << (static_cast<double>(fTiming.findTruthSpillParticleIndexEntriesScanned) / fTiming.findTruthSpillParticleIndex.calls);
      std::cerr << "\n";
      std::cerr << "  ResolveRecoTrackInteractionID: vertex_scans=" << fTiming.resolveRecoTrackInteractionVertexScans
                << ", exact_matches=" << fTiming.resolveRecoTrackInteractionExactMatches
                << ", nearest_fallbacks=" << fTiming.resolveRecoTrackInteractionNearestFallbacks
                << ", multi_matches=" << fTiming.resolveRecoTrackInteractionMultiMatches;
      if (fTiming.resolveRecoTrackInteractionID.calls > 0)
        std::cerr << ", avg_vertices/call="
                  << (static_cast<double>(fTiming.resolveRecoTrackInteractionVertexScans) / fTiming.resolveRecoTrackInteractionID.calls);
      std::cerr << "\n";
      std::cerr << "  _FillRecoBranches: spills_processed=" << fTiming.fillRecoBranchesSpillsProcessed
                << ", tracks_processed=" << fTiming.fillRecoBranchesTracksProcessed;
      if (fTiming.fillRecoBranches.calls > 0)
        std::cerr << ", avg_tracks/call="
                  << (static_cast<double>(fTiming.fillRecoBranchesTracksProcessed) / fTiming.fillRecoBranches.calls);
      std::cerr << "\n";
      std::cerr << "  GetTriggers: entries_scanned=" << fTiming.getTriggersEntriesScanned
                << ", triggers_created=" << fTiming.getTriggersCreated << "\n";
    }

    delete TMSRecoTree;    
    delete TMSLCTree;
    fTMSRecoFile->Close();
    delete fTMSRecoFile;
    TMSRecoTree = NULL;
    TMSLCTree = NULL;
    fTMSRecoFile = NULL;
  }

  void TMSRecoBranchFiller::LoadTruthSpillEntry(int spillNo) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.loadTruthSpillEntry);

    if (fTruthSpillEntryBySpillNo.empty())
    {
      ScopedTiming buildTiming(fTiming.enabled, fTiming.buildTruthSpillEntryMap);
      for (Long64_t entry = 0; entry < TMSTrueSpill->GetEntries(); ++entry)
      {
        ++fTiming.buildTruthSpillEntryMapEntries;
        TMSTrueSpill->GetEntry(entry);
        auto [it, inserted] = fTruthSpillEntryBySpillNo.emplace(_TruthSpillSpillNo, entry);
        if (!inserted && it->second != entry)
        {
          std::stringstream ss;
          ss << "Truth_Spill has multiple entries for SpillNo " << _TruthSpillSpillNo << "\n";
          throw std::runtime_error(ss.str());
        }
      }
    }

    auto it = fTruthSpillEntryBySpillNo.find(spillNo);
    if (it == fTruthSpillEntryBySpillNo.end())
    {
      ++fTiming.loadTruthSpillEntryMisses;
      std::stringstream ss;
      ss << "Could not find Truth_Spill entry for SpillNo " << spillNo << "\n";
      throw std::runtime_error(ss.str());
    }

    TMSTrueSpill->GetEntry(it->second);
  }

  unsigned long TMSRecoBranchFiller::ResolveTrueInteractionIDFromVertexIndex(const TruthMatcher * truthMatch, int trueVtxIdx) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.resolveTrueInteractionIDFromVertexIndex);

    if (!truthMatch)
      throw std::runtime_error("TMSRecoBranchFiller requires TruthMatcher to resolve true interaction IDs");
    if (trueVtxIdx < 0 || trueVtxIdx >= _TruthSpillTrueVtxN)
      throw std::runtime_error("Requested Truth_Spill vertex index is out of range");

    return truthMatch->ResolveVertexIDFromRunAndPosition(static_cast<unsigned long>(_TruthSpillRunNo),
                                                         _TruthSpillTrueVtxX[trueVtxIdx],
                                                         _TruthSpillTrueVtxY[trueVtxIdx],
                                                         _TruthSpillTrueVtxZ[trueVtxIdx]);
  }

  int TMSRecoBranchFiller::ResolveRecoTrackTruthParticleIndex(int recoTrackIdx) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.resolveRecoTrackTruthParticleIndex);

    if (recoTrackIdx < 0)
      throw std::runtime_error("Requested reco track index is negative");

    const int particleIdx = _RecoTruePartId[recoTrackIdx];
    if (particleIdx < 0 || particleIdx >= _TruthSpillNTrueParticles)
    {
      std::stringstream ss;
      ss << "Reco track " << recoTrackIdx << " maps to invalid Truth_Spill particle index " << particleIdx << "\n";
      throw std::runtime_error(ss.str());
    }

    return particleIdx;
  }

  int TMSRecoBranchFiller::FindTruthSpillParticleIndex(int vertexId, int trackId) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.findTruthSpillParticleIndex);

    for (int i = 0; i < _TruthSpillNTrueParticles; ++i)
    {
      ++fTiming.findTruthSpillParticleIndexEntriesScanned;
      if (_TruthSpillParticleVertexID[i] == vertexId && _TruthSpillTrackID[i] == trackId)
      {
        ++fTiming.findTruthSpillParticleIndexHits;
        return i;
      }
    }
    return -1;
  }

  int TMSRecoBranchFiller::ResolvePrimaryTruthParticleIndex(int particleIdx, int recoTrackIdx) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.resolvePrimaryTruthParticleIndex);

    std::vector<int> visited;
    while (_TruthSpillParent[particleIdx] != -1)
    {
      ++fTiming.resolvePrimaryParentSteps;
      if (std::find(visited.begin(), visited.end(), particleIdx) != visited.end())
      {
        ++fTiming.resolvePrimaryLoopFallbacks;
        LOG.WARNING() << "TMS reco track " << recoTrackIdx
                      << " encountered a loop while resolving Truth_Spill parents;"
                      << " falling back to particle index " << particleIdx << "\n";
        return particleIdx;
      }
      visited.push_back(particleIdx);

      const int parentTrackId = _TruthSpillParent[particleIdx];
      const int vertexId = _TruthSpillParticleVertexID[particleIdx];
      const int parentIdx = FindTruthSpillParticleIndex(vertexId, parentTrackId);
      if (parentIdx < 0)
      {
        ++fTiming.resolvePrimaryMissingParentFallbacks;
        LOG.WARNING() << "TMS reco track " << recoTrackIdx
                      << " could not resolve parent track ID " << parentTrackId
                      << " from Truth_Spill particle index " << particleIdx
                      << " (vertex ID " << vertexId << ");"
                      << " falling back to the current particle instead of crashing\n";
        return particleIdx;
      }
      particleIdx = parentIdx;
    }

    return particleIdx;
  }

  unsigned long TMSRecoBranchFiller::ResolveRecoTrackInteractionID(const TruthMatcher * truthMatch, int recoTrackIdx) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.resolveRecoTrackInteractionID);

    const int particleIdx = ResolvePrimaryTruthParticleIndex(ResolveRecoTrackTruthParticleIndex(recoTrackIdx), recoTrackIdx);

    const int trueVtxId = _TruthSpillParticleVertexID[particleIdx];
    const double birthX = _TruthSpillBirthPosition[particleIdx][0];
    const double birthY = _TruthSpillBirthPosition[particleIdx][1];
    const double birthZ = _TruthSpillBirthPosition[particleIdx][2];

    constexpr double kVertexMatchToleranceMm = 100.0;
    int matchedVtxIdx = -1;
    int nearestVtxIdx = -1;
    double nearestDist2 = std::numeric_limits<double>::infinity();
    for (int i = 0; i < _TruthSpillTrueVtxN; ++i)
    {
      ++fTiming.resolveRecoTrackInteractionVertexScans;
      const double dx = _TruthSpillTrueVtxX[i] - birthX;
      const double dy = _TruthSpillTrueVtxY[i] - birthY;
      const double dz = _TruthSpillTrueVtxZ[i] - birthZ;
      const double dist2 = dx*dx + dy*dy + dz*dz;
      if (dist2 < nearestDist2)
      {
        nearestDist2 = dist2;
        nearestVtxIdx = i;
      }

      if (_TruthSpillTrueVtxID[i] != trueVtxId)
        continue;
      if (dist2 > kVertexMatchToleranceMm * kVertexMatchToleranceMm)
        continue;

      if (matchedVtxIdx >= 0)
      {
        ++fTiming.resolveRecoTrackInteractionMultiMatches;
        LOG.WARNING() << "Reco track " << recoTrackIdx
                      << " matched multiple Truth_Spill vertices for vertex ID " << trueVtxId
                      << "; using the first in-tolerance match instead of crashing\n";
        continue;
      }
      matchedVtxIdx = i;
    }

    if (matchedVtxIdx >= 0)
    {
      ++fTiming.resolveRecoTrackInteractionExactMatches;
      return ResolveTrueInteractionIDFromVertexIndex(truthMatch, matchedVtxIdx);
    }

    if (nearestVtxIdx >= 0)
    {
      ++fTiming.resolveRecoTrackInteractionNearestFallbacks;
      LOG.WARNING() << "Reco track " << recoTrackIdx
                    << " could not match primary particle birth position to a Truth_Spill vertex for vertex ID "
                    << trueVtxId << " within " << kVertexMatchToleranceMm << " mm;"
                    << " falling back to nearest Truth_Spill vertex index " << nearestVtxIdx
                    << " at distance " << std::sqrt(nearestDist2) << " mm\n";
      return ResolveTrueInteractionIDFromVertexIndex(truthMatch, nearestVtxIdx);
    }

    throw std::runtime_error("Truth_Spill contains no vertices to resolve reco track interaction ID");
  }

  // here we copy all the TMS reco into the SRTMS branch of the StandardRecord object.
  void TMSRecoBranchFiller::_FillRecoBranches(const Trigger &trigger,
                                              caf::StandardRecord &sr,
                                              const cafmaker::Params &par,
                                              const TruthMatcher *truthMatcher) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.fillRecoBranches);

    sr.meta.tms.enabled = true;

    // Nicked from the MINVERvA example:
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

    int i = trigger.evtID; // pseudo-itterator for ixn

    int LastSpillNo = std::numeric_limits<int>::lowest(); // Starting value, small number so next spill number is larger
    TMSRecoTree->GetEntry(i); // Load first entry for now
    LastSpillNo = _SpillNo;

    caf::SRTMSInt *interaction;

    ++fTiming.fillRecoBranchesSpillsProcessed;
    LoadTruthSpillEntry(LastSpillNo);
    for (int i_tvtx = 0; i_tvtx < _TruthSpillTrueVtxN; ++i_tvtx)
    {
      auto neutrino_event_id = ResolveTrueInteractionIDFromVertexIndex(truthMatcher, i_tvtx);
      {
        ScopedTiming subTiming(fTiming.enabled, fTiming.getTrueInteraction);
        (void)truthMatcher->GetTrueInteraction(sr, neutrino_event_id, true); // called for side effect of registering the interaction; cast suppresses unused-return-value warning
      }
    }

    // the i index is incremented at the end of the following while()
    TMSRecoTree->GetEntry(i); // Load each subsequent entry in the spill, start from original i
    TMSTrueTree->GetEntry(i); // Keep Truth tree in sync with Reco
    TMSLCTree->GetEntry(i); 
    while (_SpillNo == LastSpillNo && i < TMSRecoTree->GetEntries()) // while we're in the spill
    {
      if (_nTracks > 0)
      {        
        for (int j = 0; j < _nTracks; ++j) {
          ++fTiming.fillRecoBranchesTracksProcessed;
          sr.nd.tms.nixn++;
          sr.nd.tms.ixn.emplace_back();

          interaction = &(sr.nd.tms.ixn.back()); // :(
          interaction->tracks.resize(1); // For now 1 track = 1 interaction; implicit assumption it's all (anti-)muons
          interaction->ntracks = 1;

          interaction->tracks[0].start   = caf::SRVector3D(_TrackStartPos[j][0]/10., _TrackStartPos[j][1]/10., _TrackStartPos[j][2]/10.);
          interaction->tracks[0].end     = caf::SRVector3D(_TrackEndPos[j][0]/10., _TrackEndPos[j][1]/10., _TrackEndPos[j][2]/10.);
          interaction->tracks[0].dir     = caf::SRVector3D(_TrackStartDirection[j][0], _TrackStartDirection[j][1] , _TrackStartDirection[j][2]);
          interaction->tracks[0].enddir  = caf::SRVector3D(_TrackEndDirection[j][0], _TrackEndDirection[j][1] , _TrackEndDirection[j][2]);

          interaction->tracks[0].time    = static_cast<double>(_TMSStartTime[j]); //Adds time of interaction // TODO: use _TrackTime after prod n4p1

          // Track info
          interaction->tracks[0].len_cm    = (_TrackLength[j]>0.0) ? _TrackLength[j]/10. : 0.0; // idk why we have negatives
          interaction->tracks[0].len_gcm2  = (_TrackArealDensity[j]>0.0) ? _TrackArealDensity[j]/10. : 0.0; // idk why we have negatives
          interaction->tracks[0].qual      = _Occupancy[j]; // TODO: Apparently this is a "track quality", nominally (hits in track)/(total hits)
          interaction->tracks[0].Evis      = _TrackEnergyDeposit[j];

          // As of Jan 2026 the Charge attribute in TMS output is the PDG value, -13 for mu+, 13 for mu-.
          // For CAF files we probably just want a +1 or -1 for the particle charge, so divide by 13 and multiply in a -
          interaction->tracks[0].charge    = -1 * _TrackCharge[j]/13;

          /*  Fill Truth
           *  The run numbers in the GHEP(?) or edep files are of the run number, followed by the event number, so we recreate that.
           * Long cos it's very long innit. Sorry. */
          const auto genieIxnID = ResolveRecoTrackInteractionID(truthMatcher, j);
          caf::SRTrueInteraction* srTrueIntPtr = nullptr;
          {
            ScopedTiming subTiming(fTiming.enabled, fTiming.getTrueInteraction);
            srTrueIntPtr = &truthMatcher->GetTrueInteraction(sr, genieIxnID);
          }
          caf::SRTrueInteraction& srTrueInt = *srTrueIntPtr;

          int srTrueIntIdx = -1;
          {
            ScopedTiming subTiming(fTiming.enabled, fTiming.findSRTrueInteractionIndex);
            srTrueIntIdx = static_cast<int>(std::distance(sr.mc.nu.begin(),
                                                          std::find_if(sr.mc.nu.begin(), sr.mc.nu.end(),
                                                                       [&srTrueInt](const caf::SRTrueInteraction& ixn)
                                                                       { return ixn.id == srTrueInt.id; })));
          }

          // TODO: Make TMS care about prim/sec tracks (check _RecoTruePartIdSec for secondaries)
          const int recoTruthParticleIdx = ResolveRecoTrackTruthParticleIndex(j);
          const int partG4ID = _TruthSpillTrackID[recoTruthParticleIdx];
          {
            ScopedTiming subTiming(fTiming.enabled, fTiming.getTrueParticle);
            truthMatcher->GetTrueParticle(sr, srTrueInt, partG4ID, true);
          }
          int truthVecIdx = -1;
          {
            ScopedTiming subTiming(fTiming.enabled, fTiming.findSRTrueParticleIndex);
            truthVecIdx = static_cast<int>(std::distance(srTrueInt.prim.begin(),
                                                         std::find_if(srTrueInt.prim.begin(), srTrueInt.prim.end(),
                                                                      [partG4ID](const caf::SRTrueParticle& part)
                                                                      { return part.G4ID == partG4ID; })));
          }

          interaction->tracks[0].truth.push_back(caf::TrueParticleID{srTrueIntIdx,
                                                                      caf::TrueParticleID::PartType::kPrimary,
                                                                      truthVecIdx});
        }
      }

      TMSRecoTree->GetEntry(++i); // Load each subsequent entry before loop test condition
      TMSTrueTree->GetEntry(  i); // Load each subsequent entry before loop test condition, i already incremented
      TMSLCTree  ->GetEntry(  i);

    }
  }

  void TMSRecoBranchFiller::FillInteractions(const TruthMatcher * truthMatch, caf::StandardRecord &sr) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.fillInteractions);

    for (int i_int = 0; i_int < _TruthSpillTrueVtxN; ++i_int)
    {
      auto neutrino_event_id = ResolveTrueInteractionIDFromVertexIndex(truthMatch, i_int);
      caf::SRTrueInteraction * srTrueIntPtr = nullptr;
      {
        ScopedTiming subTiming(fTiming.enabled, fTiming.getTrueInteraction);
        srTrueIntPtr = &truthMatch->GetTrueInteraction(sr, neutrino_event_id);
      }
      caf::SRTrueInteraction & srTrueInt = *srTrueIntPtr;
      LOG.VERBOSE() << "    --> resulting SRTrueInteraction has the following particles in it:\n";
      for (const caf::SRTrueParticle & part : srTrueInt.prim)
          LOG.VERBOSE() << "    (prim) id = " << part.G4ID << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";
      for (const caf::SRTrueParticle & part : srTrueInt.prefsi)
          LOG.VERBOSE() << "    (prefsi) id = " << part.G4ID << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";
      for (const caf::SRTrueParticle & part : srTrueInt.sec)
          LOG.VERBOSE() << "    (sec) id = " << part.G4ID  << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";
    }
  }


  // TODO: In future this nastiness will be handled by TMS
  std::deque<Trigger> TMSRecoBranchFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    ScopedTiming timing(fTiming.enabled, fTiming.getTriggers);
    std::deque<Trigger> triggers;
    int lastSpillNo = std::numeric_limits<int>::lowest(); // Starting value, small number so next spill number is larger

    if (fTriggers.empty())
    {
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "' from " << TMSRecoTree->GetEntries() << " TMS Reco_Tree:\n";
      fTriggers.reserve(TMSRecoTree->GetEntries());

      for (int entry = 0; entry < TMSRecoTree->GetEntries(); entry++)
      {
        ++fTiming.getTriggersEntriesScanned;
        TMSRecoTree->GetEntry(entry);

        if (_SpillNo == lastSpillNo)
          continue; // Only first 'event' in each spill populates a trigger

        lastSpillNo = _SpillNo;

        if (!fTriggers.empty())
        {
          Trigger & prev_trig = fTriggers.back(); // trigger before 'trig'
          fTriggers.emplace_back();               // add new trigger entry (unfilled)
          Trigger & trig      = fTriggers.back(); // trigger we're working on

          trig.evtID = entry;
          trig.triggerType = 1; // TODO real number?

          trig.triggerTime_ns = prev_trig.triggerTime_ns + 2E8 ;
          trig.triggerTime_s = prev_trig.triggerTime_s + 1; // TODO: Pull the 1.2 from correct place in file
          if (trig.triggerTime_ns >= 1E9) // If we have 1s worth of ns then add 1s and remove 1s worth of ns
          {
            trig.triggerTime_s += 1;
            trig.triggerTime_ns -= 1E9;
          }

          LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
                        << ", triggerType=" << trig.triggerType
                        << ", triggerTime_s=" << trig.triggerTime_s
                        << ", triggerTime_ns=" << trig.triggerTime_ns
                        << "\n";
        }
        else
        {
          fTriggers.emplace_back();
          Trigger & trig = fTriggers.back();
          trig.evtID = entry;
          trig.triggerType = 1; // TODO real number?
          trig.triggerTime_ns = 0;
          trig.triggerTime_s = 0;
          LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
                        << ", triggerType=" << trig.triggerType
                        << ", triggerTime_s=" << trig.triggerTime_s
                        << ", triggerTime_ns=" << trig.triggerTime_ns
                        << "\n";
        }
        ++fTiming.getTriggersCreated;

      }
      fLastTriggerReqd = fTriggers.end();  // since we just modified the list, any iterators have been invalidated
    }

    for (const Trigger & trigger : fTriggers)
    {
      if (triggerType < 0 || triggerType == fTriggers.back().triggerType)
      {
        triggers.push_back(trigger);
      }
    }

    return triggers;
  }

} // end namespace
