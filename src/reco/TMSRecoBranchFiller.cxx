#include "TMSRecoBranchFiller.h"
#include "truth/FillTruth.h"

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

      TMSRecoTree->SetBranchAddress("EventNo",               &_EventNo);
      TMSRecoTree->SetBranchAddress("SliceNo",               &_SliceNo);
      TMSRecoTree->SetBranchAddress("SpillNo",               &_SpillNo);
      TMSRecoTree->SetBranchAddress("RunNo",                 &_RunNo);
      TMSRecoTree->SetBranchAddress("nTracks",               &_nTracks);
      TMSRecoTree->SetBranchAddress("nHits",                 _nHitsInTrack);
      TMSRecoTree->SetBranchAddress("Length",                _TrackLength);
      TMSRecoTree->SetBranchAddress("Momentum",              _TrackMomentum);
      //TMSRecoTree->SetBranchAddress("Charge",                _TrackCharge); // TODO: Uncomment when Occupancy filled by TMS
      TMSRecoTree->SetBranchAddress("EnergyRange",           _TrackTotalEnergy);
      TMSRecoTree->SetBranchAddress("EnergyDeposit",         _TrackEnergyDeposit);
      //TMSRecoTree->SetBranchAddress("Occupancy",             _Occupancy); // TODO: Uncomment when Occupancy filled by TMS

      TMSRecoTree->SetBranchAddress("TrackHitPos",            _TrackRecoHitPos);
      TMSRecoTree->SetBranchAddress("StartPos",              _TrackStartPos);
      TMSRecoTree->SetBranchAddress("KalmanPos",             _TrackHitPos);
      TMSRecoTree->SetBranchAddress("EndPos",                _TrackEndPos);
      TMSRecoTree->SetBranchAddress("StartDirection",        _TrackStartDirection);
      TMSRecoTree->SetBranchAddress("EndDirection",          _TrackEndDirection);
      TMSLCTree->SetBranchAddress("TMSStartTime",            &_TMSStartTime);
    // Add Truth tree for the index of the true primary particles
      TMSTrueTree->SetBranchAddress("RecoTrackPrimaryParticleVtxId", _RecoTrueVtxId);
      TMSTrueTree->SetBranchAddress("RecoTrackPrimaryParticleIndex", _RecoTruePartId);
      TMSTrueTree->SetBranchAddress("RecoTrackSecondaryParticleIndex", _RecoTruePartIdSec);

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
    delete TMSRecoTree;    
    delete TMSLCTree;
    fTMSRecoFile->Close();
    delete fTMSRecoFile;
    TMSRecoTree = NULL;
    TMSLCTree = NULL;
    fTMSRecoFile = NULL;
  }

  // ---------------------------------------------------------------------------

  // here we copy all the TMS reco into the SRTMS branch of the StandardRecord object.
  void TMSRecoBranchFiller::_FillRecoBranches(const Trigger &trigger,
                                              caf::StandardRecord &sr,
                                              const cafmaker::Params &par,
                                              const TruthMatcher *truthMatcher) const
  {

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
    // Get nth entry from tree

    int LastSpillNo = -999999; //_SpillNo;
    TMSRecoTree->GetEntry(i); // Load first entry for now
    TMSLCTree->GetEntry(i);
    LastSpillNo = _SpillNo;

    sr.nd.tms.ixn.emplace_back();
    caf::SRTMSInt& interaction = sr.nd.tms.ixn.back();

    sr.nd.tms.nixn += 1; //Make sure to update nixn
    caf::TrueParticleID truePartID;
    caf::SRTrueParticle *srTruePart;
    caf::SRTrueInteraction *srTrueInt;

    unsigned total = 0; // Total number of tracks in the interaction
    interaction.ntracks = 0;
    TMSRecoTree->GetEntry(i); // Load each subsequent entry in the spill, start from original i
    TMSTrueTree->GetEntry(i); // Keep Truth tree in sync with Reco
    TMSRecoTree->GetEntry(i); 
    while (_SpillNo == LastSpillNo && i < TMSRecoTree->GetEntries()) // while we're in the spill
    {
      if (_nTracks > 0)
      {        

        total = interaction.tracks.size();
        interaction.tracks.resize(_nTracks + interaction.tracks.size());
        for (int j = 0; j < _nTracks; ++j) {
          interaction.ntracks++;
          interaction.tracks[total+j].start   = caf::SRVector3D(_TrackStartPos[j][0]/10., _TrackStartPos[j][1]/10., _TrackStartPos[j][2]/10.);;
          interaction.tracks[total+j].end     = caf::SRVector3D(_TrackEndPos[j][0]/10., _TrackEndPos[j][1]/10., _TrackEndPos[j][2]/10.);
          interaction.tracks[total+j].dir     = caf::SRVector3D(_TrackStartDirection[j][0], _TrackStartDirection[j][1] , _TrackStartDirection[j][2]);
          interaction.tracks[total+j].enddir  = caf::SRVector3D(_TrackEndDirection[j][0], _TrackEndDirection[j][1] , _TrackEndDirection[j][2]);

          interaction.tracks[total+j].time    = _TMSStartTime[j]; //Adds time of interaction

          // Track info
          //interaction.tracks[total+j].len_cm    = tmpLength_cm; //trackVec->Mag(); // TODO: Coming Soon™
          interaction.tracks[total+j].len_gcm2  = (_TrackLength[j]>0.0) ? _TrackLength[j]/10. : 0.0; // idk why we have negatives
          interaction.tracks[total+j].qual      = _Occupancy[j]; // TODO: Apparently this is a "track quality", nominally (hits in track)/(total hits)
          interaction.tracks[total+j].Evis      = _TrackEnergyDeposit[j];

          // Fill Truth
          // TODO: (unsigned long) (_RunNo*1E6 + _RecoTruePartId[j]) ... what am I smoking.
          // The run numbers in the GHEP(?) or edep files are of the run number, followed by the event number, so we recreate that. Long cos it's very long innit. Sorry.

          srTrueInt = &(truthMatcher->GetTrueInteraction(sr, (unsigned long) (_RunNo*1E6 + _RecoTruePartId[j]), true)); // Pointer to the object
          truePartID.ixn  = (long int) (_RunNo*1E6 + _RecoTrueVtxId[j]);
          //truePartID.type = is_primary ? caf::TrueParticleID::kPrimary : caf::TrueParticleID::kSecondary; // TODO: Make TMS care about prim/sec tracks
          truePartID.type = caf::TrueParticleID::kPrimary;

          interaction.tracks[total+j].truth.push_back(std::move(truePartID));
        }
      }

      TMSRecoTree->GetEntry(++i); // Load each subsequent entry before loop test condition
      TMSTrueTree->GetEntry(  i); // Load each subsequent entry before loop test condition
      TMSLCTree->GetEntry(  i);

    }
  }

  // TODO: In future this nastiness will be handled by TMS
  std::deque<Trigger> TMSRecoBranchFiller::GetTriggers(int triggerType, bool beamOnly) const
  {
    std::deque<Trigger> triggers;
    int lastSpillNo = -99999999;

    if (fTriggers.empty())
    {
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "' from " << TMSRecoTree->GetEntries() << " TMS Reco_Tree:\n";
      fTriggers.reserve(TMSRecoTree->GetEntries());

      for (int entry = 0; entry < TMSRecoTree->GetEntries(); entry++)
      {
        TMSRecoTree->GetEntry(entry);

        if (_SpillNo == lastSpillNo)
          continue; // Only first 'event' in each spill populates a trigger

        lastSpillNo = _SpillNo;

        Trigger & prev_trig = fTriggers.back(); // trigger before 'trig'
        fTriggers.emplace_back();               // add new trigger entry (unfilled)
        Trigger & trig      = fTriggers.back(); // trigger we're working on

        trig.evtID = entry;
        trig.triggerType = 1; // TODO real number?

        if (entry == 0) // TODO do this less bad
          trig.triggerTime_ns = 0;
        else
          trig.triggerTime_ns = prev_trig.triggerTime_ns + 2E8 ;

        if (entry == 0) // TODO do this less bad
          trig.triggerTime_s = 0;
        else
        {
          trig.triggerTime_s = prev_trig.triggerTime_s + 1; // TODO: Pull the 1.2 from correct place in file
          if (trig.triggerTime_ns >= 1E9)
          {
            trig.triggerTime_s += 1;
            trig.triggerTime_ns -= 1E9;
          }
        }

        LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
                      << ", triggerType=" << trig.triggerType
                      << ", triggerTime_s=" << trig.triggerTime_s
                      << ", triggerTime_ns=" << trig.triggerTime_ns
                      << "\n";
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
