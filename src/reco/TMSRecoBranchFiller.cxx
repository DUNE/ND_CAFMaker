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
      if (!TMSRecoTree) {
        std::cerr << "Did not find TMS reco tree Reco_Tree in input file " << tmsRecoFilename << std::endl;
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
      TMSRecoTree->SetBranchAddress("Time",                  _TrackTime);
      //TMSRecoTree->SetBranchAddress("Occupancy",             _Occupancy); // TODO: Uncomment when Occupancy filled by TMS

      TMSRecoTree->SetBranchAddress("TrackHitPos",           _TrackRecoHitPos);
      TMSRecoTree->SetBranchAddress("StartPos",              _TrackStartPos);
      TMSRecoTree->SetBranchAddress("KalmanPos",             _TrackHitPos);
      TMSRecoTree->SetBranchAddress("EndPos",                _TrackEndPos);
      TMSRecoTree->SetBranchAddress("StartDirection",        _TrackStartDirection);
      TMSRecoTree->SetBranchAddress("EndDirection",          _TrackEndDirection);

      // Add Truth tree for the index of the true primary particles
      TMSTrueTree->SetBranchAddress("TrueVtxX",                       _TrueVtxX);
      TMSTrueTree->SetBranchAddress("TrueVtxY",                       _TrueVtxY);
      TMSTrueTree->SetBranchAddress("TrueVtxZ",                       _TrueVtxZ);
      TMSTrueTree->SetBranchAddress("RecoTrackPrimaryParticleVtxId",  _RecoTrueVtxId);
      TMSTrueTree->SetBranchAddress("RecoTrackPrimaryParticleIndex",  _RecoTruePartId);
      TMSTrueTree->SetBranchAddress("RecoTrackSecondaryParticleIndex", _RecoTruePartIdSec);

      TMSTrueSpill->SetBranchAddress("VertexID",                       _TrueVtxId);
      TMSTrueSpill->SetBranchAddress("TrueVtxN",                       &_TrueVtxN);
      TMSTrueSpill->SetBranchAddress("RunNo",                          &_TrueRunNo);

    } else {
      fTMSRecoFile = NULL;
      std::cerr << "The TMS reco file you provided: " << tmsRecoFilename 
                << " appears to be a Zombie 🧟" << std::endl;
      throw;
    }
  }


  TMSRecoBranchFiller::~TMSRecoBranchFiller() {
    delete TMSRecoTree;
    fTMSRecoFile->Close();
    delete fTMSRecoFile;
    TMSRecoTree = NULL;
    fTMSRecoFile = NULL;
  }

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

    int LastSpillNo = -999999; // Starting value, very negative number so next spill number is larger
    TMSRecoTree->GetEntry(i); // Load first entry for now
    LastSpillNo = _SpillNo;

    caf::SRTMSInt *interaction;
    caf::TrueParticleID *truePartID;
    caf::SRTrueParticle *srTruePart;
    caf::SRTrueInteraction  srTrueInt;

    // Fill Truth parts first?
    for (int i_tru=0; i_tru< TMSTrueSpill->GetEntries(); i_tru++)
    {
      TMSTrueSpill->GetEntry(i_tru);

      for (int i_tvtx=0; i_tvtx<_TrueVtxN; i_tvtx++)
      {
        unsigned long int neutrino_event_id = (unsigned long) ((_TrueRunNo)*1E6 + _TrueVtxId[i_tvtx]);
        srTrueInt = truthMatcher->GetTrueInteraction(sr, neutrino_event_id, true);
      }
    }

    // the i index is incremented at the end of the following while()
    TMSRecoTree->GetEntry(i); // Load each subsequent entry in the spill, start from original i
    TMSTrueTree->GetEntry(i); // Keep Truth tree in sync with Reco

    // TODO: Add true info at the face of TMS?
    while (_SpillNo == LastSpillNo && i < TMSRecoTree->GetEntries()) // while we're in the spill
    {
      if (_nTracks > 0) // and we have reco tracks
      {
        for (int j = 0; j < _nTracks; ++j) {
          sr.nd.tms.nixn++;
          sr.nd.tms.ixn.emplace_back();

          interaction = &(sr.nd.tms.ixn.back()); // :(
          interaction->tracks.resize(1); // For now 1 track = 1 interaction; implicit assumption it's all (anti-)muons
          interaction->ntracks = 1;

          truePartID = new caf::TrueParticleID();
          srTruePart = new caf::SRTrueParticle();

          interaction->tracks[0].start   = caf::SRVector3D(_TrackStartPos[j][0]/10., _TrackStartPos[j][1]/10., _TrackStartPos[j][2]/10.);
          interaction->tracks[0].end     = caf::SRVector3D(_TrackEndPos[j][0]/10., _TrackEndPos[j][1]/10., _TrackEndPos[j][2]/10.);
          interaction->tracks[0].dir     = caf::SRVector3D(_TrackStartDirection[j][0], _TrackStartDirection[j][1] , _TrackStartDirection[j][2]);
          interaction->tracks[0].enddir  = caf::SRVector3D(_TrackEndDirection[j][0], _TrackEndDirection[j][1] , _TrackEndDirection[j][2]);

          // Track info
          interaction->tracks[0].len_cm    = (_TrackLength[j]>0.0) ? _TrackLength[j]/10. : 0.0; // idk why we have negatives
          interaction->tracks[0].len_gcm2  = (_TrackArealDensity[j]>0.0) ? _TrackArealDensity[j]/10. : 0.0; // idk why we have negatives
          interaction->tracks[0].qual      = _Occupancy[j]; // TODO: Apparently this is a "track quality", nominally (hits in track)/(total hits)
          interaction->tracks[0].Evis      = _TrackEnergyDeposit[j];

          // As of Jan 2026 the Charge attribute in TMS output is the PDG value, -13 for mu+, 13 for mu-.
          // For CAF files we probably just want a +1 or -1 for the particle charge, so divide by 13 and multiply in a -
          //interaction->tracks[0].charge    = -1 * _TrackCharge[j]/13; // TODO: UNCOMMENT BEFORE MERGE, REQUIRES NEW DUNEANAOBJ BUILD

          /*  Fill Truth
           *  The run numbers in the GHEP(?) or edep files are of the run number, followed by the event number, so we recreate that.
           * Long cos it's very long innit. Sorry. */
          srTrueInt = (truthMatcher->GetTrueInteraction(sr, (unsigned long) ((_RunNo%100000)*1E6 + _RecoTrueVtxId[j])));//, false)); // Pointer to the object
          truePartID->type = caf::TrueParticleID::kPrimary;
          //truePartID.type = is_primary ? caf::TrueParticleID::kPrimary : caf::TrueParticleID::kSecondary; // TODO: Make TMS care about prim/sec tracks
          truePartID->ixn  = (long int) (_RunNo*1E6 + _RecoTrueVtxId[j]);
          srTruePart->interaction_id = (long int) (_RunNo*1E6 + _RecoTrueVtxId[j]);

          interaction->tracks[0].truth.push_back(std::move(*truePartID)); // TODO Unfuck
        }
      }

      TMSRecoTree->GetEntry(++i); // Load each subsequent entry before loop test condition
      TMSTrueTree->GetEntry(  i); // Load each subsequent entry before loop test condition, i already incremented
    }
  }

  void TMSRecoBranchFiller::FillInteractions(const TruthMatcher * truthMatch, caf::StandardRecord &sr) const
  {

    for (int i_int = 0; i_int<_TrueVtxN; i_int++)
    {
      //Long_t neutrino_event_id = _TrueVtxId[i_int];//mc_int_edepsimId[i_int];
      unsigned long int neutrino_event_id = (unsigned long) (_RunNo*1E6 + _TrueVtxId[i_int]);
      caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, neutrino_event_id);
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

        if (entry == 0)
          trig.triggerTime_ns = 0;
        else
          trig.triggerTime_ns = prev_trig.triggerTime_ns + 2E8 ;

        if (entry == 0)
          trig.triggerTime_s = 0;
        else
        {
          trig.triggerTime_s = prev_trig.triggerTime_s + 1; // TODO: Pull the 1.2 from correct place in file
          if (trig.triggerTime_ns >= 1E9) // If we have 1s worth of ns then add 1s and remove 1s worth of ns
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
