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

      // Add Truth tree for the index of the true primary particles
      //TMSTrueTree->SetBranchAddress("VertexID",                       _TrueVtxId);
      //TMSTrueTree->SetBranchAddress("TrueVtxID",                      _TrueVtxId);
      //TMSTrueTree->SetBranchAddress("TrueVtxN",                       &_TrueVtxN);
      TMSTrueTree->SetBranchAddress("TrueVtxX",                       _TrueVtxX);
      TMSTrueTree->SetBranchAddress("TrueVtxY",                       _TrueVtxY);
      TMSTrueTree->SetBranchAddress("TrueVtxZ",                       _TrueVtxZ);
      TMSTrueTree->SetBranchAddress("RecoTrackPrimaryParticleVtxId",  _RecoTrueVtxId);
      TMSTrueTree->SetBranchAddress("RecoTrackPrimaryParticleIndex",  _RecoTruePartId);
      TMSTrueTree->SetBranchAddress("RecoTrackSecondaryParticleIndex", _RecoTruePartIdSec);

      TMSTrueSpill->SetBranchAddress("VertexID",                       _TrueVtxId);
      TMSTrueSpill->SetBranchAddress("TrueVtxN",                       &_TrueVtxN);
      TMSTrueSpill->SetBranchAddress("RunNo",                          &_TrueRunNo);
      //TMSTrueSpill->SetBranchAddress("TrueVtxX",                       _TrueVtxX);
      //TMSTrueSpill->SetBranchAddress("TrueVtxY",                       _TrueVtxY);
      //TMSTrueSpill->SetBranchAddress("TrueVtxZ",                       _TrueVtxZ);
      //TMSTrueSpill->SetBranchAddress("RecoTrackPrimaryParticleVtxId",  _RecoTrueVtxId);
      //TMSTrueSpill->SetBranchAddress("RecoTrackPrimaryParticleIndex",  _RecoTruePartId);
      //TMSTrueSpill->SetBranchAddress("RecoTrackSecondaryParticleIndex", _RecoTruePartIdSec);

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

  // ---------------------------------------------------------------------------

  void TMSRecoBranchFiller::FillTrueInteraction(caf::SRTrueInteraction & srTrueInt,
                                                long int trkid) const
  {
    LOG.DEBUG() << "    now copying truth info from TMS TrueInteraction to SRTrueInteraction...\n";

    const auto NaN = std::numeric_limits<float>::signaling_NaN();

    // here we are converting from mm (units from MINERvA) to cm
    //ValidateOrCopy((int) trkid, (int) srTrueInt.id, NaN, "SRTrueInteraction::id");
    // TODO: uncomment these?
    //ValidateOrCopy(_TrueVtxX[trkid]/10., srTrueInt.vtx.x, NaN, "SRTrueInteraction::vtx::x");
    //ValidateOrCopy(_TrueVtxY[trkid]/10., srTrueInt.vtx.y, NaN, "SRTrueInteraction::vtx::y");
    //ValidateOrCopy(_TrueVtxZ[trkid]/10., srTrueInt.vtx.z, NaN, "SRTrueInteraction::vtx::z");
    // Try to manually fill instead of this copy crap
    //srTrueInt.vtx.x = _TrueVtxX[trkid]/10.;
    //srTrueInt.vtx.y = _TrueVtxY[trkid]/10.;
    //srTrueInt.vtx.z = _TrueVtxZ[trkid]/10.;


  }

  void TMSRecoBranchFiller::FillTrueParticle(caf::SRTrueParticle & srTruePart,
                                             long int trkid) const
  {
    const auto NaN = std::numeric_limits<float>::signaling_NaN();
    //ValidateOrCopy(mc_traj_pdg[max_trkid], srTruePart.pdg, 0, "pdg_code");

    //ValidateOrCopy(mc_traj_edepsim_trkid[trkid], srTruePart.G4ID, -1, "SRTrueParticle::track_id");
    //ValidateOrCopy(mc_traj_parentid[trkid], srTruePart.parent, -1, "SRTrueParticle::parent");
    //ValidateOrCopy(mc_traj_point_px[trkid][0]/1000., srTruePart.p.px, NaN, "SRTrueParticle::p.px");
    //ValidateOrCopy(mc_traj_point_py[trkid][0]/1000., srTruePart.p.py, NaN, "SRTrueParticle::p.py");
    //ValidateOrCopy(mc_traj_point_pz[trkid][0]/1000., srTruePart.p.pz, NaN, "SRTrueParticle::p.pz");
    //srTruePart.p.px =
    //srTruePart.p.py =
    //srTruePart.p.pz =

    try
    {
      //ValidateOrCopy(mc_traj_point_E[trkid][0] / 1000., srTruePart.p.E, NaN, "SRTrueParticle::p.E");
      0;
    }
    catch (std::runtime_error & e)
    {
      auto diff = 10 ; //(mc_traj_point_E[trkid][0] / 1000. - srTruePart.p.E);
      if (diff < 1) // < 1 MeV
      {
        LOG.WARNING() << "True particle energy (track id=" << srTruePart.G4ID << ", pdg=" << srTruePart.pdg << ", stored E=" << srTruePart.p.E << ")"
                      << " differs by " << diff << " MeV between stored (GENIE?) and ML-reco pass-through values";
      }
      else
        throw e;
    }



    // todo: Things do not match yet the exact Genie output, need to work on the Minerva reconstruction output.
    // For now will assume that if it's a primary and we gor the eventID and trackid right, Genie will have fill it properly
    // And if it's an important secondary shared by both detectors, MLReco will have filled it.
    /*
    ValidateOrCopy(mc_traj_point_x[trkid][0]/10. - offsetX/10., srTruePart.start_pos.x, NaN, "SRTrueParticle::start_pos.x");
    ValidateOrCopy(mc_traj_point_y[trkid][0]/10. - offsetY/10., srTruePart.start_pos.y, NaN, "SRTrueParticle::start_pos.y");
    ValidateOrCopy(mc_traj_point_z[trkid][0]/10. - offsetZ/10., srTruePart.start_pos.z, NaN, "SRTrueParticle::start_pos.z");

    ValidateOrCopy(mc_traj_point_x[trkid][1]/10. - offsetX/10., srTruePart.end_pos.x, NaN, "SRTrueParticle::end_pos.x");
    ValidateOrCopy(mc_traj_point_y[trkid][1]/10. - offsetY/10., srTruePart.end_pos.y, NaN, "SRTrueParticle::end_pos.y");
    ValidateOrCopy(mc_traj_point_z[trkid][1]/10. - offsetZ/10., srTruePart.end_pos.z, NaN, "SRTrueParticle::end_pos.z");


    ValidateOrCopy(mc_traj_point_px[trkid][0], srTruePart.p.px, NaN, "SRTrueParticle::end_pos.x");
    ValidateOrCopy(mc_traj_point_py[trkid][0], srTruePart.p.py, NaN, "SRTrueParticle::end_pos.y");
    ValidateOrCopy(mc_traj_point_pz[trkid][0], srTruePart.p.pz, NaN, "SRTrueParticle::end_pos.z");

    */
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
    // Get nth entry from tree

    int LastSpillNo = -999999; //_SpillNo;
    TMSRecoTree->GetEntry(i); // Load first entry for now
    LastSpillNo = _SpillNo;

    caf::SRTMSInt *interaction;// = sr.nd.tms.ixn.back();
    caf::TrueParticleID *truePartID;
    caf::SRTrueParticle *srTruePart;
    caf::SRTrueInteraction  srTrueInt;

    const auto NaN = std::numeric_limits<float>::signaling_NaN();

    // Fill Truth parts first?
    for (int i_tru=0; i_tru< TMSTrueSpill->GetEntries(); i_tru++)
    {
      TMSTrueSpill->GetEntry(i_tru); // Keep Truth spill tree in sync with Reco

      for (int i_tvtx=0; i_tvtx<_TrueVtxN; i_tvtx++)
      {
        unsigned long int neutrino_event_id = (unsigned long) ((_TrueRunNo)*1E6 + _TrueVtxId[i_tvtx]);//mc_int_edepsimId[i_int];
        srTrueInt = truthMatcher->GetTrueInteraction(sr, neutrino_event_id, true);
        FillTrueInteraction(srTrueInt, i_tvtx); // Does nothing atm anyway
      }
    }

    unsigned total = 0; // Total number of tracks in the interaction
    TMSRecoTree->GetEntry(i); // Load each subsequent entry in the spill, start from original i
    TMSTrueTree->GetEntry(i); // Keep Truth tree in sync with Reco

    // TODO: Add true info entering TMS
    while (_SpillNo == LastSpillNo && i < TMSRecoTree->GetEntries()) // while we're in the spill
    {
      if (_nTracks > 0)
      {
        for (int j = 0; j < _nTracks; ++j) {

          sr.nd.tms.ixn.emplace_back();
          interaction = &(sr.nd.tms.ixn.back());
          interaction->tracks.resize(1);//_nTracks);
          interaction->ntracks = 1;

          truePartID = new caf::TrueParticleID();
          srTruePart = new caf::SRTrueParticle();

          //interaction.ntracks++;
          interaction->tracks[0].start   = caf::SRVector3D(_TrackStartPos[j][0]/10., _TrackStartPos[j][1]/10., _TrackStartPos[j][2]/10.);
          interaction->tracks[0].end     = caf::SRVector3D(_TrackEndPos[j][0]/10., _TrackEndPos[j][1]/10., _TrackEndPos[j][2]/10.);
          interaction->tracks[0].dir     = caf::SRVector3D(_TrackStartDirection[j][0], _TrackStartDirection[j][1] , _TrackStartDirection[j][2]);
          interaction->tracks[0].enddir  = caf::SRVector3D(_TrackEndDirection[j][0], _TrackEndDirection[j][1] , _TrackEndDirection[j][2]);

          // Track info
          interaction->tracks[0].len_gcm2  = (_TrackLength[j]>0.0) ? _TrackLength[j]/10. : 0.0; // idk why we have negatives
          interaction->tracks[0].qual      = _Occupancy[j]; // TODO: Apparently this is a "track quality", nominally (hits in track)/(total hits)
          interaction->tracks[0].Evis      = _TrackEnergyDeposit[j];

          // Fill Truth
          // TODO: (unsigned long) (_RunNo*1E6 + _RecoTruePartId[j]) ... what am I smoking.
          // The run numbers in the GHEP(?) or edep files are of the run number, followed by the event number, so we recreate that. Long cos it's very long innit. Sorry.
          //srTrueInt = &(truthMatcher->GetTrueInteraction(sr, (unsigned long) (_RunNo*1E6 + _RecoTruePartId[j]), false)); // Pointer to the object
          srTrueInt = (truthMatcher->GetTrueInteraction(sr, (unsigned long) ((_RunNo%100000)*1E6 + _RecoTrueVtxId[j])));//, false)); // Pointer to the object
	        //FillTrueParticle(*srTruePart, (unsigned long) (_RunNo*1E6 + _RecoTruePartId[j])); // Does nothing anyway?
          truePartID->type = caf::TrueParticleID::kPrimary;
          //truePartID.type = is_primary ? caf::TrueParticleID::kPrimary : caf::TrueParticleID::kSecondary; // TODO: Make TMS care about prim/sec tracks
          truePartID->ixn  = (long int) (_RunNo*1E6 + _RecoTrueVtxId[j]);
          srTruePart->interaction_id = (long int) (_RunNo*1E6 + _RecoTrueVtxId[j]);

          interaction->tracks[0].truth.push_back(std::move(*truePartID)); // TODO Unfuck
        }
      }

      TMSRecoTree->GetEntry(++i); // Load each subsequent entry before loop test condition
      TMSTrueTree->GetEntry(  i); // Load each subsequent entry before loop test condition
    }
  }

  void TMSRecoBranchFiller::FillInteractions(const TruthMatcher * truthMatch, caf::StandardRecord &sr) const
  {

    for (int i_int = 0; i_int<_TrueVtxN; i_int++)
    {
      //Long_t neutrino_event_id = _TrueVtxId[i_int];//mc_int_edepsimId[i_int];
      unsigned long int neutrino_event_id = (unsigned long) (_RunNo*1E6 + _TrueVtxId[i_int]);//mc_int_edepsimId[i_int];
      caf::SRTrueInteraction & srTrueInt = truthMatch->GetTrueInteraction(sr, neutrino_event_id);//, false);
      LOG.VERBOSE() << "    --> resulting SRTrueInteraction has the following particles in it:\n";
      for (const caf::SRTrueParticle & part : srTrueInt.prim)
          LOG.VERBOSE() << "    (prim) id = " << part.G4ID << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";
      for (const caf::SRTrueParticle & part : srTrueInt.prefsi)
          LOG.VERBOSE() << "    (prefsi) id = " << part.G4ID << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";
      for (const caf::SRTrueParticle & part : srTrueInt.sec)
          LOG.VERBOSE() << "    (sec) id = " << part.G4ID  << " pdg = " << part.pdg << ", energy = " << part.p.E << "\n";

      // here we need to fill in any additional info
      // that GENIE didn't know about: e.g., secondary particles made by GEANT4
      FillTrueInteraction(srTrueInt, i_int);
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
