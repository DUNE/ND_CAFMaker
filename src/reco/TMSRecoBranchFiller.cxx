#include "TMSRecoBranchFiller.h"
#include "truth/FillTruth.h"

/*
 * Liam O'Sullivan <liam.osullivan@uni-mainz.de>  -  Mar 2024
 * Put together mostly from the previous example and the MINERvA RecoBranchFiller
 * This code currently assumes one track <-> one interaction; to first order
 * this is sane, as the primary use of TMS is muon collection for LAr events.
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

      TMSRecoTree->SetBranchAddress("EventNo",               &_EventNo);
      TMSRecoTree->SetBranchAddress("SliceNo",               &_SliceNo);
      TMSRecoTree->SetBranchAddress("SpillNo",               &_SpillNo);
      TMSRecoTree->SetBranchAddress("nTracks",               &_nTracks);
      TMSRecoTree->SetBranchAddress("nHits",                 _nHitsInTrack);
      TMSRecoTree->SetBranchAddress("Length",                _TrackLength);
      TMSRecoTree->SetBranchAddress("Momentum",              _TrackMomentum);
      //TMSRecoTree->SetBranchAddress("Charge",                _TrackCharge);
      TMSRecoTree->SetBranchAddress("EnergyRange",           _TrackTotalEnergy);
      TMSRecoTree->SetBranchAddress("EnergyDeposit",         _TrackEnergyDeposit);
      TMSRecoTree->SetBranchAddress("Occupancy",             _Occupancy);

      TMSRecoTree->SetBranchAddress("TrackHitPos",            _TrackRecoHitPos);
      TMSRecoTree->SetBranchAddress("StartPos",              _TrackStartPos);
      TMSRecoTree->SetBranchAddress("KalmanPos",             _TrackHitPos); // TODO: how get hits from da vectur??
      //TMSRecoTree->SetBranchAddress("HitPos",             _TrackKalmanPos); // TODO: how get hits from da vectur??
      TMSRecoTree->SetBranchAddress("EndPos",                _TrackEndPos);
      TMSRecoTree->SetBranchAddress("StartDirection",        _TrackStartDirection);
      TMSRecoTree->SetBranchAddress("EndDirection",          _TrackEndDirection);
      TMSRecoTree->SetBranchAddress("Direction",          _TrackStartDirection);
    } else {
      fTMSRecoFile = NULL;
      TMSRecoTree  = NULL;
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
    TMSRecoTree->GetEntry(i);
    LastSpillNo = _SpillNo;
    //TMSRecoTree->GetEntry(idx);//trigger.evtID);

    //caf::SRTMSInt* interaction;
    //interaction = new caf::SRTMSInt();

    sr.nd.tms.ixn.emplace_back();
    caf::SRTMSInt& interaction = sr.nd.tms.ixn.back();

    interaction.ntracks = 0;
    while (_nTracks <=1 && _SpillNo == LastSpillNo && i < TMSRecoTree->GetEntries())
    {
      TMSRecoTree->GetEntry(++i);
      interaction.tracks.resize(_nTracks + interaction.tracks.size());
      for (int j = 0; j < _nTracks; ++j) {
        interaction.ntracks++;
        interaction.tracks[j].start   = caf::SRVector3D(_TrackStartPos[j][0]/10., _TrackStartPos[j][1]/10., _TrackStartPos[j][2]/10.);;
        interaction.tracks[j].end     = caf::SRVector3D(_TrackEndPos[j][0]/10., _TrackEndPos[j][1]/10., _TrackEndPos[j][2]/10.);
        interaction.tracks[j].dir     = caf::SRVector3D(_TrackStartDirection[j][0], _TrackStartDirection[j][1] , _TrackStartDirection[j][2]);
        interaction.tracks[j].enddir  = caf::SRVector3D(_TrackEndDirection[j][0], _TrackEndDirection[j][1] , _TrackEndDirection[j][2]);

        // Calculate length by summing up the distances from the kalman reco positions
        double tmpLength_cm = 0.0;
        for (int k=0; k<_nHitsInTrack[j]-1; k++)
          tmpLength_cm += sqrt( pow(_TrackHitPos[j][k][0] - _TrackHitPos[j][k+1][0], 2)
                              + pow(_TrackHitPos[j][k][1] - _TrackHitPos[j][k+1][1], 2)
                              + pow(_TrackHitPos[j][k][2] - _TrackHitPos[j][k+1][2], 2) );

        // Track info
        interaction.tracks[j].len_cm    = tmpLength_cm; //trackVec->Mag();
        interaction.tracks[j].len_gcm2  = (_TrackLength[j]>0.0) ? _TrackLength[j]/10. : 0.0; // idk why we have negatives
        interaction.tracks[j].qual      = _Occupancy[j]; // TODO: Apparently this is a "track quality", nominally (hits in track)/(total hits)
        interaction.tracks[j].Evis      = _TrackEnergyDeposit[j];
      }
    }




    /*
    //sr.meta.tms.run = ev_run;
    //sr.meta.tms.subrun = ev_sub_run;
    //sr.meta.tms.event = ev_gate;
    sr.meta.tms.readoutstart_s = trigger.triggerTime_s;
    sr.meta.tms.readoutstart_ns = trigger.triggerTime_ns;
    //FillInteractions(truthMatcher, sr);

    //std::cout << "Trigger: " << trigger.evtID << "\t" << trigger.triggerType << "\t" << trigger.triggerTime_s << " " << trigger.triggerTime_ns << std::endl;

    sr.nd.tms.nixn = fTriggers.size(); //_nTracks;
    sr.nd.tms.ixn.resize(sr.nd.tms.nixn);
    //sr.nd.tms.ixn.reserve(3000);//TMSRecoTree->GetEntries());

    // Fill in the track info 
    //if (_nTracks > 0) std::cout << "idx: " << idx << "\ttracks: " << _nTracks << std::endl;
    std::cout << "idx: " << idx << "\ttracks: " << _nTracks << std::endl;


    std::cout << "loop ";
    for (int j = 0; j < _nTracks; ++j) {
      std::cout << j << " ";

      sr.nd.tms.ixn[i].ntracks = _nTracks; // One reco track per interaction (ixn) TODO: allow tracks to have relation in future
      //sr.nd.tms.ixn[i].tracks.resize(_nTracks);//sr.nd.tms.ixn[i].ntracks);

      // Save first and last hit in track
      // TMS Reco info is saved in mm whereas CAFs use CM as default -> do conversion here
      sr.nd.tms.ixn[i].tracks[j].start = caf::SRVector3D(_TrackStartPos[j][0]/10., _TrackStartPos[j][1]/10., _TrackStartPos[j][2]/10.);
      sr.nd.tms.ixn[i].tracks[j].end   = caf::SRVector3D(_TrackEndPos[j][0]/10., _TrackEndPos[j][1]/10., _TrackEndPos[j][2]/10.);
      // Get the directions
      // TODO: At present tracks are completely straight objects, so dir is the same here for both
      sr.nd.tms.ixn[i].tracks[j].dir     = caf::SRVector3D(_TrackStartDirection[j][0], _TrackStartDirection[j][1] , _TrackStartDirection[j][2]);
      sr.nd.tms.ixn[i].tracks[j].enddir  = caf::SRVector3D(_TrackEndDirection[j][0], _TrackEndDirection[j][1] , _TrackEndDirection[j][2]);

      //TVector3* trackVec = new TVector3( (sr.nd.tms.ixn[j].tracks[0].end - sr.nd.tms.ixn[j].tracks[0].start) );

      // Calculate length by summing up the distances from the kalman reco positions
      double tmpLength_cm = 0.0;
      for (int k=0; k<_nHitsInTrack[j]-1; k++)
        tmpLength_cm += sqrt( pow(_TrackHitPos[j][k][0] - _TrackHitPos[j][k+1][0], 2)
                            + pow(_TrackHitPos[j][k][1] - _TrackHitPos[j][k+1][1], 2)
                            + pow(_TrackHitPos[j][k][2] - _TrackHitPos[j][k+1][2], 2) );

      // Track info
      sr.nd.tms.ixn[i].tracks[j].len_cm    = tmpLength_cm; //trackVec->Mag();
      sr.nd.tms.ixn[i].tracks[j].len_gcm2  = (_TrackLength[j]>0.0) ? _TrackLength[j]/10. : 0.0; // idk why we have negatives
      sr.nd.tms.ixn[i].tracks[j].qual      = _Occupancy[j]; // TODO: Apparently this is a "track quality", nominally (hits in track)/(total hits)
      sr.nd.tms.ixn[i].tracks[j].Evis      = _TrackEnergyDeposit[j];
      sr.nd.tms.ixn[i].tracks[j].E         = _TrackTotalEnergy[j];

      //delete trackVec;
    }
         */
  }



  std::deque<Trigger> TMSRecoBranchFiller::GetTriggers(int triggerType) const
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
        //std::cout << _TrackRecoHitPos[0][0][3] << " " << (int) _TrackRecoHitPos[0][0][3] << " " << (int) (_TrackRecoHitPos[0][0][3]*1E9)%1000 << std::endl;
        trig.triggerType = 2147483647; // TODO real number

        if (entry == 0) // TODO do this less bad
          trig.triggerTime_ns = 0;
        else
          trig.triggerTime_ns = prev_trig.triggerTime_ns + 2E8 ;//+ (int) _TrackRecoHitPos[0][0][3];

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

            //std::cout << "  added trigger:  evtID=" << trig.evtID
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
