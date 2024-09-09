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
        std::cerr << "Did not find TMS reco tree Line_Candidates in input file " << tmsRecoFilename << std::endl;
        std::cerr << "Are you sure this is a TMS reco file?" << std::endl;
        throw;
      }

      TMSRecoTree->SetBranchAddress("EventNo",               &_EventNo);
      TMSRecoTree->SetBranchAddress("SliceNo",               &_SliceNo);
      TMSRecoTree->SetBranchAddress("SpillNo",               &_SpillNo);
      TMSRecoTree->SetBranchAddress("nTracks",               &_nTracks);
      TMSRecoTree->SetBranchAddress("nHits",                 _nHitsInTrack);
      TMSRecoTree->SetBranchAddress("Length",                _TrackLength);
      TMSRecoTree->SetBranchAddress("Charge",                _TrackCharge);
      TMSRecoTree->SetBranchAddress("EnergyRange",           _TrackTotalEnergy);
      TMSRecoTree->SetBranchAddress("EnergyDeposit",         _TrackEnergyDeposit);
      TMSRecoTree->SetBranchAddress("Occupancy",             _Occupancy);

      //TMSRecoTree->SetBranchAddress("HitPos",           _TrackHitPos); // TODO: how get hits from da vectur??
      //TMSRecoTree->SetBranchAddress("RecoHitPos",       _TrackRecoHitPos);
      TMSRecoTree->SetBranchAddress("StartPos",              _TrackStartPos);
      TMSRecoTree->SetBranchAddress("EndPos",                _TrackEndPos);
      TMSRecoTree->SetBranchAddress("Direction",             _TrackDirection);
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
#ifndef DISABLE_TMS

    sr.meta.tms.enabled = true;

    // Nicked from the MINVERvA example:
    // figure out where in our list of triggers this event index is.
    // we should always be looking forwards, since we expect to be traversing in that direction
//    auto it_start = (fLastTriggerReqd == fTriggers.end()) ? fTriggers.cbegin() : fLastTriggerReqd;
//    auto itTrig = std::find(it_start, fTriggers.cend(), trigger);
//    if (itTrig == fTriggers.end())
//    {
//      LOG.FATAL() << "Reco branch filler '" << GetName() << "' could not find trigger with evtID == " << trigger.evtID << "!  Abort.\n";
//      abort();
//    }
//    std::size_t idx = std::distance(fTriggers.cbegin(), itTrig);
//    LOG.VERBOSE() << "    Reco branch filler '" << GetName() << "', trigger.evtID == " << trigger.evtID << ", internal evt idx = " << idx << ".\n";


    // Get nth entry from tree
    TMSRecoTree->GetEntry(trigger.evtID);

    // First set number of tracks
    sr.nd.tms.nixn = _nTracks;
    sr.nd.tms.ixn.resize(sr.nd.tms.nixn);

    //sr.nd.tms.nixn = _nLines;
    //sr.SRTMSInt.ntracks= _nLines;

    // Fill in the track info 
    for (int i = 0; i < _nTracks; ++i) {
      /* Currently we only really care about catching (anti-)muons from LAr,
       * so we assume that each track is probably a(n) (anti-)muon originating
       * from a single interaction upstream in LAr.
       *     Liam     */

      sr.nd.tms.ixn[i].ntracks = 1; // One reco track per interaction (ixn) TODO: allow tracks to have relation in future
      sr.nd.tms.ixn[i].tracks.resize(sr.nd.tms.ixn[i].ntracks);

      // Save first and last hit in track
      // TMS Reco info is saved in mm whereas CAFs use CM as default -> do conversion here
      sr.nd.tms.ixn[i].tracks[0].start = caf::SRVector3D(_TrackStartPos[i][0]/10., _TrackStartPos[i][1]/10., _TrackStartPos[i][2]/10.);
      sr.nd.tms.ixn[i].tracks[0].end   = caf::SRVector3D(_TrackEndPos[i][0]/10., _TrackEndPos[i][1]/10., _TrackEndPos[i][2]/10.);
      // Get the directions
      // TODO: At present tracks are completely straight objects, so dir is the same here for both
      sr.nd.tms.ixn[i].tracks[0].dir     = caf::SRVector3D(_TrackDirection[i][0], _TrackDirection[i][1] , _TrackDirection[i][2]);
      sr.nd.tms.ixn[i].tracks[0].enddir  = caf::SRVector3D(_TrackDirection[i][0], _TrackDirection[i][1] , _TrackDirection[i][2]);

      TVector3* trackVec = new TVector3( (sr.nd.tms.ixn[i].tracks[0].end - sr.nd.tms.ixn[i].tracks[0].start) );

      // Track info
      sr.nd.tms.ixn[i].tracks[0].len_cm    = trackVec->Mag();
      sr.nd.tms.ixn[i].tracks[0].len_gcm2  = _TrackLength[i]/10.;
      sr.nd.tms.ixn[i].tracks[0].qual      = _Occupancy[i]; // TODO: Apparently this is a "track quality", nominally (hits in track)/(total hits)
      sr.nd.tms.ixn[i].tracks[0].Evis      = _TrackEnergyDeposit[i];
      sr.nd.tms.ixn[i].tracks[0].E         = _TrackTotalEnergy[i];

      delete trackVec;
    }
  }



  std::deque<Trigger> TMSRecoBranchFiller::GetTriggers(int triggerType) const
  {
    std::deque<Trigger> triggers;
    if (fTriggers.empty())
    {
      LOG.DEBUG() << "Loading triggers with type " << triggerType << " within branch filler '" << GetName() << "' from " << TMSRecoTree->GetEntries() << " TMS Reco_Tree:\n";
      fTriggers.reserve(TMSRecoTree->GetEntries());

      for (int entry = 0; entry < TMSRecoTree->GetEntries(); entry++)
      {
        // TODO: BIIIIIIG TODO

        TMSRecoTree->GetEntry(entry);

        fTriggers.emplace_back();
        Trigger & trig = fTriggers.back();

        trig.evtID = Long_t(_EventNo);


        // todo: these are placeholder values until we can propagate enough info through the reco files
        LOG.VERBOSE() << "  added trigger:  evtID=" << trig.evtID
                      << "\n";

      }
      fLastTriggerReqd = fTriggers.end();  // since we just modified the list, any iterators have been invalidated
    }

    for (const Trigger & trigger : fTriggers)
    {
      if (triggerType < 0 || triggerType == fTriggers.back().triggerType)
        triggers.push_back(trigger);
    }

    return triggers;
  }

#endif  // DISABLE_TMS

} // end namespace
