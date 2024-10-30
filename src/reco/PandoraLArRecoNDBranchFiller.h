/// Fill Pandora LArRecoND branches
///
/// \author  John Back <J.J.Back@warwick.ac.uk>
/// \date    Aug 2024

#ifndef ND_CAFMAKER_PandoraLArRecoNDBranchFiller_H
#define ND_CAFMAKER_PandoraLArRecoNDBranchFiller_H

#include <iostream>
#include <vector>
#include <algorithm>

// The virtual base class
#include "reco/IRecoBranchFiller.h"
#include "truth/FillTruth.h"

// ROOT headers
#include "TFile.h"
#include "TTree.h"

// duneanaobj
#include "duneanaobj/StandardRecord/StandardRecord.h"

namespace cafmaker
{
  class PandoraLArRecoNDBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      PandoraLArRecoNDBranchFiller(const std::string & pandoraLArRecoNDFilename);

      std::deque<Trigger> GetTriggers(int triggerType) const  override;

      RecoFillerType FillerType() const override { return RecoFillerType::BaseReco; }

      ~PandoraLArRecoNDBranchFiller();

    private:
      void _FillRecoBranches(const Trigger &trigger,
			     caf::StandardRecord &sr,
			     const cafmaker::Params &par,
			     const TruthMatcher *truthMatch= nullptr) const override;

      void FillTracks(caf::StandardRecord &sr, const int nClusters, const TruthMatcher *truthMatch) const;
      void FillShowers(caf::StandardRecord &sr, const int nClusters, const TruthMatcher *truthMatch) const;
      
      TFile *m_LArRecoNDFile;
      TTree *m_LArRecoNDTree;

      int m_eventId;
      int m_run;
      int m_subRun;
      int m_unixTime;
      int m_startTime;
      std::vector<int> *m_isShowerVect = nullptr;
      std::vector<int> *m_sliceIdVect = nullptr;
      std::vector<float> *m_startXVect = nullptr;
      std::vector<float> *m_startYVect = nullptr;
      std::vector<float> *m_startZVect = nullptr;
      std::vector<float> *m_endXVect = nullptr;
      std::vector<float> *m_endYVect = nullptr;
      std::vector<float> *m_endZVect = nullptr;
      std::vector<float> *m_dirXVect = nullptr;
      std::vector<float> *m_dirYVect = nullptr;
      std::vector<float> *m_dirZVect = nullptr;
      std::vector<float> *m_energyVect = nullptr;
      std::vector<int> *m_n3DHitsVect = nullptr;
      std::vector<long> *m_mcNuIdVect = nullptr;
      std::vector<long> *m_mcLocalIdVect = nullptr;
      std::vector<int> *m_isPrimaryVect = nullptr;
      std::vector<float> *m_completenessVect = nullptr;

      float m_LArRho{1.3973f}; // LAr density (g/cm3)

      int m_nuIdOffset{100000000};
      int m_maxMCId{1000000};

      mutable std::vector<cafmaker::Trigger> m_Triggers;
      mutable decltype(m_Triggers)::const_iterator  m_LastTriggerReqd; ///< the last trigger requested using _FillRecoBranches
  };

}
#endif //ND_CAFMAKER_PandoraLArRecoNDBranchFiller_H
