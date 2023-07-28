/// \file FillTruth.h
///
/// Fill truth branches.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>, based on code by C. Marshall <chris.marshall@rochester.edu>
/// \date    Jan. 2022

#ifndef ND_CAFMAKER_FILLTRUTH_H
#define ND_CAFMAKER_FILLTRUTH_H

#include "fwd.h"

// fixme: this will need to be put back to the actual response_helper type when DIRT-II finishes model recommendations
#include <string>
namespace nusyst
{
  using response_helper = std::string;
}

// forward declarations
namespace caf
{
  class StandardRecord;
  class SRTrueParticle;
  class SRTrueInteraction;
}

namespace cafmaker
{
  class TruthMatcher
  {
    public:
      TruthMatcher(TTree * gTree, const genie::NtpMCEventRecord * gEvt);

      /// Find a TrueParticle within a given StandardRecord, or, if it doesn't exist, optionally make a new one
      ///
      /// \param sr         The caf::StandardRecord in question
      /// \param ixnID      Interaction ID (should match the one coming from upstream, i.e., edep-sim)
      /// \param G4ID       TrackID of the particle from GEANT4 (or, if not propagated by GEANT4, GENIE)
      /// \param isPrimary  Was this a "primary" particle (i.e., came out of the true neutrino interaction)?
      /// \param createNew  Should a new SRTrueParticle be made if one corresponding to the given characteristics is not found?
      /// \return           The caf::SRTrueParticle that was found, or if none found and createNew is true, a new instance
      caf::SRTrueParticle &
      GetTrueParticle(caf::StandardRecord &sr, int ixnID, int G4ID, bool isPrimary, bool createNew = true);

      /// Find a TrueInteraction within  a given StandardRecord, or, if it doesn't exist, optionally make a new one
      ///
      /// \param sr         The caf::StandardRecord in question
      /// \param ixnID      Interaction ID (should match the one coming from upstream, i.e., edep-sim)
      /// \param createNew  Should a new SRTrueInteraction be made if one corresponding to the given ID is not found?
      /// \return           The caf::SRTrueParticle that was found, or if none found and createNew is true, a new instance
      caf::SRTrueInteraction &GetTrueInteraction(caf::StandardRecord &sr, int ixnID, bool createNew = true);

    private:
      static void FillInteraction(caf::SRTrueInteraction& nu, const genie::NtpMCEventRecord * gEvt);

      TTree * fGTree;
      const genie::NtpMCEventRecord * fGEvt;
  };
}
#endif //ND_CAFMAKER_FILLTRUTH_H
