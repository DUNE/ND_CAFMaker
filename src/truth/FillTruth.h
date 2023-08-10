/// \file FillTruth.h
///
/// Fill truth branches.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>, based on code by C. Marshall <chris.marshall@rochester.edu>
/// \date    Jan. 2022

#ifndef ND_CAFMAKER_FILLTRUTH_H
#define ND_CAFMAKER_FILLTRUTH_H

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>

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
  /// Convenience method for filling truth branches that does two things:
  ///  - Checks if a value contains the expected default value, and if so, copies the new value in
  ///  - If value does not contain the default, verifies that the provided new value matches the one already there
  ///
  /// \param input     The value that would be copied in if unfilled
  /// \param target    The destination value
  /// \param unsetVal  The default value expected
  template <typename InputType, typename OutputType>
  void ValidateOrCopy(const InputType & input, OutputType & target, const OutputType & unsetVal)
  {
    const auto defaultComp = [](const decltype(input) & a, const decltype(target) &b) -> bool { return static_cast<OutputType>(a) == b; };
    const auto defaultAssgn = [](const decltype(input) & a, decltype(target) &b) {  b = a; };
    return ValidateOrCopy(input, target, unsetVal, defaultComp, defaultAssgn);
  }


  /// Similar to the other variant of ValidateOrCopy(),
  /// but allowing the user to specify a function that determines if input and target are equal
  ///
  /// \param input     The value that would be copied in if unfilled
  /// \param target    The destination value
  /// \param unsetVal  The default value expected
  /// \param compFn    Function that returns true if target == input, or false otherwise
  /// \param assgnFn   Function that assigns the value of input to target
  template <typename InputType, typename OutputType>
  void ValidateOrCopy(const InputType & input, OutputType & target, const OutputType & unsetVal,
                      std::function<bool(const decltype(input) &, const decltype(target) &)> compFn,
                      std::function<void(const decltype(input) &, decltype(target) &)> assgnFn)
  {
    // vals match?  nothing more to do
    if (compFn(input, target))
     return;

    // is this the default val?
    // note that NaN isn't equal to anything, even itself, so we have to do the check differently
    bool isNan = false;
    if constexpr (std::numeric_limits<OutputType>::has_signaling_NaN)
    {
      isNan = std::isnan(target) && std::isnan(unsetVal);
    }
    if (target == unsetVal || isNan)
    {
      assgnFn(input, target);
      return;
    }

    // if neither of the above conditions were met,
    // we have a discrepancy.  bail loudly
    std::cerr << "Mismatch between branch value (" << target << ") and supplied value (" << input << ")!  Abort.\n";
    abort();
  }
  // --------------------------------------------------------------

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
      GetTrueParticle(caf::StandardRecord &sr, int ixnID, int G4ID, bool isPrimary, bool createNew = true) const;

      /// Find a TrueParticle within a given StandardRecord, or, if it doesn't exist, optionally make a new one
      ///
      /// \param sr         The caf::StandardRecord in question
      /// \param ixn        Interaction object (if you only have its ID, use the other signature of GetTrueParticle() instead)
      /// \param G4ID       TrackID of the particle from GEANT4 (or, if not propagated by GEANT4, GENIE)
      /// \param isPrimary  Was this a "primary" particle (i.e., came out of the true neutrino interaction)?
      /// \param createNew  Should a new SRTrueParticle be made if one corresponding to the given characteristics is not found?
      /// \return           The caf::SRTrueParticle that was found, or if none found and createNew is true, a new instance
      caf::SRTrueParticle &
      GetTrueParticle(caf::StandardRecord &sr, caf::SRTrueInteraction& ixn, int G4ID, bool isPrimary, bool createNew = true) const;

      /// Find a TrueInteraction within  a given StandardRecord, or, if it doesn't exist, optionally make a new one
      ///
      /// \param sr         The caf::StandardRecord in question
      /// \param ixnID      Interaction ID (should match the one coming from upstream, i.e., edep-sim)
      /// \param createNew  Should a new SRTrueInteraction be made if one corresponding to the given ID is not found?
      /// \return           The caf::SRTrueParticle that was found, or if none found and createNew is true, a new instance
      caf::SRTrueInteraction & GetTrueInteraction(caf::StandardRecord &sr, int ixnID, bool createNew = true) const;

      /// Find a TrueInteraction within  a given StandardRecord, or, if it doesn't exist, optionally make a new one.
      /// Use the 'interaction ID' variant of GetTrueInteraction() if at all possible.
      /// This version just exists to supply necessary hacks when interaction IDs from upstream are broken...
      ///
      /// \param sr           The caf::StandardRecord in question
      /// \param srTrueIxnCmp Function used to decide whether a SRTrueInteraction already in the StandardRecord matches desired criteria
      /// \param genieCmp     Function used to match a GENIE record to desired criteria
      /// \param createNew    Should a new SRTrueInteraction be made if one corresponding to the given ID is not found?
      /// \return             The caf::SRTrueParticle that was found, or if none found and createNew is true, a new instance
      caf::SRTrueInteraction & GetTrueInteraction(caf::StandardRecord &sr,
                                                  std::function<bool(const caf::SRTrueInteraction&)> srTrueIxnCmp,
                                                  std::function<bool(const genie::NtpMCEventRecord*)> genieCmp,
                                                  bool createNew = true) const;

      bool HaveGENIE() const { return fGTree != nullptr; }

    private:
      static void FillInteraction(caf::SRTrueInteraction& nu, const genie::NtpMCEventRecord * gEvt);

      mutable TTree * fGTree;
      const genie::NtpMCEventRecord * fGEvt;
  };
}
#endif //ND_CAFMAKER_FILLTRUTH_H
