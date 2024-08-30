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
#include <map>
#include <memory>
#include <sstream>

#include "fwd.h"
#include "util/Loggable.h"
#include "util/FloatMath.h"
  //TG4Event
#include "TG4Event.h"

  // fixme: this will need to be put back to the actual response_helper type when DIRT-II finishes model recommendations
#include <string>

  namespace nusyst
  {
    using response_helper = std::string;
  }

  // forward declarations
  class TTree;
  class TFile;

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
    void ValidateOrCopy(const InputType & input, OutputType & target, const OutputType & unsetVal, const std::string & fieldName="")
    {
      const auto defaultComp = [](std::add_const_t<std::remove_reference_t<decltype(input)>> & a,
                                 std::add_const_t<std::remove_reference_t<decltype(target)>> &b) -> bool { return static_cast<OutputType>(a) == b; };
      const auto defaultAssgn = [](std::add_const_t<std::remove_reference_t<decltype(input)>> & a, decltype(target) &b) {  b = a; };
      ValidateOrCopy(input, target, unsetVal, defaultComp, defaultAssgn, fieldName);
    }

   // --------------------------------------------------------------

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
                        std::function<bool(std::add_const_t<std::remove_reference_t<decltype(input)>> &,
                                           std::add_const_t<std::remove_reference_t<decltype(target)>> &)> compFn,
                        std::function<void(std::add_const_t<std::remove_reference_t<decltype(input)>> &,
                                           decltype(target) &)> assgnFn,
                        const std::string & fieldName="")
    {
      LOG_S("ValidateOrCopy()").VERBOSE() << "     " << (!fieldName.empty() ? "field='" + fieldName + "';" : "")
                                          << " supplied val=" << input << "; previous branch val=" << target << "; default=" << unsetVal << "\n";

      // is the target value already the desired value?
      // or was the supplied value the default value (which implies nothing should be set)?
      // then nothing more to do.
      if (compFn(input, target) || compFn(input, unsetVal)) return;

      // note that NaN and inf aren't equal to anything, even themselves, so we have check that differently
      if constexpr (std::numeric_limits<InputType>::has_signaling_NaN && std::numeric_limits<OutputType>::has_signaling_NaN)
        if ( (std::isnan(input) && std::isnan(target)) || (std::isnan(input) && std::isnan(unsetVal)) ) return;
      if constexpr ( std::numeric_limits<InputType>::has_infinity && std::numeric_limits<OutputType>::has_infinity )
        if ( (std::isinf(input) && std::isinf(target)) || (std::isinf(input) && std::isinf(unsetVal)) ) return;

      // is this the default val?
      bool isNanInf = false;
      if constexpr (std::numeric_limits<OutputType>::has_signaling_NaN)
        isNanInf = (std::isnan(target) && std::isnan(unsetVal));
      if constexpr ( std::numeric_limits<OutputType>::has_infinity )
        isNanInf = isNanInf || (std::isinf(target) && std::isinf(unsetVal));
      if (target == unsetVal || isNanInf)
      {
        assgnFn(input, target);
        return;
      }

      // if neither of the above conditions were met,
      // we have a discrepancy.  bail loudly
      std::stringstream ss;
      ss << (!fieldName.empty() ? "For field name '" + fieldName + "': " : "")
         << "Mismatch between branch value (" << target << ") and supplied value (" << input << ")!  Abort.\n";
      throw std::runtime_error(ss.str());
    }

    // --------------------------------------------------------------
    // specialize the template for double-> float conversions, which we do a lot,
    // and which have roundoff problems in the comparison operator otherwise
    template <>
  void ValidateOrCopy<double, float>(const double & input, float & target, const float & unsetVal, const std::string& fieldName);
 
    template <>
  void ValidateOrCopy<float, float>(const float & input, float & target, const float & unsetVal, const std::string& fieldName); 
  //same for long int and int (true interaction id)
  template <>
  void ValidateOrCopy<int, long int>(const int & input, long int & target, const long int & unsetVal, const std::string & fieldName);

  // --------------------------------------------------------------

  class TruthMatcher : public cafmaker::Loggable
  {
    public:
      TruthMatcher(const std::vector<std::string> & ghepFilenames,
                  std::string edepsimFilename,
                   const genie::NtpMCEventRecord *gEvt,
                   std::function<int(const genie::NtpMCEventRecord *)> genieFillerCallback);

      /// Find a TrueParticle within a given StandardRecord, or, if it doesn't exist, optionally make a new one
      ///
      /// \param sr         The caf::StandardRecord in question
      /// \param ixnID      Interaction ID (should match the one coming from upstream, i.e., edep-sim)
      /// \param G4ID       TrackID of the particle from GEANT4 (or, if not propagated by GEANT4, GENIE)
      /// \param isPrimary  Was this a "primary" particle (i.e., came out of the true neutrino interaction)?
      /// \param createNew  Should a new SRTrueParticle be made if one corresponding to the given characteristics is not found?
      /// \return           The caf::SRTrueParticle that was found, or if none found and createNew is true, a new instance
      caf::SRTrueParticle &
      GetTrueParticle(caf::StandardRecord& sr, int ixnID, int G4ID, bool isPrimary, bool createNew = true) const;

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

      /// Find a TrueParticle within a given StandardRecord, or, if it doesn't exist, optionally make a new one
      ///
      /// \param sr         The caf::StandardRecord in question
      /// \param ixn        Interaction object (if you only have its ID, use the other signature of GetTrueParticle() instead)
      /// \param G4ID       TrackID of the particle from GEANT4 (or, if not propagated by GEANT4, GENIE)
      /// \param cmp        Function used to decide whether a SRTrueParticle already in the SRTrueInteraction matches desired criteria
      /// \param isPrimary  Was this a "primary" particle (i.e., came out of the true neutrino interaction)?
      /// \param createNew  Should a new SRTrueParticle be made if one corresponding to the given characteristics is not found?
      /// \return           The caf::SRTrueParticle that was found, or if none found and createNew is true, a new instance
      caf::SRTrueParticle &
      GetTrueParticle(caf::StandardRecord &sr,
                      caf::SRTrueInteraction& ixn,
                      int G4ID,
                      std::function<bool(const caf::SRTrueParticle&)> cmp,
                      bool isPrimary,
                      bool createNew = true) const;

      /// Find a TrueInteraction within  a given StandardRecord, or, if it doesn't exist, optionally make a new one
      ///
      /// \param sr         The caf::StandardRecord in question
      /// \param ixnID      Interaction ID (should match the one coming from upstream, i.e., edep-sim)
      /// \param createNew  Should a new SRTrueInteraction be made if one corresponding to the given ID is not found?
      /// \return           The caf::SRTrueParticle that was found, or if none found and createNew is true, a new instance
      caf::SRTrueInteraction & GetTrueInteraction(caf::StandardRecord & sr, unsigned long ixnID, bool createNew = true) const;
      bool HaveGENIE() const;

      void SetLogThrehsold(cafmaker::Logger::THRESHOLD thresh) override;

    private:
    static void FillInteraction(caf::SRTrueInteraction& nu, const genie::NtpMCEventRecord * gEvt, const TG4Event * g4event, int nixn);
    // static void FillParticle(caf::SRTrueParticle * part, std::size_t nixn, const TG4Event * g4event);
    static  int FillParticle(caf::SRTrueInteraction &ixn, std::size_t nixn, int G4ID, std::vector<caf::SRTrueParticle> & collection, int & counter, const TG4Event * g4event);



      /// Internal class organizing the GENIE trees by run to make them more easily accessible
      class GTreeContainer : public cafmaker::Loggable
      {
        public:
          GTreeContainer(const std::vector<std::string> & filenames, const genie::NtpMCEventRecord * gEvt=nullptr);

          // container interface
          const auto begin()  { return fGTrees.begin(); }
          const auto end()    { return fGTrees.end(); }

          /// Select the GENIE event in the known trees corresponding to a particular run and entry number.
          /// If no such event is found, throws an exception.
          void SelectEvent(unsigned long int runNum, unsigned int evtNum);

          const genie::NtpMCEventRecord * GEvt() const;
          void SetGEvtAddr(const genie::NtpMCEventRecord * evt);

        private:
          const genie::NtpMCEventRecord * fGEvt;
          std::map<unsigned long int, TTree*> fGTrees;
          std::vector<std::unique_ptr<TFile>> fGFiles;
      };

      mutable GTreeContainer fGTrees;
      std::function<int(const genie::NtpMCEventRecord *)> fGENIEWriterCallback;  ///< Callback function that'll write a copy of a GENIE event out to storage

      // Geant4 file (To read secondaries)

      class EdepSimTreeContainer : public cafmaker::Loggable
      {
        public:

          EdepSimTreeContainer(std::string filename);
          void SelectEvent(unsigned long int runNum, unsigned int evtNum);
          void SelectEvent(unsigned long int vertex_id);
          const TG4Event * G4Event() const;
          const TTree * GetEdepTree() const;
          void  LoadTree();

        private:
          TFile * fEdepFile;
          TTree * fEdepTree;
          std::map<unsigned long int, int> fEdepEntries;
          const TG4Event * fG4Event;
          bool f_isTreeLoaded;
      };
      mutable EdepSimTreeContainer fEdepSimTree;

  };
}
#endif //ND_CAFMAKER_FILLTRUTH_H
