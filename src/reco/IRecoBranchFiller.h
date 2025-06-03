/// \file IRecoBranchFiller.h
///
///  Base class for reco branch fillers for the CAFMaker.
///

#ifndef ND_CAFMAKER_IRECOBRANCHFILLER_H
#define ND_CAFMAKER_IRECOBRANCHFILLER_H

#include <deque>
#include <stdexcept>

#include "fwd.h"
#include "util/Loggable.h"

namespace cafmaker
{
  class TruthMatcher;

  struct Trigger{


    long int          evtID;
    int               triggerType;
    unsigned long int triggerTime_s;
    unsigned int      triggerTime_ns;

    bool operator==(const Trigger & other) const  { return evtID == other.evtID && triggerType == other.triggerType; }
  };

  enum class RecoFillerType
  {
    Unknown,
    BaseReco,  ///<  Full reconstruction stack like SPINE or Pandora
    Matcher,   ///<  Post-hoc matching run across detectors, but not creating new trigger entries etc.
  };

  class IRecoBranchFiller: public Loggable
  {
    public:
      IRecoBranchFiller(std::string n, cafmaker::Logger::THRESHOLD logThresh=cafmaker::Logger::THRESHOLD::WARNING)
        : Loggable(n, logThresh), name(n)
      {}

      /// Public interface for filling.
      ///
      /// Just checks if the derived class is configured and hands off to the internal method.
      /// (Configuration might include things like setting another input file corresponding
      ///  to a specific reconstruction, for instance.)
      /// \param sr  The StandardRecord whose branches should be filled.
      void FillRecoBranches(const Trigger &trigger,
                            caf::StandardRecord &sr,
                            const Params &par,
                            const TruthMatcher * truthMatcher= nullptr) const
      {
        if (!isConfigured)
          throw std::runtime_error("Reco branch filler hasn't been configured!");

        _FillRecoBranches(trigger, sr, par, truthMatcher);
      }

      std::string GetName() const { return name; };

      /// \brief Return a list of the triggers in the input for this branch filler
      ///
      /// \param  triggerType   (Detector-specific) type of trigger to select.  <0 means "all"
      /// \return List of selected triggers (a std::deque because we're always working at the beginning or end)
      virtual std::deque<Trigger> GetTriggers(int triggerType=-1) const = 0;


      /// What type of IRecoBranchFiller is this?
      virtual RecoFillerType  FillerType() const = 0;

    protected:
      /// Actual implementation of reco branch filling.  Derived classes should override this.
      virtual void _FillRecoBranches(const Trigger &trigger,
                                     caf::StandardRecord &sr,
                                     const cafmaker::Params &par,
                                     const TruthMatcher * truthMatcher) const = 0;

      void SetConfigured(bool configured = true)    { isConfigured = configured; };

      std::string name;
    private:
      bool isConfigured;
  };
}

#endif //ND_CAFMAKER_IRECOBRANCHFILLER_H
