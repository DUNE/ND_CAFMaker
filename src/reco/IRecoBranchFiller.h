/// \file IRecoBranchFiller.h
///
///  Base class for reco branch fillers for the CAFMaker.
///

#ifndef ND_CAFMAKER_IRECOBRANCHFILLER_H
#define ND_CAFMAKER_IRECOBRANCHFILLER_H

#include <stdexcept>

#include "fwd.h"

namespace cafmaker
{
  class TruthMatcher;

  class IRecoBranchFiller
  {
    public:
      /// Public interface for filling.
      ///
      /// Just checks if the derived class is configured and hands off to the internal method.
      /// (Configuration might include things like setting another input file corresponding
      ///  to a specific reconstruction, for instance.)
      /// \param sr  The StandardRecord whose branches should be filled.
      void FillRecoBranches(std::size_t evtIdx,
                            caf::StandardRecord &sr,
                            const Params &par,
                            const TruthMatcher * truthMatcher=nullptr) const
      {
        if (!isConfigured)
          throw std::runtime_error("Reco branch filler hasn't been configured!");

        _FillRecoBranches(evtIdx, sr, par, truthMatcher);
      }
      std::string GetName() { return name; };

    protected:
      /// Actual implementation of reco branch filling.  Derived classes should override this.
      virtual void _FillRecoBranches(std::size_t evtIdx,
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
