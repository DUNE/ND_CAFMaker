/// \file IRecoBranchFiller.h
///
///  Base class for reco branch fillers for the CAFMaker.
///

#include <stdexcept>

#ifndef ND_CAFMAKER_IRECOBRANCHFILLER_H
#define ND_CAFMAKER_IRECOBRANCHFILLER_H

namespace caf
{
  class StandardRecord;
}

namespace cafmaker
{
  class dumpTree;
  class params;

  class IRecoBranchFiller
  {
    public:
      /// Public interface for filling.
      ///
      /// Just checks if the derived class is configured and hands off to the internal method.
      /// (Configuration might include things like setting another input file corresponding
      ///  to a specific reconstruction, for instance.)
      /// \param sr  The StandardRecord whose branches should be filled.
      void FillRecoBranches(caf::StandardRecord& sr, const cafmaker::dumpTree & dt, const cafmaker::params &par) const
      {
        if (!isConfigured)
          throw std::runtime_error("Reco branch filler hasn't been configured!");

        _FillRecoBranches(sr, dt, par);
      }

    protected:
      /// Actual implementation of reco branch filling.  Derived classes should override this.
      virtual void _FillRecoBranches(caf::StandardRecord& sr, const cafmaker::dumpTree & dt, const cafmaker::params &par) const = 0;

      void SetConfigured(bool configured = true)    { isConfigured = configured; };

    private:
      bool isConfigured;
  };
}

#endif //ND_CAFMAKER_IRECOBRANCHFILLER_H
