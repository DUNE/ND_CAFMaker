/// \file ParameterizedBranchFiller.h
///
/// Fill reco branches using parameterized reconstruction.
///

#ifndef ND_CAFMAKER_PARAMETERIZEDRECOBRANCHFILLER_H
#define ND_CAFMAKER_PARAMETERIZEDRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

class TLorentzVector;
class TRandom;
class TVector3;

namespace cafmaker
{

  class ParameterizedRecoBranchFiller : public IRecoBranchFiller
  {
    public:
      explicit ParameterizedRecoBranchFiller(TRandom * rando)
        : fRando(rando)
      {
        // we don't have any further files to load or anything,
        // so we're all set
        SetConfigured(true);
        name = "ParameterizedReco";
      }

    protected:
      void _FillRecoBranches(std::size_t evtIdx, caf::StandardRecord &sr, const cafmaker::dumpTree &dt,
                             const cafmaker::Params &par) const override;

    private:
      void decayPi0( const TLorentzVector & pi0, TVector3 &gamma1, TVector3 &gamma2 ) const;

      void recoElectron(caf::StandardRecord & sr, const cafmaker::Params &par ) const;
      void recoMuonECAL(caf::StandardRecord & sr, const cafmaker::Params &par ) const;
      void recoMuonLAr(caf::StandardRecord & sr, const cafmaker::Params &par ) const;
      void recoMuonTracker(caf::StandardRecord & sr, const cafmaker::Params &par ) const;

      TRandom * fRando;
  };

}

#endif //ND_CAFMAKER_PARAMETERIZEDRECOBRANCHFILLER_H
