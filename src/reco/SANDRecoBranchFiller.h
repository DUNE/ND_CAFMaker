///
/// Fill reco branches using SAND data.
///

#ifndef ND_CAFMAKER_SANDRECOBRANCHFILLER_H
#define ND_CAFMAKER_SANDRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

class TLorentzVector;
class TRandom;
class TVector3;

namespace cafmaker
{

  class SANDRecoBranchFiller : public IRecoBranchFiller
  {
    public:
      explicit SANDRecoBranchFiller(const std::string &sand_reco_file)
        : fSandfile(sand_reco_file)
      {
        // we don't have any further files to load or anything,
        // so we're all set
        SetConfigured(true);
      }

    protected:
      void _FillRecoBranches(std::size_t evtIdx, caf::StandardRecord &sr, const cafmaker::dumpTree &dt,
                             const cafmaker::Params &par) const override;

    private:
  //    void decayPi0( const TLorentzVector & pi0, TVector3 &gamma1, TVector3 &gamma2 ) const;

    //  void recoElectron(caf::StandardRecord & sr, const cafmaker::Params &par ) const;
    //  void recoMuonECAL(caf::StandardRecord & sr, const cafmaker::Params &par ) const;
    //  void recoMuonLAr(caf::StandardRecord & sr, const cafmaker::Params &par ) const;
    //  void recoMuonTracker(caf::StandardRecord & sr, const cafmaker::Params &par ) const;

      std::string fSandfile;
  };

}

#endif //ND_CAFMAKER_SANDRECOBRANCHFILLER_H
