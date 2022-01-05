/// \file Params.h
///
///  Parameters extracted from command line and passed around


#ifndef ND_CAFMAKER_PARAMS_H
#define ND_CAFMAKER_PARAMS_H

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

namespace cafmaker
{
  struct RunConfig
  {
    fhicl::Atom<int> run    { fhicl::Name("Run"),    fhicl::Comment("Run number"), 1 };    // CAFAna doesn't like run number 0
    fhicl::Atom<int> subrun { fhicl::Name("Subrun"), fhicl::Comment("Subrun number"), 0 };

    fhicl::Atom<bool> fhc        { fhicl::Name("IsFHC"),    fhicl::Comment("Is this an FHC run?"), true };
    fhicl::Atom<bool> IsGasTPC   { fhicl::Name("IsGasTPC"), fhicl::Comment("Was GArTPC geometry used? (If not, TMS)"), false };

    fhicl::Atom<double> OA_xcoord  { fhicl::Name("OffAxisXCoord"), fhicl::Comment("Off-axis position of MPD in cm"), 0. };  // on-axis by default

  };

  struct ControlConfig
  {
    fhicl::Atom<int>  first  { fhicl::Name("FirstEvt"), fhicl::Comment("Start processing from this event number"), 0 };

    fhicl::Atom<int>  seed   { fhicl::Name("Seed"), fhicl::Comment("Random seed to use"), 7 };  // a very random number
  };

  struct PseudoRecoParams
  {

    fhicl::Atom<double> trk_muRes       { fhicl::Name("trk_muRes"),       fhicl::Comment("Fractional muon energy resolution of HP GAr TPC"),             0.02 };
    fhicl::Atom<double> LAr_muRes       { fhicl::Name("LAr_muRes"),       fhicl::Comment("Fractional muon energy resolution of muons contained in LAr"), 0.05 };
    fhicl::Atom<double> ECAL_muRes      { fhicl::Name("ECAL_muRes"),      fhicl::Comment("Fractional muon energy resolution of muons ending in ECAL"),   0.1 };
    fhicl::Atom<double> em_const        { fhicl::Name("em_const"),        fhicl::Comment("EM energy resolution constant term: A + B/sqrt(E) (GeV)"),     0.03 };
    fhicl::Atom<double> em_sqrtE        { fhicl::Name("em_sqrtE"),        fhicl::Comment("EM energy resolution 1/sqrt(E) term: A + B/sqrt(E) (GeV)"),    0.1 };
    fhicl::Atom<double> michelEff       { fhicl::Name("michelEff"),       fhicl::Comment("Michel finder efficiency"),                                    0.75 };
    fhicl::Atom<double> CC_trk_length   { fhicl::Name("CC_trk_length"),   fhicl::Comment("Minimum track length for CC (cm)"),                            100. };
    fhicl::Atom<double> pileup_frac     { fhicl::Name("pileup_frac"),     fhicl::Comment("Fraction of events with non-zero pile-up"),                    0.1 };
    fhicl::Atom<double> pileup_max      { fhicl::Name("pileup_max"),      fhicl::Comment("Maximum energy assumed for pileup events (GeV)"),              0.5 };
    fhicl::Atom<double> gastpc_len      { fhicl::Name("gastpc_len"),      fhicl::Comment("Gas TPC track length cut (cm)"),                               6. };
    fhicl::Atom<double> gastpc_B        { fhicl::Name("gastpc_B"),        fhicl::Comment("Gas TPC B field strength (Tesla)"),                            0.4 };
    fhicl::Atom<double> gastpc_padPitch { fhicl::Name("gastpc_padPitch"), fhicl::Comment("(Fixed) pad pitch of gas TPC (cm)"),                           0.1 };  // Actual pad pitch varies, which is going to be impossible to implement
    fhicl::Atom<double> gastpc_X0       { fhicl::Name("gastpc_X0"),       fhicl::Comment("Gas TPC radiation length (cm)"),                               1300. };
  };

  /// FHICL table specifying which params are accepted
  struct FhiclConfig
  {
    fhicl::Table<ControlConfig> cafmaker     { fhicl::Name("CAFMakerSettings") };
    fhicl::Table<RunConfig>     runInfo      { fhicl::Name("RunParams") };

    fhicl::Atom<bool>  grid  { fhicl::Name("OnGrid"), false };

    fhicl::Table<PseudoRecoParams> pseudoReco { fhicl::Name("PseudoRecoParams") };

  };
  using Params = fhicl::Table<cafmaker::FhiclConfig>;

  // params will be extracted from command line, and passed to the reconstruction
  struct params {
    double OA_xcoord;
    bool fhc, grid, IsGasTPC;
    int seed, run, subrun, first, n, nfiles;
    double trk_muRes, LAr_muRes, ECAL_muRes;
    double em_const, em_sqrtE;
    double michelEff;
    double CC_trk_length;
    double pileup_frac, pileup_max;
    double gastpc_len, gastpc_B, gastpc_padPitch, gastpc_X0;
  };
}

#endif //ND_CAFMAKER_PARAMS_H
