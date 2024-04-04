/// \file Params.h
///
///  Parameters extracted from command line and passed around


#ifndef ND_CAFMAKER_PARAMS_H
#define ND_CAFMAKER_PARAMS_H

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
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
    // these are mandatory and have no default values
    fhicl::OptionalSequence<std::string> contNuGHEPFiles     { fhicl::Name{"ContainedNuGHEPFiles"},   fhicl::Comment("Input .ghep (GENIE) file(s) for in-detector neutrinos") };
//    fhicl::OptionalAtom<std::string> contNuGHEPFile     {fhicl::Name{"ContainedNuGHEPFile"},   fhicl::Comment("Input .ghep (GENIE) file for in-detector neutrinos") };
    fhicl::OptionalSequence<std::string> uncontNuGHEPFiles   {fhicl::Name{"UncontainedNuGHEPFiles"}, fhicl::Comment("Input .ghep (GENIE) file(s) for rock/hall neutrinos") };
//    fhicl::OptionalAtom<std::string> uncontNuGHEPFile   {fhicl::Name{"UncontainedNuGHEPFile"}, fhicl::Comment("Input .ghep (GENIE) file for rock/hall neutrinos") };
    fhicl::Atom<std::string> outputFile                      { fhicl::Name{"OutputFile"},           fhicl::Comment("Filename for output CAF") };

    // this one is mandatory but has a default.  (the 'fhicl.fcl' file is provided in the 'sim_inputs' directory).
    // fixme: this file is currently not used for anything, but will be once DIRT-II is done and re-enables the interaction systematics
    fhicl::Atom<std::string> nusystsFcl   { fhicl::Name{"NuSystsFCLFile"}, fhicl::Comment(".fcl configuration file for nusystematics"), "fhicl.fcl" };

    // these are optional, but will change the contents of the output CAF if supplied
    fhicl::OptionalAtom<std::string> ndlarRecoFile  { fhicl::Name{"NDLArRecoFile"}, fhicl::Comment("Input ND-LAr (ML) reco .h5 file") };
    fhicl::OptionalAtom<std::string> tmsRecoFile  { fhicl::Name{"TMSRecoFile"}, fhicl::Comment("Input TMS reco .root file") };
    fhicl::OptionalAtom<std::string> sandRecoFile  { fhicl::Name{"SANDRecoFile"}, fhicl::Comment("Input SAND reco .root file") };

    // this is optional by way of the default value. Will result in an extra output file if enabled
    fhicl::Atom<bool> makeFlatCAF { fhicl::Name{"MakeFlatCAF"}, fhicl::Comment("Make 'flat' CAF in addition to structured CAF?"), true };

    // these are optional and have defaults
    fhicl::Atom<int>  first   { fhicl::Name("FirstEvt"), fhicl::Comment("Start processing from this event number"), 0 };
    fhicl::Atom<int>  numevts { fhicl::Name("NumEvts"), fhicl::Comment("Number of events to process (-1 means 'all')"), -1 };
    fhicl::Atom<int>  seed    { fhicl::Name("Seed"), fhicl::Comment("Random seed to use"), -1 };  // use the run number by default

    // options are VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL
    fhicl::Atom<std::string> verbosity { fhicl::Name("Verbosity"), fhicl::Comment("Verbosity level of output"), "WARNING" };
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

}

#endif //ND_CAFMAKER_PARAMS_H
