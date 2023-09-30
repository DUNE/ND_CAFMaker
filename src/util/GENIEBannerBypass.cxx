/// \file GENIEBannerBypass.cxx
///
/// Work around the obnoxious GENIE splash screen
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    Sep. 2023

#include <string>

typedef unsigned int UInt_t;

using namespace std;

namespace genie
{
  namespace utils
  {
    namespace print
    {
      /// Null implementation of the function that spits out the GENIE splash screen.
      /// Because the CAFMaker library is loaded before the GENIE ones,
      /// this implementation is discovered first and used instead of the official GENIE one.
      /// This way, we don't get the GENIE printout cluttering our output.
      void PrintBanner(string, UInt_t)
      {}
    }
  }
}
