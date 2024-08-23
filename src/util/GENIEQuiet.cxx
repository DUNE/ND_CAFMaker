
#include "Framework/Messenger/Messenger.h"

#include "util/GENIEQuiet.h"

namespace cafmaker
{
  void QuietGENIE()
  {
    genie::Messenger::Instance()->SetPrioritiesFromXmlFile("Messenger_whisper.xml");
  }
}