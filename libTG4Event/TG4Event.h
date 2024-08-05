//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Jul 18 21:14:56 2024 by ROOT version 6.26/07)
//      from the StreamerInfo in file /exp/dune/data/users/noeroy/mywork/MiniRun5/edepsim/MiniRun5_1E19_RHC.spill.0000000.EDEPSIM_SPILLS.root
//////////////////////////////////////////////////////////


#ifndef TG4Event_h
#define TG4Event_h
class TG4Event;

#include "Rtypes.h"
#include "TObject.h"
#include "Riostream.h"
#include <vector>
#include "TG4PrimaryVertex.h"
#include "TG4Trajectory.h"
#include <map>
namespace std {} using namespace std;
#include "TG4HitSegment.h"
#ifdef __MAKECINT__
#pragma link C++ class pair<string,vector<TG4HitSegment> >+;
#endif

class TG4Event : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   int         RunId;       //
   int         EventId;     //
   std::vector< TG4PrimaryVertex> Primaries;    // (TG4PrimaryVertex)
   std::vector< TG4Trajectory>    Trajectories;    // (TG4Trajectory)
   std::vector<std::pair<std::string,std::vector<TG4HitSegment> > > SegmentDetectors;    // (pair<string,vector<TG4HitSegment> >)

   TG4Event();
   TG4Event(TG4Event && ) = default;
   TG4Event &operator=(const TG4Event & );
   TG4Event(const TG4Event & );
   virtual ~TG4Event();

   ClassDef(TG4Event,2); // Generated by MakeProject.
};
#endif
