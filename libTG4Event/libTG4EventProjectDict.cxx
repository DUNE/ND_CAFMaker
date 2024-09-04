// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME libTG4EventProjectDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "libTG4EventProjectHeaders.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *pairlEstringcOvectorlETG4HitSegmentgRsPgR_Dictionary();
   static void pairlEstringcOvectorlETG4HitSegmentgRsPgR_TClassManip(TClass*);
   static void *new_pairlEstringcOvectorlETG4HitSegmentgRsPgR(void *p = nullptr);
   static void *newArray_pairlEstringcOvectorlETG4HitSegmentgRsPgR(Long_t size, void *p);
   static void delete_pairlEstringcOvectorlETG4HitSegmentgRsPgR(void *p);
   static void deleteArray_pairlEstringcOvectorlETG4HitSegmentgRsPgR(void *p);
   static void destruct_pairlEstringcOvectorlETG4HitSegmentgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,vector<TG4HitSegment> >*)
   {
      pair<string,vector<TG4HitSegment> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,vector<TG4HitSegment> >));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,vector<TG4HitSegment> >", "string", 208,
                  typeid(pair<string,vector<TG4HitSegment> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOvectorlETG4HitSegmentgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,vector<TG4HitSegment> >) );
      instance.SetNew(&new_pairlEstringcOvectorlETG4HitSegmentgRsPgR);
      instance.SetNewArray(&newArray_pairlEstringcOvectorlETG4HitSegmentgRsPgR);
      instance.SetDelete(&delete_pairlEstringcOvectorlETG4HitSegmentgRsPgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOvectorlETG4HitSegmentgRsPgR);
      instance.SetDestructor(&destruct_pairlEstringcOvectorlETG4HitSegmentgRsPgR);

      ::ROOT::AddClassAlternate("pair<string,vector<TG4HitSegment> >","std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::vector<TG4HitSegment, std::allocator<TG4HitSegment> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const pair<string,vector<TG4HitSegment> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOvectorlETG4HitSegmentgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const pair<string,vector<TG4HitSegment> >*)nullptr)->GetClass();
      pairlEstringcOvectorlETG4HitSegmentgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOvectorlETG4HitSegmentgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_TG4PrimaryParticle(void *p = nullptr);
   static void *newArray_TG4PrimaryParticle(Long_t size, void *p);
   static void delete_TG4PrimaryParticle(void *p);
   static void deleteArray_TG4PrimaryParticle(void *p);
   static void destruct_TG4PrimaryParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TG4PrimaryParticle*)
   {
      ::TG4PrimaryParticle *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TG4PrimaryParticle >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TG4PrimaryParticle", ::TG4PrimaryParticle::Class_Version(), "TG4PrimaryParticle.h", 18,
                  typeid(::TG4PrimaryParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TG4PrimaryParticle::Dictionary, isa_proxy, 4,
                  sizeof(::TG4PrimaryParticle) );
      instance.SetNew(&new_TG4PrimaryParticle);
      instance.SetNewArray(&newArray_TG4PrimaryParticle);
      instance.SetDelete(&delete_TG4PrimaryParticle);
      instance.SetDeleteArray(&deleteArray_TG4PrimaryParticle);
      instance.SetDestructor(&destruct_TG4PrimaryParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TG4PrimaryParticle*)
   {
      return GenerateInitInstanceLocal((::TG4PrimaryParticle*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TG4PrimaryParticle*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TG4PrimaryVertex(void *p = nullptr);
   static void *newArray_TG4PrimaryVertex(Long_t size, void *p);
   static void delete_TG4PrimaryVertex(void *p);
   static void deleteArray_TG4PrimaryVertex(void *p);
   static void destruct_TG4PrimaryVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TG4PrimaryVertex*)
   {
      ::TG4PrimaryVertex *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TG4PrimaryVertex >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TG4PrimaryVertex", ::TG4PrimaryVertex::Class_Version(), "TG4PrimaryVertex.h", 21,
                  typeid(::TG4PrimaryVertex), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TG4PrimaryVertex::Dictionary, isa_proxy, 4,
                  sizeof(::TG4PrimaryVertex) );
      instance.SetNew(&new_TG4PrimaryVertex);
      instance.SetNewArray(&newArray_TG4PrimaryVertex);
      instance.SetDelete(&delete_TG4PrimaryVertex);
      instance.SetDeleteArray(&deleteArray_TG4PrimaryVertex);
      instance.SetDestructor(&destruct_TG4PrimaryVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TG4PrimaryVertex*)
   {
      return GenerateInitInstanceLocal((::TG4PrimaryVertex*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TG4PrimaryVertex*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TG4TrajectoryPoint(void *p = nullptr);
   static void *newArray_TG4TrajectoryPoint(Long_t size, void *p);
   static void delete_TG4TrajectoryPoint(void *p);
   static void deleteArray_TG4TrajectoryPoint(void *p);
   static void destruct_TG4TrajectoryPoint(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TG4TrajectoryPoint*)
   {
      ::TG4TrajectoryPoint *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TG4TrajectoryPoint >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TG4TrajectoryPoint", ::TG4TrajectoryPoint::Class_Version(), "TG4TrajectoryPoint.h", 17,
                  typeid(::TG4TrajectoryPoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TG4TrajectoryPoint::Dictionary, isa_proxy, 4,
                  sizeof(::TG4TrajectoryPoint) );
      instance.SetNew(&new_TG4TrajectoryPoint);
      instance.SetNewArray(&newArray_TG4TrajectoryPoint);
      instance.SetDelete(&delete_TG4TrajectoryPoint);
      instance.SetDeleteArray(&deleteArray_TG4TrajectoryPoint);
      instance.SetDestructor(&destruct_TG4TrajectoryPoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TG4TrajectoryPoint*)
   {
      return GenerateInitInstanceLocal((::TG4TrajectoryPoint*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TG4TrajectoryPoint*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TG4Trajectory(void *p = nullptr);
   static void *newArray_TG4Trajectory(Long_t size, void *p);
   static void delete_TG4Trajectory(void *p);
   static void deleteArray_TG4Trajectory(void *p);
   static void destruct_TG4Trajectory(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TG4Trajectory*)
   {
      ::TG4Trajectory *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TG4Trajectory >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TG4Trajectory", ::TG4Trajectory::Class_Version(), "TG4Trajectory.h", 20,
                  typeid(::TG4Trajectory), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TG4Trajectory::Dictionary, isa_proxy, 4,
                  sizeof(::TG4Trajectory) );
      instance.SetNew(&new_TG4Trajectory);
      instance.SetNewArray(&newArray_TG4Trajectory);
      instance.SetDelete(&delete_TG4Trajectory);
      instance.SetDeleteArray(&deleteArray_TG4Trajectory);
      instance.SetDestructor(&destruct_TG4Trajectory);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TG4Trajectory*)
   {
      return GenerateInitInstanceLocal((::TG4Trajectory*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TG4Trajectory*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TG4HitSegment(void *p = nullptr);
   static void *newArray_TG4HitSegment(Long_t size, void *p);
   static void delete_TG4HitSegment(void *p);
   static void deleteArray_TG4HitSegment(void *p);
   static void destruct_TG4HitSegment(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TG4HitSegment*)
   {
      ::TG4HitSegment *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TG4HitSegment >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TG4HitSegment", ::TG4HitSegment::Class_Version(), "TG4HitSegment.h", 18,
                  typeid(::TG4HitSegment), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TG4HitSegment::Dictionary, isa_proxy, 4,
                  sizeof(::TG4HitSegment) );
      instance.SetNew(&new_TG4HitSegment);
      instance.SetNewArray(&newArray_TG4HitSegment);
      instance.SetDelete(&delete_TG4HitSegment);
      instance.SetDeleteArray(&deleteArray_TG4HitSegment);
      instance.SetDestructor(&destruct_TG4HitSegment);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TG4HitSegment*)
   {
      return GenerateInitInstanceLocal((::TG4HitSegment*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TG4HitSegment*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TG4Event(void *p = nullptr);
   static void *newArray_TG4Event(Long_t size, void *p);
   static void delete_TG4Event(void *p);
   static void deleteArray_TG4Event(void *p);
   static void destruct_TG4Event(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TG4Event*)
   {
      ::TG4Event *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TG4Event >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TG4Event", ::TG4Event::Class_Version(), "TG4Event.h", 25,
                  typeid(::TG4Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TG4Event::Dictionary, isa_proxy, 4,
                  sizeof(::TG4Event) );
      instance.SetNew(&new_TG4Event);
      instance.SetNewArray(&newArray_TG4Event);
      instance.SetDelete(&delete_TG4Event);
      instance.SetDeleteArray(&deleteArray_TG4Event);
      instance.SetDestructor(&destruct_TG4Event);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TG4Event*)
   {
      return GenerateInitInstanceLocal((::TG4Event*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TG4Event*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TG4PrimaryParticle::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TG4PrimaryParticle::Class_Name()
{
   return "TG4PrimaryParticle";
}

//______________________________________________________________________________
const char *TG4PrimaryParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4PrimaryParticle*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TG4PrimaryParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4PrimaryParticle*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TG4PrimaryParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4PrimaryParticle*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TG4PrimaryParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4PrimaryParticle*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TG4PrimaryVertex::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TG4PrimaryVertex::Class_Name()
{
   return "TG4PrimaryVertex";
}

//______________________________________________________________________________
const char *TG4PrimaryVertex::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4PrimaryVertex*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TG4PrimaryVertex::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4PrimaryVertex*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TG4PrimaryVertex::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4PrimaryVertex*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TG4PrimaryVertex::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4PrimaryVertex*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TG4TrajectoryPoint::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TG4TrajectoryPoint::Class_Name()
{
   return "TG4TrajectoryPoint";
}

//______________________________________________________________________________
const char *TG4TrajectoryPoint::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4TrajectoryPoint*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TG4TrajectoryPoint::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4TrajectoryPoint*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TG4TrajectoryPoint::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4TrajectoryPoint*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TG4TrajectoryPoint::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4TrajectoryPoint*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TG4Trajectory::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TG4Trajectory::Class_Name()
{
   return "TG4Trajectory";
}

//______________________________________________________________________________
const char *TG4Trajectory::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4Trajectory*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TG4Trajectory::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4Trajectory*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TG4Trajectory::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4Trajectory*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TG4Trajectory::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4Trajectory*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TG4HitSegment::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TG4HitSegment::Class_Name()
{
   return "TG4HitSegment";
}

//______________________________________________________________________________
const char *TG4HitSegment::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4HitSegment*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TG4HitSegment::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4HitSegment*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TG4HitSegment::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4HitSegment*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TG4HitSegment::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4HitSegment*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TG4Event::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TG4Event::Class_Name()
{
   return "TG4Event";
}

//______________________________________________________________________________
const char *TG4Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4Event*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TG4Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TG4Event*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TG4Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4Event*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TG4Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TG4Event*)nullptr)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOvectorlETG4HitSegmentgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,vector<TG4HitSegment> > : new pair<string,vector<TG4HitSegment> >;
   }
   static void *newArray_pairlEstringcOvectorlETG4HitSegmentgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,vector<TG4HitSegment> >[nElements] : new pair<string,vector<TG4HitSegment> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOvectorlETG4HitSegmentgRsPgR(void *p) {
      delete ((pair<string,vector<TG4HitSegment> >*)p);
   }
   static void deleteArray_pairlEstringcOvectorlETG4HitSegmentgRsPgR(void *p) {
      delete [] ((pair<string,vector<TG4HitSegment> >*)p);
   }
   static void destruct_pairlEstringcOvectorlETG4HitSegmentgRsPgR(void *p) {
      typedef pair<string,vector<TG4HitSegment> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class pair<string,vector<TG4HitSegment> >

//______________________________________________________________________________
void TG4PrimaryParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class TG4PrimaryParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TG4PrimaryParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(TG4PrimaryParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TG4PrimaryParticle(void *p) {
      return  p ? new(p) ::TG4PrimaryParticle : new ::TG4PrimaryParticle;
   }
   static void *newArray_TG4PrimaryParticle(Long_t nElements, void *p) {
      return p ? new(p) ::TG4PrimaryParticle[nElements] : new ::TG4PrimaryParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_TG4PrimaryParticle(void *p) {
      delete ((::TG4PrimaryParticle*)p);
   }
   static void deleteArray_TG4PrimaryParticle(void *p) {
      delete [] ((::TG4PrimaryParticle*)p);
   }
   static void destruct_TG4PrimaryParticle(void *p) {
      typedef ::TG4PrimaryParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TG4PrimaryParticle

//______________________________________________________________________________
void TG4PrimaryVertex::Streamer(TBuffer &R__b)
{
   // Stream an object of class TG4PrimaryVertex.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TG4PrimaryVertex::Class(),this);
   } else {
      R__b.WriteClassBuffer(TG4PrimaryVertex::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TG4PrimaryVertex(void *p) {
      return  p ? new(p) ::TG4PrimaryVertex : new ::TG4PrimaryVertex;
   }
   static void *newArray_TG4PrimaryVertex(Long_t nElements, void *p) {
      return p ? new(p) ::TG4PrimaryVertex[nElements] : new ::TG4PrimaryVertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_TG4PrimaryVertex(void *p) {
      delete ((::TG4PrimaryVertex*)p);
   }
   static void deleteArray_TG4PrimaryVertex(void *p) {
      delete [] ((::TG4PrimaryVertex*)p);
   }
   static void destruct_TG4PrimaryVertex(void *p) {
      typedef ::TG4PrimaryVertex current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TG4PrimaryVertex

//______________________________________________________________________________
void TG4TrajectoryPoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class TG4TrajectoryPoint.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TG4TrajectoryPoint::Class(),this);
   } else {
      R__b.WriteClassBuffer(TG4TrajectoryPoint::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TG4TrajectoryPoint(void *p) {
      return  p ? new(p) ::TG4TrajectoryPoint : new ::TG4TrajectoryPoint;
   }
   static void *newArray_TG4TrajectoryPoint(Long_t nElements, void *p) {
      return p ? new(p) ::TG4TrajectoryPoint[nElements] : new ::TG4TrajectoryPoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_TG4TrajectoryPoint(void *p) {
      delete ((::TG4TrajectoryPoint*)p);
   }
   static void deleteArray_TG4TrajectoryPoint(void *p) {
      delete [] ((::TG4TrajectoryPoint*)p);
   }
   static void destruct_TG4TrajectoryPoint(void *p) {
      typedef ::TG4TrajectoryPoint current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TG4TrajectoryPoint

//______________________________________________________________________________
void TG4Trajectory::Streamer(TBuffer &R__b)
{
   // Stream an object of class TG4Trajectory.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TG4Trajectory::Class(),this);
   } else {
      R__b.WriteClassBuffer(TG4Trajectory::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TG4Trajectory(void *p) {
      return  p ? new(p) ::TG4Trajectory : new ::TG4Trajectory;
   }
   static void *newArray_TG4Trajectory(Long_t nElements, void *p) {
      return p ? new(p) ::TG4Trajectory[nElements] : new ::TG4Trajectory[nElements];
   }
   // Wrapper around operator delete
   static void delete_TG4Trajectory(void *p) {
      delete ((::TG4Trajectory*)p);
   }
   static void deleteArray_TG4Trajectory(void *p) {
      delete [] ((::TG4Trajectory*)p);
   }
   static void destruct_TG4Trajectory(void *p) {
      typedef ::TG4Trajectory current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TG4Trajectory

//______________________________________________________________________________
void TG4HitSegment::Streamer(TBuffer &R__b)
{
   // Stream an object of class TG4HitSegment.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TG4HitSegment::Class(),this);
   } else {
      R__b.WriteClassBuffer(TG4HitSegment::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TG4HitSegment(void *p) {
      return  p ? new(p) ::TG4HitSegment : new ::TG4HitSegment;
   }
   static void *newArray_TG4HitSegment(Long_t nElements, void *p) {
      return p ? new(p) ::TG4HitSegment[nElements] : new ::TG4HitSegment[nElements];
   }
   // Wrapper around operator delete
   static void delete_TG4HitSegment(void *p) {
      delete ((::TG4HitSegment*)p);
   }
   static void deleteArray_TG4HitSegment(void *p) {
      delete [] ((::TG4HitSegment*)p);
   }
   static void destruct_TG4HitSegment(void *p) {
      typedef ::TG4HitSegment current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TG4HitSegment

//______________________________________________________________________________
void TG4Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class TG4Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TG4Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(TG4Event::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TG4Event(void *p) {
      return  p ? new(p) ::TG4Event : new ::TG4Event;
   }
   static void *newArray_TG4Event(Long_t nElements, void *p) {
      return p ? new(p) ::TG4Event[nElements] : new ::TG4Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_TG4Event(void *p) {
      delete ((::TG4Event*)p);
   }
   static void deleteArray_TG4Event(void *p) {
      delete [] ((::TG4Event*)p);
   }
   static void destruct_TG4Event(void *p) {
      typedef ::TG4Event current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TG4Event

namespace ROOT {
   static TClass *vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR_Dictionary();
   static void vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(void *p);
   static void deleteArray_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(void *p);
   static void destruct_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<string,vector<TG4HitSegment> > >*)
   {
      vector<pair<string,vector<TG4HitSegment> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<string,vector<TG4HitSegment> > >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<string,vector<TG4HitSegment> > >", -2, "vector", 386,
                  typeid(vector<pair<string,vector<TG4HitSegment> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<pair<string,vector<TG4HitSegment> > >) );
      instance.SetNew(&new_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<string,vector<TG4HitSegment> > > >()));

      ::ROOT::AddClassAlternate("vector<pair<string,vector<TG4HitSegment> > >","std::vector<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::vector<TG4HitSegment, std::allocator<TG4HitSegment> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::vector<TG4HitSegment, std::allocator<TG4HitSegment> > > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<pair<string,vector<TG4HitSegment> > >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<pair<string,vector<TG4HitSegment> > >*)nullptr)->GetClass();
      vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<string,vector<TG4HitSegment> > > : new vector<pair<string,vector<TG4HitSegment> > >;
   }
   static void *newArray_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<string,vector<TG4HitSegment> > >[nElements] : new vector<pair<string,vector<TG4HitSegment> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(void *p) {
      delete ((vector<pair<string,vector<TG4HitSegment> > >*)p);
   }
   static void deleteArray_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(void *p) {
      delete [] ((vector<pair<string,vector<TG4HitSegment> > >*)p);
   }
   static void destruct_vectorlEpairlEstringcOvectorlETG4HitSegmentgRsPgRsPgR(void *p) {
      typedef vector<pair<string,vector<TG4HitSegment> > > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<pair<string,vector<TG4HitSegment> > >

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 386,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)nullptr)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlETG4TrajectoryPointgR_Dictionary();
   static void vectorlETG4TrajectoryPointgR_TClassManip(TClass*);
   static void *new_vectorlETG4TrajectoryPointgR(void *p = nullptr);
   static void *newArray_vectorlETG4TrajectoryPointgR(Long_t size, void *p);
   static void delete_vectorlETG4TrajectoryPointgR(void *p);
   static void deleteArray_vectorlETG4TrajectoryPointgR(void *p);
   static void destruct_vectorlETG4TrajectoryPointgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TG4TrajectoryPoint>*)
   {
      vector<TG4TrajectoryPoint> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TG4TrajectoryPoint>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TG4TrajectoryPoint>", -2, "vector", 386,
                  typeid(vector<TG4TrajectoryPoint>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETG4TrajectoryPointgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TG4TrajectoryPoint>) );
      instance.SetNew(&new_vectorlETG4TrajectoryPointgR);
      instance.SetNewArray(&newArray_vectorlETG4TrajectoryPointgR);
      instance.SetDelete(&delete_vectorlETG4TrajectoryPointgR);
      instance.SetDeleteArray(&deleteArray_vectorlETG4TrajectoryPointgR);
      instance.SetDestructor(&destruct_vectorlETG4TrajectoryPointgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TG4TrajectoryPoint> >()));

      ::ROOT::AddClassAlternate("vector<TG4TrajectoryPoint>","std::vector<TG4TrajectoryPoint, std::allocator<TG4TrajectoryPoint> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TG4TrajectoryPoint>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETG4TrajectoryPointgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TG4TrajectoryPoint>*)nullptr)->GetClass();
      vectorlETG4TrajectoryPointgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETG4TrajectoryPointgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETG4TrajectoryPointgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4TrajectoryPoint> : new vector<TG4TrajectoryPoint>;
   }
   static void *newArray_vectorlETG4TrajectoryPointgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4TrajectoryPoint>[nElements] : new vector<TG4TrajectoryPoint>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETG4TrajectoryPointgR(void *p) {
      delete ((vector<TG4TrajectoryPoint>*)p);
   }
   static void deleteArray_vectorlETG4TrajectoryPointgR(void *p) {
      delete [] ((vector<TG4TrajectoryPoint>*)p);
   }
   static void destruct_vectorlETG4TrajectoryPointgR(void *p) {
      typedef vector<TG4TrajectoryPoint> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TG4TrajectoryPoint>

namespace ROOT {
   static TClass *vectorlETG4TrajectorygR_Dictionary();
   static void vectorlETG4TrajectorygR_TClassManip(TClass*);
   static void *new_vectorlETG4TrajectorygR(void *p = nullptr);
   static void *newArray_vectorlETG4TrajectorygR(Long_t size, void *p);
   static void delete_vectorlETG4TrajectorygR(void *p);
   static void deleteArray_vectorlETG4TrajectorygR(void *p);
   static void destruct_vectorlETG4TrajectorygR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TG4Trajectory>*)
   {
      vector<TG4Trajectory> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TG4Trajectory>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TG4Trajectory>", -2, "vector", 386,
                  typeid(vector<TG4Trajectory>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETG4TrajectorygR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TG4Trajectory>) );
      instance.SetNew(&new_vectorlETG4TrajectorygR);
      instance.SetNewArray(&newArray_vectorlETG4TrajectorygR);
      instance.SetDelete(&delete_vectorlETG4TrajectorygR);
      instance.SetDeleteArray(&deleteArray_vectorlETG4TrajectorygR);
      instance.SetDestructor(&destruct_vectorlETG4TrajectorygR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TG4Trajectory> >()));

      ::ROOT::AddClassAlternate("vector<TG4Trajectory>","std::vector<TG4Trajectory, std::allocator<TG4Trajectory> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TG4Trajectory>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETG4TrajectorygR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TG4Trajectory>*)nullptr)->GetClass();
      vectorlETG4TrajectorygR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETG4TrajectorygR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETG4TrajectorygR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4Trajectory> : new vector<TG4Trajectory>;
   }
   static void *newArray_vectorlETG4TrajectorygR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4Trajectory>[nElements] : new vector<TG4Trajectory>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETG4TrajectorygR(void *p) {
      delete ((vector<TG4Trajectory>*)p);
   }
   static void deleteArray_vectorlETG4TrajectorygR(void *p) {
      delete [] ((vector<TG4Trajectory>*)p);
   }
   static void destruct_vectorlETG4TrajectorygR(void *p) {
      typedef vector<TG4Trajectory> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TG4Trajectory>

namespace ROOT {
   static TClass *vectorlETG4PrimaryVertexgR_Dictionary();
   static void vectorlETG4PrimaryVertexgR_TClassManip(TClass*);
   static void *new_vectorlETG4PrimaryVertexgR(void *p = nullptr);
   static void *newArray_vectorlETG4PrimaryVertexgR(Long_t size, void *p);
   static void delete_vectorlETG4PrimaryVertexgR(void *p);
   static void deleteArray_vectorlETG4PrimaryVertexgR(void *p);
   static void destruct_vectorlETG4PrimaryVertexgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TG4PrimaryVertex>*)
   {
      vector<TG4PrimaryVertex> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TG4PrimaryVertex>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TG4PrimaryVertex>", -2, "vector", 386,
                  typeid(vector<TG4PrimaryVertex>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETG4PrimaryVertexgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TG4PrimaryVertex>) );
      instance.SetNew(&new_vectorlETG4PrimaryVertexgR);
      instance.SetNewArray(&newArray_vectorlETG4PrimaryVertexgR);
      instance.SetDelete(&delete_vectorlETG4PrimaryVertexgR);
      instance.SetDeleteArray(&deleteArray_vectorlETG4PrimaryVertexgR);
      instance.SetDestructor(&destruct_vectorlETG4PrimaryVertexgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TG4PrimaryVertex> >()));

      ::ROOT::AddClassAlternate("vector<TG4PrimaryVertex>","std::vector<TG4PrimaryVertex, std::allocator<TG4PrimaryVertex> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TG4PrimaryVertex>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETG4PrimaryVertexgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TG4PrimaryVertex>*)nullptr)->GetClass();
      vectorlETG4PrimaryVertexgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETG4PrimaryVertexgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETG4PrimaryVertexgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4PrimaryVertex> : new vector<TG4PrimaryVertex>;
   }
   static void *newArray_vectorlETG4PrimaryVertexgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4PrimaryVertex>[nElements] : new vector<TG4PrimaryVertex>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETG4PrimaryVertexgR(void *p) {
      delete ((vector<TG4PrimaryVertex>*)p);
   }
   static void deleteArray_vectorlETG4PrimaryVertexgR(void *p) {
      delete [] ((vector<TG4PrimaryVertex>*)p);
   }
   static void destruct_vectorlETG4PrimaryVertexgR(void *p) {
      typedef vector<TG4PrimaryVertex> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TG4PrimaryVertex>

namespace ROOT {
   static TClass *vectorlETG4PrimaryParticlegR_Dictionary();
   static void vectorlETG4PrimaryParticlegR_TClassManip(TClass*);
   static void *new_vectorlETG4PrimaryParticlegR(void *p = nullptr);
   static void *newArray_vectorlETG4PrimaryParticlegR(Long_t size, void *p);
   static void delete_vectorlETG4PrimaryParticlegR(void *p);
   static void deleteArray_vectorlETG4PrimaryParticlegR(void *p);
   static void destruct_vectorlETG4PrimaryParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TG4PrimaryParticle>*)
   {
      vector<TG4PrimaryParticle> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TG4PrimaryParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TG4PrimaryParticle>", -2, "vector", 386,
                  typeid(vector<TG4PrimaryParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETG4PrimaryParticlegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TG4PrimaryParticle>) );
      instance.SetNew(&new_vectorlETG4PrimaryParticlegR);
      instance.SetNewArray(&newArray_vectorlETG4PrimaryParticlegR);
      instance.SetDelete(&delete_vectorlETG4PrimaryParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlETG4PrimaryParticlegR);
      instance.SetDestructor(&destruct_vectorlETG4PrimaryParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TG4PrimaryParticle> >()));

      ::ROOT::AddClassAlternate("vector<TG4PrimaryParticle>","std::vector<TG4PrimaryParticle, std::allocator<TG4PrimaryParticle> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TG4PrimaryParticle>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETG4PrimaryParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TG4PrimaryParticle>*)nullptr)->GetClass();
      vectorlETG4PrimaryParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETG4PrimaryParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETG4PrimaryParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4PrimaryParticle> : new vector<TG4PrimaryParticle>;
   }
   static void *newArray_vectorlETG4PrimaryParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4PrimaryParticle>[nElements] : new vector<TG4PrimaryParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETG4PrimaryParticlegR(void *p) {
      delete ((vector<TG4PrimaryParticle>*)p);
   }
   static void deleteArray_vectorlETG4PrimaryParticlegR(void *p) {
      delete [] ((vector<TG4PrimaryParticle>*)p);
   }
   static void destruct_vectorlETG4PrimaryParticlegR(void *p) {
      typedef vector<TG4PrimaryParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TG4PrimaryParticle>

namespace ROOT {
   static TClass *vectorlETG4HitSegmentgR_Dictionary();
   static void vectorlETG4HitSegmentgR_TClassManip(TClass*);
   static void *new_vectorlETG4HitSegmentgR(void *p = nullptr);
   static void *newArray_vectorlETG4HitSegmentgR(Long_t size, void *p);
   static void delete_vectorlETG4HitSegmentgR(void *p);
   static void deleteArray_vectorlETG4HitSegmentgR(void *p);
   static void destruct_vectorlETG4HitSegmentgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TG4HitSegment>*)
   {
      vector<TG4HitSegment> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TG4HitSegment>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TG4HitSegment>", -2, "vector", 386,
                  typeid(vector<TG4HitSegment>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETG4HitSegmentgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TG4HitSegment>) );
      instance.SetNew(&new_vectorlETG4HitSegmentgR);
      instance.SetNewArray(&newArray_vectorlETG4HitSegmentgR);
      instance.SetDelete(&delete_vectorlETG4HitSegmentgR);
      instance.SetDeleteArray(&deleteArray_vectorlETG4HitSegmentgR);
      instance.SetDestructor(&destruct_vectorlETG4HitSegmentgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TG4HitSegment> >()));

      ::ROOT::AddClassAlternate("vector<TG4HitSegment>","std::vector<TG4HitSegment, std::allocator<TG4HitSegment> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TG4HitSegment>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETG4HitSegmentgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TG4HitSegment>*)nullptr)->GetClass();
      vectorlETG4HitSegmentgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETG4HitSegmentgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETG4HitSegmentgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4HitSegment> : new vector<TG4HitSegment>;
   }
   static void *newArray_vectorlETG4HitSegmentgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TG4HitSegment>[nElements] : new vector<TG4HitSegment>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETG4HitSegmentgR(void *p) {
      delete ((vector<TG4HitSegment>*)p);
   }
   static void deleteArray_vectorlETG4HitSegmentgR(void *p) {
      delete [] ((vector<TG4HitSegment>*)p);
   }
   static void destruct_vectorlETG4HitSegmentgR(void *p) {
      typedef vector<TG4HitSegment> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TG4HitSegment>

namespace {
  void TriggerDictionaryInitialization_libTG4EventProjectDict_Impl() {
    static const char* headers[] = {
"libTG4EventProjectHeaders.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/fhiclcpp/v4_15_03/include",
"/cvmfs/larsoft.opensciencegrid.org/products/cetlib/v3_13_04/slf7.x86_64.e20.debug/include",
"/cvmfs/larsoft.opensciencegrid.org/products/hep_concurrency/v1_07_04/slf7.x86_64.e20.debug/include",
"/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2021_1_1/Linux64bit+3.10-2.17-e20/include",
"/cvmfs/larsoft.opensciencegrid.org/products/cetlib_except/v1_07_04/include",
"/cvmfs/larsoft.opensciencegrid.org/products/boost/v1_75_0/Linux64bit+3.10-2.17-e20-debug/include",
"/cvmfs/larsoft.opensciencegrid.org/products/hdf5/v1_10_5a/Linux64bit+3.10-2.17-e20/include",
"/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include",
"/cvmfs/larsoft.opensciencegrid.org/products/xrootd/v5_4_3b/Linux64bit+3.10-2.17-e20-p3913-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/pythia/v6_4_28w/Linux64bit+3.10-2.17-gcc930-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/postgresql/v14_5/Linux64bit+3.10-2.17-p3913/include",
"/cvmfs/larsoft.opensciencegrid.org/products/openblas/v0_3_21/Linux64bit+3.10-2.17-e20/include",
"/cvmfs/larsoft.opensciencegrid.org/products/fftw/v3_3_10/Linux64bit+3.10-2.17/include",
"/cvmfs/larsoft.opensciencegrid.org/products/dk2nugenie/v01_10_01k/Linux64bit+3.10-2.17-e20-debug/include",
"/cvmfs/larsoft.opensciencegrid.org/products/dk2nudata/v01_10_01d/Linux64bit+3.10-2.17-e20-debug/include",
"/cvmfs/larsoft.opensciencegrid.org/products/genie/v3_04_00d/Linux64bit+3.10-2.17-e20-debug/include/GENIE",
"/cvmfs/larsoft.opensciencegrid.org/products/genie/v3_04_00d/Linux64bit+3.10-2.17-e20-debug/include",
"/cvmfs/larsoft.opensciencegrid.org/products/clhep/v2_4_5_1a/Linux64bit+3.10-2.17-e20-debug/include",
"/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/include",
"/cvmfs/larsoft.opensciencegrid.org/products/genie/v3_04_00d/Linux64bit+3.10-2.17-e20-debug/include/GENIE",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/etc/",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/etc//cling",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/etc//cling/plugins/include",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/include/",
"/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2021_1_1/Linux64bit+3.10-2.17-e20/include",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/include/",
"/exp/dune/app/users/noeroy/CAFS/edepsim_bugfix/reduced_edepsim_v2/ND_CAFMaker/libTG4Event/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libTG4EventProjectDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{template <class _CharT> struct __attribute__((annotate("$clingAutoload$bits/char_traits.h")))  __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace std{template <typename > class __attribute__((annotate("$clingAutoload$bits/memoryfwd.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TG4PrimaryParticle.h")))  __attribute__((annotate("$clingAutoload$libTG4EventProjectHeaders.h")))  TG4PrimaryParticle;
class __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TG4PrimaryVertex.h")))  __attribute__((annotate("$clingAutoload$libTG4EventProjectHeaders.h")))  TG4PrimaryVertex;
class __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TG4TrajectoryPoint.h")))  __attribute__((annotate("$clingAutoload$libTG4EventProjectHeaders.h")))  TG4TrajectoryPoint;
class __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TG4Trajectory.h")))  __attribute__((annotate("$clingAutoload$libTG4EventProjectHeaders.h")))  TG4Trajectory;
class __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TG4HitSegment.h")))  __attribute__((annotate("$clingAutoload$libTG4EventProjectHeaders.h")))  TG4HitSegment;
class __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Generated by MakeProject.)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TG4Event.h")))  __attribute__((annotate("$clingAutoload$libTG4EventProjectHeaders.h")))  TG4Event;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libTG4EventProjectDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "libTG4EventProjectHeaders.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TG4Event", payloadCode, "@",
"TG4HitSegment", payloadCode, "@",
"TG4PrimaryParticle", payloadCode, "@",
"TG4PrimaryVertex", payloadCode, "@",
"TG4Trajectory", payloadCode, "@",
"TG4TrajectoryPoint", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libTG4EventProjectDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libTG4EventProjectDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libTG4EventProjectDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libTG4EventProjectDict() {
  TriggerDictionaryInitialization_libTG4EventProjectDict_Impl();
}
