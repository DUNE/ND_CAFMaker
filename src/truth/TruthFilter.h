#ifndef ND_CAFMAKER_TRUTHFILTER_H
#define ND_CAFMAKER_TRUTHFILTER_H

namespace caf
{
  class StandardRecord;
}

namespace cafmaker
{
  /// Remove from sr.mc the secondary particles below a kinetic energy threshold and/or beyond a generation,
  /// keeping the particle structure consistent:
  ///   - primaries are never removed,
  ///   - a secondary that a reco object links to is always kept,
  ///   - every kept secondary keeps its parent, and so its whole chain of ancestors up to the primary.
  /// Must be called after all the reco fillers (matchers included) have filled the StandardRecord.
  ///
  /// \param sr             The StandardRecord to filter
  /// \param minKE_GeV      Secondaries with a kinetic energy below this [GeV] are removed; <= 0 disables the cut
  /// \param maxGeneration  Secondaries of a generation beyond this are removed (the daughters of a primary are
  ///                       generation 1, their daughters 2, ...); < 0 disables the cut
  void FilterTrueParticles(caf::StandardRecord & sr, double minKE_GeV, int maxGeneration);
}

#endif //ND_CAFMAKER_TRUTHFILTER_H
