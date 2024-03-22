#ifndef ND_CAFMAKER_FLOATMATH_H
#define ND_CAFMAKER_FLOATMATH_H

#include <limits>

namespace cafmaker
{
  namespace util
  {
    /// Adapted from https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    /// See also https://github.com/randomascii/blogstuff/blob/main/FloatingPoint/CompareAsInt/CompareAsInt.cpp
    bool AreEqual(double a, float b, double maxDiff = 1e-6, float maxRelDiff = 1e-5);
  }
}

#endif //ND_CAFMAKER_FLOATMATH_H
