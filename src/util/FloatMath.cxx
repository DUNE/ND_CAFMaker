#include "FloatMath.h"

#include <cmath>

namespace cafmaker
{
  namespace util
  {
     bool AreEqual(double a, float b, double maxDiff, float maxRelDiff)
     {
       double A = a;
       double B = b;  // promote the float to a double.  Guaranteed by IEEE standard that this will not lose any info.

       // Check if the numbers are really close -- needed
       // when comparing numbers near zero.
       auto diff = std::abs(A - B);
       if (diff <= maxDiff)
         return true;

       A = std::abs(A);
       B = std::abs(B);
       double largest = (B > A) ? B : A;

       if (diff <= largest * maxRelDiff)
         return true;
       return false;
     }
  } // namespace util
} // namespace cafmaker
