#include <pricer/stdafx.h>
#include "UComparator.h"

namespace Pricer
{
  namespace UComparator
  {
    bool AlmostEqual(double x, double y)
    {
      return std::abs(x - y) < COMPAR_ERROR;
    }
  }
}