#include <pricer/stdafx.h>
#include "UFctInter.h"

namespace Pricer
{
  namespace UFctInter
  {
    fctInter2Points InterInVar()
    {
      return[](double x, double x1, double x2, double y1, double y2) -> double
      {
        // TODO currently kind of handle the case when x is negative
        // we have to do better now 
        if (x2 < 0.0)
          throw std::exception("X2 have to be positive");

        // TODO Change equal by almost (x1 == x) by almost equal
        if ((x <= 0.0) || (x1 == x))
          return y1;

        // TODO do we check x2 == x1
        // can we make different level of security by defining 
        // macro #ifdef ...
        return std::sqrt((x - x1) / (x2 - x1)
          * (y2 * y2 * x2 + (x2 - x) / (x - x1) * y1 * y1 * x1)
          / x);
      };
    }

    fctInter2Points InterLin()
    {
      return[](double x, double x1, double x2, double y1, double y2) -> double
      {
        return (y2 - y1) / (x2 - x1) * (x - x1) + y1;
      };
    }
  }
}