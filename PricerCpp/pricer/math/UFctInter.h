#pragma once

// Come from the numerical recipes
#include <functional>

namespace Pricer
{
  namespace UFctInter
  {
    typedef std::function<double(double, double, double, double, double)> fctInter2Points;
    fctInter2Points InterInVar();
    fctInter2Points InterLin();
  }
}
