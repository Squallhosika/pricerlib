#ifndef pricer_math_ufunctioninter
#define pricer_math_ufunctioninter

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

#endif
