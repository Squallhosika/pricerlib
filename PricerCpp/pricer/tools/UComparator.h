#ifndef pricer_tools_ucomparator
#define pricer_tools_ucomparator

namespace Pricer
{
  namespace UComparator
  {
    const double COMPAR_ERROR = std::numeric_limits<double>::epsilon() / 100;

    bool AlmostEqual(double x, double y);
  }
}

#endif
