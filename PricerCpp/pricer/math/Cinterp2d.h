#ifndef pricer_math_cinterp2d
#define pricer_math_cinterp2d

// Come from the numerical recipes
#include <vector>
#include <algorithm>
#include "CInterp.h"
#include <pricer/UType.h>

namespace Pricer
{
  class CBase_interp2d
  {
  public:
    double virtual interp(double x, double y) = 0;
  };

  // interpolation in variance for x
  // and in spline for y
  class CVar_spline_interp2d : public CBase_interp2d
  {
  public:
    CVar_spline_interp2d(const std::vector<double>& x, const std::vector<double>& y,
    const DMatrix& z);

    double virtual interp(double x, double y);

  private:
    std::vector<double> m_oX;
    std::vector<double> m_oY;
    DMatrix m_oZ;
    // TODO replace Spline_interp by 
    // Base interp to be more generic
    std::vector<Spline_interp> m_oSplines;

    void setSpline();
    int locateX(double x);
  };
}

#endif