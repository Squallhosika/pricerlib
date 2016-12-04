#ifndef pricer_math_cinterp2d
#define pricer_math_cinterp2d

// Come from the numerical recipes
#include <vector>
#include <algorithm>
#include "CInterp.h"
#include <pricer/UType.h>
#include "UFctInter.h"

namespace Pricer
{
  class CBase_interp2d
  {
  public:
    ~CBase_interp2d() {};
    double virtual interp(double x, double y) = 0;
  };
  
  // TODO uncomment that and make CSplineInY_interp2d 
  // or CVar_spline_interp2d heritate from
  // this class do we have to keep CSplineInY_interp2d

  //class CYThenX_interp2d : CBase_interp2d
  //{
  //public:
  //  CYThenX_interp2d(const std::vector<double>& x, const std::vector<double>& y,
  //  const DMatrix& z);
  //  ~CYThenX_interp2d() {};

  //private:
  //  std::vector<double> m_oX;
  //  std::vector<double> m_oY;
  //  DMatrix m_oZ;
  //  std::vector<ptr<Base_interp>> m_inters;
  //};

  // interpolation in variance for x
  // and in spline for y
  class CSplineInY_interp2d : public CBase_interp2d
  {
  public:
    CSplineInY_interp2d(const std::vector<double>& x, const std::vector<double>& y, const DMatrix& z, 
      const UFctInter::fctInter2Points& p_interInX);

    double virtual interp(double x, double y);

  protected:
    std::vector<double> m_oX;
    std::vector<double> m_oY;
    DMatrix m_oZ;
    // TODO replace Spline_interp by 
    // Base interp to be more generic
    std::vector<Spline_interp> m_oSplines;
    UFctInter::fctInter2Points m_interInX;
    bool m_bInfomExtrapol = true;

    void setSpline();
    int locateX(double x);
  };

  // TODO remove when sure that We do not use it anymore
  //class CVar_spline_interp2d : public CSplineInY_interp2d
  //{
  //public:
  //  CVar_spline_interp2d(const std::vector<double>& x, const std::vector<double>& y,
  //  const DMatrix& z);

  //  double virtual interp(double x, double y);

  //private:
  //  std::vector<double> m_oX;
  //  std::vector<double> m_oY;
  //  DMatrix m_oZ;
  //  // TODO replace Spline_interp by 
  //  // Base interp to be more generic
  //  std::vector<Spline_interp> m_oSplines;
  //};
}

#endif