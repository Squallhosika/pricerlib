#pragma once

// Come from the numerical recipes
#include <vector>
#include <algorithm>
#include "CInterp.h"
#include <pricer/UType.h>
#include "UFctInter.h"
#include <pricer/math/CMatrix.h>

namespace Pricer
{
  class CBase_interp2d
  {
  public:
    virtual ~CBase_interp2d();
    double virtual interp(double x, double y) = 0;
  };
  
  // TODO make sure that all the intelligence
  // of the Y then X is well give to
  // CYThenX_interp2d to CSplineInY_interp2d

  class CYThenX_interp2d : public CBase_interp2d
  {
  public:
    CYThenX_interp2d(const vectd& x, const vectd& y, const CMatrix& z);
    ~CYThenX_interp2d() {};

  protected:
    vectd m_oX;
    vectd m_oY;
    CMatrix m_oZ;
    std::vector<ptr<Base_interp>> m_inters;
  };

  // interpolation in variance for x
  // and in spline for y
  class CSplineInY_interp2d : public CYThenX_interp2d
  {
  public:
    CSplineInY_interp2d(const vectd& x, const vectd& y, const CMatrix& z,
      const UFctInter::fctInter2Points& p_interInX, bool p_bInfomExtrapol);

    double virtual interp(double x, double y);

  protected:
    // TODO replace Spline_interp by 
    // Base interp to be more generic
    std::vector<Spline_interp> m_oSplines;
    UFctInter::fctInter2Points m_interInX;
    bool m_bInfomExtrapol;

    void setSpline();
    int locateX(double x);
  };

  // TODO remove when sure that We do not use it anymore
  //class CVar_spline_interp2d : public CSplineInY_interp2d
  //{
  //public:
  //  CVar_spline_interp2d(const std::vector<double>& x, const std::vector<double>& y,
  //  const CMatrix& z);

  //  double virtual interp(double x, double y);

  //private:
  //  std::vector<double> m_oX;
  //  std::vector<double> m_oY;
  //  CMatrix m_oZ;
  //  // TODO replace Spline_interp by 
  //  // Base interp to be more generic
  //  std::vector<Spline_interp> m_oSplines;
  //};
}