#include <pricer/stdafx.h>
#include "CInterp2d.h"

#include <algorithm>
#include <iostream>

namespace Pricer
{
  CSplineInY_interp2d::CSplineInY_interp2d(const std::vector<double>& x,
    const std::vector<double>& y, const DMatrix& z, 
    const UFctInter::fctInter2Points& p_interInX,
    bool p_bInfomExtrapol)
    : m_oX(x), m_oY(y), m_oZ(z), m_interInX(p_interInX),
    m_bInfomExtrapol(p_bInfomExtrapol)
  {
    m_oSplines.reserve(y.size());
    setSpline();
  }

  void CSplineInY_interp2d::setSpline()
  {
    // TODO add the spline option when we set it
    // Look at the Spline_interp constructor for
    // more information
    for (size_t i = 0; i < m_oX.size(); i++)
      m_oSplines.push_back(Spline_interp(m_oY, m_oZ[i]));
  }

  int CSplineInY_interp2d::locateX(double x)
  {
    int ju, jm, jl;
    if (m_oX.size() < 2) throw("locate size error");
    bool ascnd = (m_oX.back() >= m_oX.front());
    jl = 0;
    ju = m_oX.size() - 1;
    while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if (x >= m_oX[jm] == ascnd)
        jl = jm;
      else
        ju = jm;
    }
    return std::max(0, std::min((int)m_oX.size() - 2, jl));
  }

  double CSplineInY_interp2d::interp(double x, double y)
  {
    if (m_bInfomExtrapol)
    {
      if (x < m_oX.front() || x > m_oX.back())
        std::cout << "extrapol in x: x=" << x << " not in [" << m_oX.front() << ", " 
        << m_oX.back() << "]" << std::endl;

      if (y < m_oY.front() || y > m_oY.back())
        std::cout << "extrapol in y: y=" << y << " not in [" << m_oY.front() << ", "
        << m_oY.back() << "]" << std::endl;
    }

    int l_jl = locateX(x);
    double l_x1 = m_oX[l_jl];
    double l_x2 = m_oX[l_jl + 1];
    double l_z1 = m_oSplines[l_jl].interp(y);
    double l_z2 = m_oSplines[l_jl + 1].interp(y);

    return m_interInX(x, l_x1, l_x2, l_z1, l_z2);
  }

  //double CVar_spline_interp2d::interp(double x, double y)
  //{
  //  int l_jl = locateX(x);
  //  double l_x1 = m_oX[l_jl];
  //  double l_x2 = m_oX[l_jl + 1];
  //  double l_z1 = m_oSplines[l_jl].interp(y);
  //  double l_z2 = m_oSplines[l_jl + 1].interp(y);

  //  return std::sqrt((l_z2 * l_z2 * l_x2 + (l_x2 - x) / (x - l_x1) * l_z1 * l_z1 * l_x1)
  //    / ((l_x2 - x) / (x - l_x1) + 1.) / x);
  //}


}