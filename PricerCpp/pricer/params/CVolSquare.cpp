#include <pricer/stdafx.h>
#include "CVolSquare.h"

#include <algorithm>

namespace Pricer
{
  CVolSquare::CVolSquare(const std::vector<double>& p_oTenors,
    const std::vector<double> p_oStrikes, const DMatrix& p_oPoints)
    : m_oTenors(p_oTenors), m_oStrikes(p_oStrikes), m_oPoints(p_oPoints),
    m_oInter(p_oTenors, p_oStrikes, p_oPoints)
  {
  }

  double CVolSquare::GetPoint(double tenor, double strike) const
  {
    return const_cast<CVar_spline_interp2d&>(m_oInter).interp(tenor, strike);
  }
}