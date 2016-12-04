#include <pricer/stdafx.h>
#include "CVolSquare.h"

#include <algorithm>

namespace Pricer
{
  CVolSquare::CVolSquare(double p_spot, const std::vector<double>& p_oTenors,
    const std::vector<double> p_oStrikes, const DMatrix& p_oPoints)
    : CVol(p_spot), m_oTenors(p_oTenors), m_oStrikes(p_oStrikes), m_oPoints(p_oPoints),
    m_oInter(p_oTenors, p_oStrikes, p_oPoints, UFctInter::InterInVar())
  {
    m_spDerivator = std::make_shared<CVolDerivSimple>(this);
  }

  double CVolSquare::GetPoint(double tenor, double strike) const
  {
    return const_cast<CSplineInY_interp2d&>(m_oInter).interp(tenor, strike);
  }
}