#include <pricer/stdafx.h>
#include "CVolSquare.h"

#include <pricer/math/CMatrix.h>

#include <algorithm>

namespace Pricer
{
  CVolSquare::CVolSquare(double p_spot, const vectd& p_oTenors,
    const vectd p_oStrikes, const CMatrix& p_oPoints)
    : CVol(p_spot), m_oTenors(p_oTenors), m_oStrikes(p_oStrikes), m_oPoints(p_oPoints),
    m_oInter(p_oTenors, p_oStrikes, p_oPoints, UFctInter::InterInVar(), false)
  {
    m_spDerivator = std::make_shared<CVolDerivSimple>(this);
  }

  CVolSquare::CVolSquare(double p_spot, const CGridWithValue& p_spGrid)
    : CVolSquare(p_spot, p_spGrid.TimeSteps(), p_spGrid.SpaceSteps(), p_spGrid.Values())
  {
  }

  double CVolSquare::GetPoint(double tenor, double strike) const
  {
    return const_cast<CSplineInY_interp2d&>(m_oInter).interp(tenor, strike);
  }

  ptr<CGridWithValue> CVolSquare::GridToMath()
  {
    return std::make_shared<CGridWithValue>(m_oTenors, m_oStrikes);
  }
}