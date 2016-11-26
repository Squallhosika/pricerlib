#include "stdafx.h"
#include "UPayoffs.h"

namespace Pricer
{
  double CPayoffId::operator()(double p_s) const
  {
    return p_s;
  }

  CPayoffConst::CPayoffConst(double p_const)
    : m_const(p_const)
  {
  }

  double CPayoffConst::operator()(double) const
  {
    return m_const;
  }

  CPayoffCall::CPayoffCall(double p_k)
    : m_k(p_k)
  {
  }

  double CPayoffCall::operator()(double p_s) const
  {
    return std::max(0.0, p_s - m_k);
  }

  CPayoffDigit::CPayoffDigit(double p_k)
    : m_k(p_k)
  {
  }

  double CPayoffDigit::operator()(double p_s) const
  {
    return p_s > m_k ? 1.0 : 0.0;
  }

}