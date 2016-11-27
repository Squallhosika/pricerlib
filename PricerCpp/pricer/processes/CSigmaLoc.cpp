#include <pricer/stdafx.h>
#include "CSigmaLoc.h"

namespace Pricer
{
  double CSigmaLoc::Value(double p_t, double p_s) const
  {
    return (*this)(p_t, p_s);
  }

  double CSigmaLocS::operator()(double p_s) const
  {
    return (*this)(0.0, p_s);
  }

  double CSigmaLocS::DerivInS(double p_s) const
  {
    return this->DerivInS(0.0, p_s);
  }

  CSigmaConst::CSigmaConst(double p_sig)
    : m_sig(p_sig)
  {}

  double CSigmaConst::operator()(double) const
  {
    return m_sig;
  }

  double CSigmaConst::DerivInS(double) const
  {
    return 0.0;
  }

  CSigmaCEV::CSigmaCEV(double p_sig0, double p_beta)
    : m_sig0(p_sig0), m_beta(p_beta)
  {}

  double CSigmaCEV::operator()(double p_s) const
  {
    return m_sig0 * std::pow(std::abs(p_s), m_beta);
  }

  double CSigmaCEV::DerivInS(double p_s) const
  {
    return m_sig0 * m_beta * std::pow(std::abs(p_s), m_beta - 1.0);
  }

  CSigmaLNShifted::CSigmaLNShifted(double p_sig0, double p_shift)
    : m_sig0(p_sig0), m_shift(p_shift)
  {}

  double CSigmaLNShifted::operator()(double p_s) const
  {
    return m_sig0 * (p_s + m_shift);
  }

  double CSigmaLNShifted::DerivInS(double p_s) const
  {
    return m_sig0;
  }

  CSigmaConvex::CSigmaConvex(double p_sig0, double p_alpha, double p_beta, double p_atm)
    : m_sig0(p_sig0), m_alpha(p_alpha), m_beta(p_beta), m_atm(p_atm)
  {}

  double CSigmaConvex::operator()(double p_s) const
  {
    return m_sig0 + m_alpha * std::pow(std::abs(p_s - m_atm), 1 + m_beta);
  }

  double CSigmaConvex::DerivInS(double p_s) const
  {
    if (p_s < m_atm)
      return -1.0 * m_alpha * (1 + m_beta) * std::pow(std::abs(p_s - m_atm), m_beta);

    return m_alpha * (1 + m_beta) * std::pow(p_s - m_atm, m_beta);
  }

}