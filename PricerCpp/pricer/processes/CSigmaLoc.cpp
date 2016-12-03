#include <pricer/stdafx.h>
#include "CSigmaLoc.h"

namespace Pricer
{
  double CSigmaLoc::Value(double p_t, double p_s) const
  {
    return (*this)(p_t, p_s);
  }

  double CSigmaLoc::DerivInS(double p_t, double p_s) const
  {
    return 0.5 * ((*this)(p_t, p_s + s_epsiForS) 
      - (*this)(p_t, p_s - s_epsiForS)) / s_epsiForS;
  }

  const double CSigmaLoc::s_epsiForS = 0.01;

  double CSigmaLocS::operator()(double p_t, double p_s) const
  {
    return (*this)(p_s);
  }

  double CSigmaLocS::DerivInS(double p_t, double p_s) const
  {
    return this->DerivInS(p_s);
  }

  // TODO check the consicuve call of this function in debug
  // particulary which operator() is call
  double CSigmaLocS::DerivInS(double p_s) const
  {
    return CSigmaLoc::DerivInS(0.0, p_s);
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