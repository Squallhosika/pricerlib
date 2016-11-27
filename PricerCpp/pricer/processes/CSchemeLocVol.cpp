#include <pricer/stdafx.h>
#include "CSchemeLocVol.h"

namespace Pricer
{
  CSchemeLocVol::CSchemeLocVol(const ptr<CSigmaLoc>& p_sigma)
    :m_sigma(p_sigma)
  {
  }

  ptr<CSchemeLocVol> CSchemeLocVol::Factory(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein)
  {
    if (p_isMillstein)
      return std::static_pointer_cast<CSchemeLocVol>(
      std::make_shared<CSchemeLocVolMilstein>(p_sigma));

    return std::static_pointer_cast<CSchemeLocVol>(
      std::make_shared<CSchemeLocVolEuler>(p_sigma));
  }

  double CSchemeLocVolEuler::evol(double p_t, double p_s, double p_step, double p_w) const
  {
    return p_s + (*m_sigma)(p_t, p_s)* std::sqrt(p_step) * p_w;
  }

  double CSchemeLocVolMilstein::evol(double p_t, double p_s, double p_step, double p_w) const
  {
    double l_sigma = (*m_sigma)(p_t, p_s);
    return p_s + l_sigma * std::sqrt(p_step) * p_w
      + 0.5 * l_sigma * m_sigma->DerivInS(p_t, p_s) * p_step * (p_w  * p_w - 1.0);
  }
}