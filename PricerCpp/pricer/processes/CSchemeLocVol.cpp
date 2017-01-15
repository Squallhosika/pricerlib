#include <pricer/stdafx.h>
#include "CSchemeLocVol.h"

namespace Pricer
{
  CSchemeLocalVol::CSchemeLocalVol(const ptr<CSigmaLoc>& p_sigma)
    :m_sigma(p_sigma)
  {
  }

  CSchemeLocVolBach::CSchemeLocVolBach(const ptr<CSigmaLoc>& p_sigma)
    : CSchemeLocalVol(p_sigma)
  {
  }

  ptr<CSchemeLocVolBach> CSchemeLocVolBach::Factory(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein)
  {
    if (p_isMillstein)
      return std::static_pointer_cast<CSchemeLocVolBach>(
      std::make_shared<CSchemeLocVolBachMilstein>(p_sigma));

    return std::static_pointer_cast<CSchemeLocVolBach>(
      std::make_shared<CSchemeLocVolBachEuler>(p_sigma));
  }

  double CSchemeLocVolBachEuler::evol(double p_t, double p_s, double p_step, double p_w) const
  {
    return p_s + m_sigma->Value(p_t, p_s)* std::sqrt(p_step) * p_w;
  }

  double CSchemeLocVolBachMilstein::evol(double p_t, double p_s, double p_step, double p_w) const
  {
    double l_sigma = m_sigma->Value(p_t, p_s);
    return p_s + l_sigma * std::sqrt(p_step) * p_w
      + 0.5 * l_sigma * m_sigma->DerivInS(p_t, p_s) * p_step * (p_w  * p_w - 1.0);
  }

  CSchemeLocVolBS::CSchemeLocVolBS(const ptr<CSigmaLoc>& p_sigma)
    : CSchemeLocalVol(p_sigma)
  {
  }

  ptr<CSchemeLocVolBS> CSchemeLocVolBS::Factory(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein)
  {
    if (p_isMillstein)
      return std::static_pointer_cast<CSchemeLocVolBS>(
      std::make_shared<CSchemeLocVolBSMilstein>(p_sigma));

    return std::static_pointer_cast<CSchemeLocVolBS>(
      std::make_shared<CSchemeLocVolBSEuler>(p_sigma));
  }

  // TODO Look for scheme in log(S_t)
  double CSchemeLocVolBSEuler::evol(double p_t, double p_s, double p_step, double p_w) const
  {
    // WARNTODO Do on a train check the formula
    return  p_s * (1.0 + m_sigma->Value(p_t, p_s)* std::sqrt(p_step) * p_w);
  }

  double CSchemeLocVolBSMilstein::evol(double p_t, double p_s, double p_step, double p_w) const
  {
    // WARNTODO Do on a train check the formula
    double l_sigma = m_sigma->Value(p_t, p_s);
    return p_s * (1.0 + l_sigma * std::sqrt(p_step) * p_w
      + 0.5 * l_sigma * m_sigma->DerivInS(p_t, p_s) * p_step * (p_w  * p_w - 1.0));
  }
}