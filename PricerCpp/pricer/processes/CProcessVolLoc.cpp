#include <pricer/stdafx.h>
#include "CProcessVolLoc.h"

namespace Pricer
{
  CProcessVolLoc::CProcessVolLoc(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein)
    : m_sigma(p_sigma)
  {
    m_spScheme = CSchemeLocVol::Factory(p_sigma, p_isMillstein);
  }

  double CProcessVolLoc::evol(double p_t, double p_s, double p_step, double p_w) const
  {
    return m_spScheme->evol(p_t, p_s, p_step, p_w);
    // return p_s + m_sigma(p_s) * std::sqrt(p_step) * p_w;
  }

  // TODO Not the StDev of the processus
  // have to be improve
  double CProcessVolLoc::StDev(double p_t, double p_s, double p_step) const
  {
    return (*m_sigma)(p_t, p_s)* std::sqrt(p_step);
  }
  
  ptr<CSigmaLoc> CProcessVolLoc::Sigma() const
  {
    return m_sigma;
  }
}