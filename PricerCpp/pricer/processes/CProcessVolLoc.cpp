#include <pricer/stdafx.h>
#include "CProcessVolLoc.h"

namespace Pricer
{
  template<typename SCHEME>
  CProcessVolLoc<SCHEME>::CProcessVolLoc(const ptr<CSigmaLoc>& p_sigma)
    : m_sigma(p_sigma)
  {
  }

  template<typename SCHEME>
  double CProcessVolLoc<SCHEME>::evol(double p_t, double p_s, double p_step, double p_w) const
  {
    return m_spScheme->evol(p_t, p_s, p_step, p_w);
  }
  
  template<typename SCHEME>
  ptr<CSigmaLoc> CProcessVolLoc<SCHEME>::Sigma() const
  {
    return m_sigma;
  }

  CProcessVolLocBach::CProcessVolLocBach(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein)
    : CProcessVolLoc<CSchemeLocVolBach>(p_sigma)
  {
    m_spScheme = CSchemeLocVolBach::Factory(p_sigma, p_isMillstein);
  }

  // TODO Not the StDev of the processus
  // have to be improve
  double CProcessVolLocBach::StDev(double p_t, double p_s, double p_step) const
  {
    return (*m_sigma)(p_t, p_s)* std::sqrt(p_step);
  }

  CProcessVolLocBS::CProcessVolLocBS(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein)
    : CProcessVolLoc<CSchemeLocVolBS>(p_sigma)
  {
      m_spScheme = CSchemeLocVolBS::Factory(p_sigma, p_isMillstein);
    }

  double CProcessVolLocBS::StDev(double p_t, double p_s, double p_step) const
  {
    // WARNTODO check the formula
    return std::sqrt(std::exp(0.5 * (*m_sigma)(p_t, p_s) * p_step));
  }
}