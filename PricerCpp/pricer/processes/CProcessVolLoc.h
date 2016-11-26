#ifndef pricer_cprocessvolloc
#define pricer_cprocessvolloc

#include "IProcess.h"
#include <algorithm>
#include "CSchemeLocVol.h"

namespace Pricer
{
  template<typename SIGMA>
  class CProcessVolLoc : public IProcess
  {
  public:
    CProcessVolLoc(const SIGMA& p_sigma, bool p_isMillstein = false);
    //~CProcessVolLoc();

    double evol(double p_s, double p_step, double p_w) const;
    inline double StDev(double p_s, double p_step) const;
    inline SIGMA Sigma() const;

  private:
    SIGMA m_sigma;
    ptr<CSchemeLocVol<SIGMA>> m_spScheme;

  };

  template<typename SIGMA>
  CProcessVolLoc<SIGMA>::CProcessVolLoc(const SIGMA& p_sigma, bool p_isMillstein)
    : m_sigma(p_sigma)
  {
    m_spScheme = CSchemeLocVol<SIGMA>::Factory(p_sigma, p_isMillstein);
  }

  template<typename SIGMA>
  double CProcessVolLoc<SIGMA>::evol(double p_s, double p_step, double p_w) const
  {
    return m_spScheme->evol(p_s, p_step, p_w);
    // return p_s + m_sigma(p_s) * std::sqrt(p_step) * p_w;
  }

  // TODO Not the StDev of the processus
  // have to be improve
  template<typename SIGMA>
  inline double CProcessVolLoc<SIGMA>::StDev(double p_s, double p_step) const
  {
    return m_sigma(p_s) * std::sqrt(p_step);
  }

  template<typename SIGMA>
  inline SIGMA CProcessVolLoc<SIGMA>::Sigma() const
  {
    return m_sigma;
  }
}
#endif
