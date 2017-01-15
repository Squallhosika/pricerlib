#ifndef pricer_cprocessvolloc
#define pricer_cprocessvolloc

#include "IProcess.h"
#include <algorithm>
#include "CSchemeLocVol.h"

namespace Pricer
{
  template<typename SCHEME>
  class CProcessVolLoc : public IProcess
  {
  public:
    CProcessVolLoc(const ptr<CSigmaLoc>& p_sigma);
    virtual ~CProcessVolLoc() {};

    double evol(double p_t, double p_s, double p_step, double p_w) const;
    virtual double StDev(double p_t, double p_s, double p_step) const = 0;
    ptr<CSigmaLoc> Sigma() const;

  protected:
    ptr<SCHEME> m_spScheme;
    ptr<CSigmaLoc> m_sigma;
  };

  class CProcessVolLocBach : public CProcessVolLoc<CSchemeLocVolBach>
  {
  public:
    CProcessVolLocBach(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein = false);
    double StDev(double p_t, double p_s, double p_step) const;
  };

  class CProcessVolLocBS : public CProcessVolLoc<CSchemeLocVolBS>
  {
  public:
    CProcessVolLocBS(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein = false);
    double StDev(double p_t, double p_s, double p_step) const;
  };
}
#endif
