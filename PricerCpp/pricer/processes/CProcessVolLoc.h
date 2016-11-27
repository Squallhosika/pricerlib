#ifndef pricer_cprocessvolloc
#define pricer_cprocessvolloc

#include "IProcess.h"
#include <algorithm>
#include "CSchemeLocVol.h"

namespace Pricer
{
  class CProcessVolLoc : public IProcess
  {
  public:
    CProcessVolLoc(const ptr<CSigmaLoc>& p_sigma, bool p_isMillstein = false);
    //~CProcessVolLoc();

    double evol(double p_t, double p_s, double p_step, double p_w) const;
    double StDev(double p_t, double p_s, double p_step) const;
    ptr<CSigmaLoc> Sigma() const;

  private:
    ptr<CSigmaLoc> m_sigma;
    ptr<CSchemeLocVol> m_spScheme;

  };
}
#endif
