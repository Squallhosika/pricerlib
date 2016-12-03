#ifndef pricer_params_cvolderivsimple
#define pricer_params_cvolderivsimple

#include "IVolDerivator.h"

namespace Pricer
{
  class CVolDerivSimple : public IVolDerivator
  {
  public:
    CVolDerivSimple(const ptr<CVol>& p_vol);
    ~CVolDerivSimple();

    double DerivInT(double T, double K) const;
    double DerivInK(double T, double K) const;
    double DerivInK2(double T, double K) const;

    double DerivCallBachInT(double T, double K) const;
    double DerivCallBachInK(double T, double K) const;
    double DerivCallBachInK2(double T, double K) const;

  private:
    ptr<CVol> m_vol;
    double m_tenorEpsi = 0.01;
    double m_strikeEpsi = 0.01;
  };
}

#endif