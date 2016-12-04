#ifndef pricer_params_cvolderivsimple
#define pricer_params_cvolderivsimple

#include "IVolDerivator.h"

namespace Pricer
{
  class CVolDerivSimple : public IVolDerivator
  {
  public:
    CVolDerivSimple(CVol* p_vol);
    ~CVolDerivSimple();

    double virtual DerivInT(double T, double K) const;
    double virtual DerivInK(double T, double K) const;
    double virtual DerivInK2(double T, double K) const;

    double virtual DerivCallBachInT(double T, double K) const;
    double virtual DerivCallBachInK(double T, double K) const;
    double virtual DerivCallBachInK2(double T, double K) const;

  private:
    // build inside the vol object so should not own the
    // object memory deleting
    // TODO should we use a weak pointer here instead of
    // a raw pointer
    CVol* m_vol;
    double m_tenorEpsi = 0.01;
    double m_strikeEpsi = 0.01;
  };
}

#endif