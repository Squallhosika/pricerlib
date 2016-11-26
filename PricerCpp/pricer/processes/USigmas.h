#ifndef pricer_usigma
#define pricer_usigma

#include <algorithm>

namespace Pricer
{
  class CSigmaConst
  {
  public:
    CSigmaConst(double p_sig);

    double operator()(double) const;
    double Deriv(double) const;

  private:
    double m_sig;
  };

  class CSigmaCEV
  {
  public:
    CSigmaCEV(double p_sig0, double p_beta);

    double operator()(double p_s) const;
    double Deriv(double p_s) const;

  private:
    double m_sig0;
    double m_beta;
  };

  class CSigmaLNShifted
  {
  public:
    CSigmaLNShifted(double p_sig0, double p_shift);

    double operator()(double p_s) const;
    double Deriv(double p_s) const;

  private:
    double m_sig0;
    double m_shift;
  };

  class CSigmaConvex
  {
  public:
    CSigmaConvex(double p_sig0, double p_alpha,  double p_beta, double p_atm);

    double operator()(double p_s) const;
    double Deriv(double p_s) const;

  private:
    double m_sig0;
    double m_alpha;
    double m_beta;
    double m_atm;
  };
}


#endif