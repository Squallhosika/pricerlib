#ifndef pricer_usigma
#define pricer_usigma

#include <algorithm>

namespace Pricer
{
  class CSigmaLoc
  {
  public:
    virtual ~CSigmaLoc();
    double virtual operator()(double p_t, double p_s) const = 0;
    double Value(double p_t, double p_s) const;
    double virtual DerivInS(double p_t, double p_s) const = 0;
  };

  class CSigmaLocS : public CSigmaLoc
  {
  public: 
    virtual ~CSigmaLocS();
    double virtual operator()(double p_t, double p_s) const;
    double virtual operator()(double p_s) const = 0;
    double virtual DerivInS(double p_t, double p_s) const;
    double virtual DerivInS(double p_s) const = 0;
  };

  class CSigmaConst
  {
  public:
    CSigmaConst(double p_sig);

    double operator()(double) const;
    double DerivInS(double) const;

  private:
    double m_sig;
  };

  class CSigmaCEV
  {
  public:
    CSigmaCEV(double p_sig0, double p_beta);

    double operator()(double p_s) const;
    double DerivInS(double p_s) const;

  private:
    double m_sig0;
    double m_beta;
  };

  class CSigmaLNShifted
  {
  public:
    CSigmaLNShifted(double p_sig0, double p_shift);

    double operator()(double p_s) const;
    double DerivInS(double p_s) const;

  private:
    double m_sig0;
    double m_shift;
  };

  class CSigmaConvex
  {
  public:
    CSigmaConvex(double p_sig0, double p_alpha,  double p_beta, double p_atm);

    double operator()(double p_s) const;
    double DerivInS(double p_s) const;

  private:
    double m_sig0;
    double m_alpha;
    double m_beta;
    double m_atm;
  };
}


#endif