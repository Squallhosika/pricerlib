#ifndef pricer_upayoffs
#define pricer_upayoffs

#include <algorithm>

namespace Pricer
{
  class CPayoffId
  {
  public:
    double operator()(double p_s) const;
  };

  class CPayoffConst
  {
  public:
    CPayoffConst(double p_const);
    double operator()(double) const;

  private:
    double m_const;
  };

  class CPayoffCall
  {
  public:
    CPayoffCall(double p_k);
    double operator()(double p_s) const;

  private:
    double m_k;
  };

  class CPayoffDigit
  {
  public:
    CPayoffDigit(double p_k);
    double operator()(double p_s) const;

  private:
    double m_k;
  };
}


#endif