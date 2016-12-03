#ifndef pricer_upayoffs
#define pricer_upayoffs

#include <algorithm>

namespace Pricer
{
  class CPayoffEUropean
  {
  public:
    virtual ~CPayoffEUropean() {};
    double virtual operator()(double p_s) const = 0;
  };

  class CPayoffId : public CPayoffEUropean
  {
  public:
    double virtual operator()(double p_s) const;
  };

  class CPayoffConst : public CPayoffEUropean
  {
  public:
    CPayoffConst(double p_const);
    double operator()(double) const;

  private:
    double m_const;
  };

  class CPayoffCall : public CPayoffEUropean
  {
  public:
    CPayoffCall(double p_k);
    double operator()(double p_s) const;

  private:
    double m_k;
  };

  class CPayoffDigit : public CPayoffEUropean
  {
  public:
    CPayoffDigit(double p_k);
    double operator()(double p_s) const;

  private:
    double m_k;
  };
}


#endif