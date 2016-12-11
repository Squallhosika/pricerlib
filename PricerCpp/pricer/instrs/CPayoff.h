#ifndef pricer_upayoffs
#define pricer_upayoffs

#include <algorithm>
#include <pricer/UType.h>

namespace Pricer
{

  class CPayoff
  {
  public:
    virtual ~CPayoff();
  };

  class CPayoffEUropean : public CPayoff
  {
  public:
    virtual ~CPayoffEUropean();
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

  // Class of path dependent payoff that 
  // will always require all the path
  // exemple barier ...
  class CPayoffPathDepFull : public CPayoff
  {
  public:
    CPayoffPathDepFull();
    virtual ~CPayoffPathDepFull();
  };

  // Require only certain value in the
  // path
  class CPayoffPathDepPartial : public CPayoff
  {
  public:
    CPayoffPathDepPartial(const timeSteps& p_oNeededTS);
    virtual ~CPayoffPathDepPartial();
    double operator()(const std::vector<double>& p_s) const;
    const timeSteps& NeededTS() const;

  protected:
    timeSteps m_spNeededTS;
    virtual double Value(const std::vector<double>& p_s) const = 0;
  };

  // Not a VolSwap but just a vol because return 
  // only volatility use to price the ATM VolSwap
  class CPayoffVolBach : public CPayoffPathDepPartial
  {
  public:
    CPayoffVolBach(const timeSteps& p_oNeededTS);

  protected:
    virtual double Value(const std::vector<double>& p_s) const;
  };

  // Not a VarSwap but just a var because return 
  // only volatility use to price the ATM VolSwap
  class CPayoffVarBach : public CPayoffPathDepPartial
  {
  public:
    CPayoffVarBach(const timeSteps& p_oNeededTS);

  protected:
    virtual double Value(const std::vector<double>& p_s) const;
  };

}


#endif