#ifndef pricer_ceurstat
#define pricer_ceurstat

#include <vector>

#include <pricer/UType.h>
#include <pricer/instrs/CPayoff.h>

namespace Pricer
{
  // TODO add interval error price
  // info before full pricing
  class CStatEurPriceOnly
  {
  public:
    typedef std::vector<double>    timeSteps;
    typedef std::vector<double>    path;
    typedef std::vector<ptr<path>> paths;
    
    CStatEurPriceOnly(const ptr<CPayoffEUropean>& p_payoff, 
      const ptr<timeSteps>& p_spTimeSteps, size_t p_nbSimu = 1000);
    //~CStatEurPriceOnly();

    double Price() const;
    void Add(ptr<path> p_spPath);

  private:
    ptr<CPayoffEUropean> m_payoff;
    ptr<timeSteps> m_spTimeSteps;
    paths          m_paths;
  };
}

#endif