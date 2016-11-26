#ifndef pricer_ceurstat
#define pricer_ceurstat

#include <vector>

#include "pricer/UType.h"

namespace Pricer
{
  // TODO add interval error price
  // info before full pricing
  template<typename PAYOFF>
  class CStatEurPriceOnly
  {
  public:
    typedef std::vector<double>    timeSteps;
    typedef std::vector<double>    path;
    typedef std::vector<ptr<path>> paths;
    
    CStatEurPriceOnly(const PAYOFF& p_payoff, const ptr<timeSteps>& p_spTimeSteps, size_t p_nbSimu = 1000);
    //~CStatEurPriceOnly();

    double Price() const;
    void Add(ptr<path> p_spPath);

  private:
    PAYOFF         m_payoff;
    ptr<timeSteps> m_spTimeSteps;
    paths          m_paths;
  };

  template<typename PAYOFF>
  CStatEurPriceOnly<PAYOFF>::CStatEurPriceOnly(const PAYOFF& p_payoff, const ptr<timeSteps>& p_spTimeSteps, size_t p_nbSimu = 1000)
    : m_payoff(p_payoff), m_spTimeSteps(p_spTimeSteps)
  {
    m_paths.reserve(p_nbSimu);
  }

  template<typename PAYOFF>
  double CStatEurPriceOnly<PAYOFF>::Price() const
  {
    double m_result = 0.0;

    for (ptr<path> l_spPath : m_paths)
      m_result += m_payoff(l_spPath->back());

    return m_result / m_paths.size();
  }

  template<typename PAYOFF>
  void CStatEurPriceOnly<PAYOFF>::Add(ptr<path> p_spPath)
  {
    m_paths.push_back(p_spPath);
  }
}

#endif