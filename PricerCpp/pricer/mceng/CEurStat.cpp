#include <pricer/stdafx.h>
#include "CEurStat.h"

namespace Pricer
{
  CStatEurPriceOnly::CStatEurPriceOnly(const ptr<CPayoffEUropean>& p_payoff,
    const ptr<timeSteps>& p_spTimeSteps, size_t p_nbSimu)
    : m_payoff(p_payoff), m_spTimeSteps(p_spTimeSteps)
  {
    m_paths.reserve(p_nbSimu);
  }

  double CStatEurPriceOnly::Price() const
  {
    double m_result = 0.0;

    for (ptr<path> l_spPath : m_paths)
      m_result += m_payoff->operator()(l_spPath->back());

    return m_result / m_paths.size();
  }

  void CStatEurPriceOnly::Add(ptr<path> p_spPath)
  {
    m_paths.push_back(p_spPath);
  }
}