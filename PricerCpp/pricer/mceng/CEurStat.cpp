#include <pricer/stdafx.h>
#include "CEurStat.h"

#include <pricer/tools/UComparator.h>

namespace Pricer
{
  CStat::CStat(const ptr<timeSteps>& p_spTimeSteps, size_t p_iNbSimu)
    : m_spTimeSteps(p_spTimeSteps)
  {
    m_paths.reserve(p_iNbSimu);
  }

  void CStat::Add(ptr<path> p_spPath)
  {
    m_paths.push_back(p_spPath);
  }

  CStatEurPriceOnly::CStatEurPriceOnly(const ptr<CPayoffEUropean>& p_payoff,
    const ptr<timeSteps>& p_spTimeSteps, size_t p_iNbSimu)
    : CStat(p_spTimeSteps, p_iNbSimu), m_payoff(p_payoff)
  {
  }

  double CStatEurPriceOnly::Price() const
  {
    double l_dRes = 0.0;

    for (ptr<path> l_spPath : m_paths)
      l_dRes += m_payoff->operator()(l_spPath->back());

    return l_dRes / m_paths.size();
  }

  CStatPathDepend::CStatPathDepend(const ptr<CPayoffPathDepPartial>& p_payoff,
    const ptr<timeSteps>& p_spTimeSteps, size_t p_iNbSimu)
    : CStat(p_spTimeSteps, p_iNbSimu), m_payoff(p_payoff)
  {
  }

  double CStatPathDepend::Price() const
  {
    double l_dRes = 0.0;
    std::vector<size_t> l_oPoints = pointToSelect();
    
    for (ptr<path> l_spPath : m_paths)
    {
      std::vector<double> l_oValues;
      l_oValues.reserve(l_oPoints.size());
      for (size_t i = 0; i < l_oPoints.size(); i++)
        l_oValues.push_back(l_spPath->at(l_oPoints[i]));

      l_dRes += m_payoff->operator()(l_oValues);
    }

    return l_dRes / m_paths.size();
  }

  // Make the assumption that the Payoff timSetps 
  // and the stat TImeeStet are sorted
  std::vector<size_t> CStatPathDepend::pointToSelect() const
  {
    std::vector<size_t> l_oRes;
    l_oRes.reserve(m_payoff->NeededTS().size());
    // TODO comapare in speed with the version we use iterator for both
    // redefine the below line and lets test after some change
    //timeSteps::const_iterator l_pPathTSIter = m_spTimeSteps->begin();
    timeSteps::const_iterator l_pPayoffTSIter = m_payoff->NeededTS().begin();
    size_t l_iIter = 0;

    while (l_iIter != m_spTimeSteps->size()
      && l_pPayoffTSIter != m_payoff->NeededTS().end())
    {
      if (UComparator::AlmostEqual(m_spTimeSteps->at(l_iIter), *l_pPayoffTSIter))
      {
        l_oRes.push_back(l_iIter);
        l_iIter++;
        l_pPayoffTSIter++;
      }
      else if (m_spTimeSteps->at(l_iIter) < *l_pPayoffTSIter)
      {
        l_iIter++;
      }
      else
      {
        std::exception("Granularity of the simulation is not enouth in regard of the payoff");
      }
    }
    return l_oRes;
  }
}