#include <pricer/stdafx.h>
#include "CGrid.h"

#include <numeric>
#include <algorithm>

namespace Pricer
{
  CGrid::CGrid(const vectd& p_timeSteps, const vectd& p_oSpaceSteps)
    : m_oTimeSteps(p_timeSteps), m_oSpaceSteps(p_oSpaceSteps)
  {
    sort();
  }

  CGrid::CGrid(double p_MinTS, double p_MaxTS, size_t p_nTS, 
    double p_MinSpaceStep, double p_MaxSpaceStep, size_t p_nSpaceStep)
    : m_oTimeSteps(p_nTS), m_oSpaceSteps(p_nSpaceStep)
  {
    std::generate(m_oTimeSteps.begin(), m_oTimeSteps.end(), 
      SGener_iter(p_MinTS, p_MaxTS, p_nTS));
    std::generate(m_oSpaceSteps.begin(), m_oSpaceSteps.end(),
      SGener_iter(p_MinSpaceStep, p_MaxSpaceStep, p_nSpaceStep));
    sort();
  }

  void CGrid::sort()
  {
    std::sort(m_oTimeSteps.begin(), m_oTimeSteps.end());
    std::sort(m_oSpaceSteps.begin(), m_oSpaceSteps.end());
  }

  vectd CGrid::TimeSteps() const
  {
    return m_oTimeSteps;
  }

  vectd CGrid::TimeStepsDiff() const
  {
    // TODO check for a beter way than fill all and erase
    // the first element. if a solution is found apply also
    // to the method SpaceStepsDiff()
    vectd l_res(m_oTimeSteps.size());
    std::adjacent_difference(m_oTimeSteps.begin(), m_oTimeSteps.end(), l_res.begin());
    l_res.erase(l_res.begin());
    return l_res;
  }
  vectd CGrid::SpaceSteps() const
  {
    return m_oSpaceSteps;
  }

  vectd CGrid::SpaceStepsDiff() const
  {
    vectd l_res(m_oSpaceSteps.size());
    std::adjacent_difference(m_oSpaceSteps.begin(), m_oSpaceSteps.end(), l_res.begin());
    l_res.erase(l_res.begin());
    return l_res;
  }

  CGridWithValue::CGridWithValue(const vectd& p_timeSteps, const vectd& p_oSpaceSteps)
    : CGrid(p_timeSteps, p_oSpaceSteps), 
    m_oValues(p_timeSteps.size(), p_oSpaceSteps.size())
  {
  }

  const CMatrix& CGridWithValue::Values() const
  {
    return m_oValues;
  }

  void CGridWithValue::SetValue(size_t p_iTS, size_t p_iSS, double value)
  {
    m_oValues[p_iTS][p_iSS] = value;
  }

  void CGridWithValue::Fill(const std::function<double(double, double)>& p_fctTimeSpace)
  {
    for (size_t iTS = 0; iTS < m_oTimeSteps.size(); iTS++)
    for (size_t iSS = 0; iSS < m_oSpaceSteps.size(); iSS++)
      SetValue(iTS, iSS, p_fctTimeSpace(m_oTimeSteps[iTS], m_oSpaceSteps[iSS]));
  }
}
