#include <pricer/stdafx.h>
#include "CGrid.h"

#include <numeric>
#include <algorithm>

namespace Pricer
{
  CGrid::CGrid(const vectd& p_timeSteps, const vectd& p_oSpaceSteps)
    : m_oTimeSteps(p_timeSteps), m_oSpaceSteps(p_oSpaceSteps)
  {
    std::sort(m_oTimeSteps.begin(), m_oTimeSteps.end());
    std::sort(m_oSpaceSteps.begin(), m_oSpaceSteps.end());
  }

  CGrid::vectd CGrid::TimeSteps() const
  {
    return m_oTimeSteps;
  }

  CGrid::vectd CGrid::TimeStepsDiff() const
  {
    // TODO check for a beter way than fill all and erase
    // the first element. if a solution is found apply also
    // to the method SpaceStepsDiff()
    vectd l_res(m_oTimeSteps.size());
    std::adjacent_difference(m_oTimeSteps.begin(), m_oTimeSteps.end(), l_res.begin());
    l_res.erase(l_res.begin());
    return l_res;
  }
  CGrid::vectd CGrid::SpaceSteps() const
  {
    return m_oSpaceSteps;
  }

  CGrid::vectd CGrid::SpaceStepsDiff() const
  {
    vectd l_res(m_oSpaceSteps.size());
    std::adjacent_difference(m_oSpaceSteps.begin(), m_oSpaceSteps.end(), l_res.begin());
    l_res.erase(l_res.begin());
    return l_res;
  }
}
