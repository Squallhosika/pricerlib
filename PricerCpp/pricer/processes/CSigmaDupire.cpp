#include <pricer/stdafx.h>
#include "CSigmaDupire.h"
#include <pricer/tools/UComparator.h>

namespace Pricer
{
  CSigmaDupire::CSigmaDupire(const ptr<CVol>& p_vol)
    : m_vol(p_vol)
  {
  }

  void CSigmaDupire::Init(const ptr<CGrid>& p_grid)
  {
    std::vector<double> l_oTimeSteps = p_grid->TimeSteps();
    std::vector<double> l_oSpaceSteps = p_grid->SpaceSteps();

    // TODO check for the efficiance in memory: The vector is deplace
    // lot of time or not
    DMatrix l_oLocVolVals;
    l_oLocVolVals.push_back(std::vector<double>(l_oSpaceSteps));
    l_oLocVolVals.reserve(l_oTimeSteps.size());

    for (auto iter = ++(l_oTimeSteps.begin()); iter != l_oTimeSteps.end(); iter++)
    {
      std::vector<double> l_oSmile;
      l_oSmile.reserve(l_oSpaceSteps.size());
      for (double sstep : l_oSpaceSteps)
        l_oSmile.push_back(evalLocalVol(*iter, sstep));

      l_oLocVolVals.push_back(l_oSmile);
    }
    l_oLocVolVals.front() = *(++l_oLocVolVals.begin());

    m_inter = std::make_shared<CSplineInY_interp2d>(l_oTimeSteps, l_oSpaceSteps,
      l_oLocVolVals, UFctInter::InterLin());
  }

  double CSigmaDupire::operator()(double p_t, double p_s) const
  {
    return m_inter->interp(p_t, p_s);
  }

  double CSigmaDupire::evalLocalVol(double p_tenor, double p_strike)
  {
    // DEBUGTODO remove the following variable :
    double derivInT = m_vol->Derivator()->DerivCallBachInT(p_tenor, p_strike);
    double derivInK2 = m_vol->Derivator()->DerivCallBachInK2(p_tenor, p_strike);

    // TODO find a better way to undle strange case
    // at the beginning of the diffusion
    if ( UComparator::AlmostEqual(derivInK2, 0.0))
      return m_vol->GetPoint(p_tenor, p_strike);

    return std::sqrt(2.0 * derivInT / derivInK2);
  }
}