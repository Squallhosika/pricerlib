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
    CMatrix l_oLocVolVals;
    l_oLocVolVals.push_back(std::vector<double>(l_oSpaceSteps));
    l_oLocVolVals.reserve(l_oTimeSteps.size());

    for (auto iter = ++(l_oTimeSteps.begin()); iter != l_oTimeSteps.end(); iter++)
    {
      std::vector<double> l_oSmile;
      l_oSmile.reserve(l_oSpaceSteps.size());
      //for (double sstep : l_oSpaceSteps)
      for (size_t i = 0; i < l_oSpaceSteps.size(); i++)
        l_oSmile.push_back(evalLocalVol(*iter, l_oSpaceSteps[i])); // sstep));

      l_oLocVolVals.push_back(l_oSmile);
    }
    l_oLocVolVals.front() = *(++l_oLocVolVals.begin());

    m_inter = std::make_shared<CSplineInY_interp2d>(l_oTimeSteps, l_oSpaceSteps,
      l_oLocVolVals, UFctInter::InterLin(), true);
  }

  double CSigmaDupire::operator()(double p_t, double p_s) const
  {
    return m_inter->interp(p_t, p_s);
  }

  CSigmaDupireBach::CSigmaDupireBach(const ptr<CVol>& p_vol)
    : CSigmaDupire(p_vol)
  {
  }

  double CSigmaDupireBach::evalLocalVol(double p_tenor, double p_strike)
  {
    // DEBUGTODO remove the following variable :
    double derivInT = m_vol->Derivator()->DerivCallBachInT(p_tenor, p_strike);
    double derivInK2 = m_vol->Derivator()->DerivCallBachInK2(p_tenor, p_strike);

    // TODO find a better way to undle strange case
    // at the beginning of the diffusion
    if (derivInK2 < 0.0 || derivInT < 0.0 || std::abs(derivInT) < 1e-5
      || std::abs(derivInK2) < 1e-5)
      return m_vol->GetPoint(p_tenor, p_strike);

    // DEBUGTODO
    double l_res = std::sqrt(2.0 * derivInT / derivInK2);

    return l_res;
  }

  double CSigmaDupireBS::evalLocalVol(double p_tenor, double p_strike)
  {
    // DEBUGTODO remove the following variable :
    double derivInT = m_vol->Derivator()->DerivCallBSInT(p_tenor, p_strike);
    double derivInK2 = m_vol->Derivator()->DerivCallBSInK2(p_tenor, p_strike);

    // TODO find a better way to undle strange case
    // at the beginning of the diffusion
    if (derivInK2 < 0.0 || derivInT < 0.0 || std::abs(derivInT) < 1e-5
      || std::abs(derivInK2 * p_strike * p_strike) < 1e-5)
      return m_vol->GetPoint(p_tenor, p_strike);

    // DEBUGTODO
    double l_res = std::sqrt(2.0 * derivInT / (p_strike * p_strike * derivInK2));

    return l_res;
  }
}