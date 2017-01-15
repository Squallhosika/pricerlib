#include <pricer/stdafx.h>
#include "CBatch.h"
#include <pricer/instrs/CPayoff.h>
#include <pricer/math/UVolInverser.h>

namespace Pricer
{
  CBatch::CBatch(SPricerOption p_opt, const ptr<IProcess>& p_process)
    : m_opt(p_opt), m_spProcess(p_process)
  {
  }

  double CBatch::GetPriceCall(double p_S0, double p_mat, double p_strike)
  {
    ptr<CPayoffEUropean> l_payoff(new CPayoffCall(p_strike));

    return PriceEuropean(p_S0, p_mat, l_payoff);
  }

  void CBatch::FillGridVIBach(double p_S0, CGridWithValue& p_grid)
  {
    std::function<double(double, double)> l_fct =
      [p_S0, this](double ts, double ss)
    {
      return UVolInverser::InversCallBach(p_S0, 0.0, ss, this->GetPriceCall(p_S0, ts, ss), 
        ts, UVolInverser::VolMin, UVolInverser::VolMin);
    };

    p_grid.Fill(l_fct);
  }
}