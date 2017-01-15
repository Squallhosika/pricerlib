#pragma once

#include <pricer/env/SPricerOption.h>
#include <pricer/UType.h>
#include <pricer/instrs/CPayoff.h>
#include <pricer/processes/IProcess.h>
#include <pricer/env/SPricerOption.h>
#include <pricer/math/CGrid.h>

namespace Pricer
{
  class CBatch
  {
  public:
    CBatch(SPricerOption p_opt, const ptr<IProcess>& p_process);
    double GetPriceCall(double p_S0, double p_mat, double p_strike);
    double GetPricePut(double p_S0, double p_mat, double p_strike);
    double virtual PriceEuropean(double p_S0, double p_mat, 
      const ptr<CPayoffEUropean>& p_spPayoff) = 0;

    void FillGridVIBach(double p_S0, CGridWithValue& p_grid);
  protected:
    SPricerOption m_opt;
    ptr<IProcess> m_spProcess;
    ptr<timeSteps> m_spTS;
  };
}