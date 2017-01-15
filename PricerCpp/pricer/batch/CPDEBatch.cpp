#include <pricer/stdafx.h>
#include "CPDEBatch.h"
#include <pricer/pdeeng/CPdeEng.h>

namespace Pricer
{
  CPDEBatch::CPDEBatch(SPricerOption p_opt, const ptr<IProcess>& p_process)
    : CBatch(p_opt, p_process)
  {}

  double CPDEBatch::PriceEuropean(double p_S0, double p_mat, const ptr<CPayoffEUropean>& p_spPayoff)
  {
    // TODO do that when we add make a better pde pricer for bach and BS local vol
    //ptr<timeSteps> m_spTS = std::make_shared<timeSteps>(m_opt.PDEsizeInT, p_mat / m_opt.PDEsizeInT);

    //CPdeEng l_pdeEng(m_spProcess, p_spPayoff, *m_spTS, p_S0 - 50.0, p_S0 + 50.0, m_opt.PDEtheta);
    //l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
    //std::cout << "m_underMin: " << l_pdeEng.m_underMin << std::endl;
    //std::cout << "m_underMax: " << l_pdeEng.m_underMax << std::endl;
    //l_pdeEng.InitH(opt.PDEsizeInH);
    //std::vector<double> l_res = l_pdeEng.Compute();
    //Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

    //return l_stat.Price();
    throw new NotImplementedException();
    return INFINITY;
  }
}