#include <pricer/stdafx.h>
#include "CMCBatch.h"
#include <pricer/tools/UComparator.h>

namespace Pricer
{
  CMCBatch::CMCBatch(SPricerOption p_opt, const ptr<IProcess>& p_process)
    : CBatch(p_opt, p_process), m_dLastMat(HUGE_VALD)
  {}

  void CMCBatch::genRnd(double p_mat)
  {
    if (UComparator::AlmostEqual(p_mat, m_dLastMat))
      return;

    m_spTS = std::make_shared<timeSteps>(m_opt.MCnbTimeStep, p_mat / m_opt.MCnbTimeStep);
    m_spPathGen = std::make_shared<CPathGen<CNormGen>>(m_spProcess, CNormGen(), m_spTS, m_opt.MCnbSimu);
    m_spPathGen->GenSequence();
    m_dLastMat = p_mat;
  }

  double CMCBatch::PriceEuropean(double p_S0, double p_mat,
    const ptr<CPayoffEUropean>& p_spPayoff)
  {
    genRnd(p_mat);

    CStatEurPriceOnly  l_stat(p_spPayoff, m_spTS);
    m_spPathGen->FillStat(p_S0, l_stat);

    return l_stat.Price();
  }
}