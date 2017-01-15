#pragma once

#include "CBatch.h"
#include <pricer/mceng/CEurStat.h>
#include <pricer/mceng/CPathGen.h>
#include <pricer/mceng/CGenerator.h>

namespace Pricer
{
  class CMCBatch : public CBatch
  {
  public:
    CMCBatch(SPricerOption p_opt, const ptr<IProcess>& p_process);
    double virtual PriceEuropean(double p_S0, double p_mat,
      const ptr<CPayoffEUropean>& p_spPayoff);
  
  private:
    ptr<CPathGen<CNormGen>> m_spPathGen;
    double m_dLastMat;

    void genRnd(double p_mat);
  };
}
