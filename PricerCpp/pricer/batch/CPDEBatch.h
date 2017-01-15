#include "CBatch.h"

namespace Pricer
{
  class CPDEBatch : public CBatch
  {
  public:
    CPDEBatch(SPricerOption p_opt, const ptr<IProcess>& p_process);
    double virtual PriceEuropean(double p_S0, double p_mat,
      const ptr<CPayoffEUropean>& p_spPayoff);

  private:
        
  };
}