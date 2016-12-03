#ifndef pricer_csigmadupire
#define pricer_csigmadupire

#include "pricer/UType.h"
#include "CSigmaLoc.h"
#include <pricer/params/CVol.h>
#include <pricer/math/Cinterp2d.h>

namespace Pricer
{
  class CSigmaDupire : public CSigmaLoc
  {
  public:
    CSigmaDupire();
    ~CSigmaDupire() {};
    void Init(const ptr<CVol>& p_vol);
    double virtual operator()(double p_t, double p_s) const;

  private:
    ptr<CBase_interp2d> m_inter;
  };
}

#endif