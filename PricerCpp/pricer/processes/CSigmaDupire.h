#ifndef pricer_csigmadupire
#define pricer_csigmadupire

#include "pricer/UType.h"
#include "CSigmaLoc.h"
#include <pricer/params/CVol.h>
#include <pricer/math/Cinterp2d.h>
#include <pricer/math/CGrid.h>

namespace Pricer
{
  class CSigmaDupire : public CSigmaLoc
  {
  public:
    CSigmaDupire(const ptr<CVol>& p_vol);
    virtual ~CSigmaDupire() {};
    void Init(const ptr<CGrid>& p_grid);
    double virtual operator()(double p_t, double p_s) const;

  protected:
    ptr<CVol> m_vol;
    ptr<CBase_interp2d> m_inter;

    virtual double evalLocalVol(double p_tenor, double p_strike) = 0;
  };

  class CSigmaDupireBach : public CSigmaDupire
  {
  public:
    CSigmaDupireBach(const ptr<CVol>& p_vol);
  private:
    double evalLocalVol(double p_tenor, double p_strike);
  };

  class CSigmaDupireBS : public CSigmaDupire
  {
  public:
    CSigmaDupireBS(const ptr<CVol>& p_vol);
  private:
    double evalLocalVol(double p_tenor, double p_strike);
  };
}

#endif