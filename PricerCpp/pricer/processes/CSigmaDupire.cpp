#include <pricer/stdafx.h>
#include "CSigmaDupire.h"

namespace Pricer
{
  void Init(const ptr<CVol>& p_vol)
  {

  }

  double CSigmaDupire::operator()(double p_t, double p_s) const
  {
    return m_inter->interp(p_t, p_s);
  }
}