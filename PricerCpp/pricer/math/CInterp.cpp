#include "pricer/stdafx.h"
#include "CInterp.h"



namespace Pricer
{
  int Base_interp::locate(double x)
  {
    int ju, jm, jl;
    if (n < 2 || mm < 2 || mm > n) throw("locate size error");
    bool ascnd = (xx[n - 1] >= xx[0]);
    jl = 0;
    ju = n - 1;
    while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if (x >= xx[jm] == ascnd)
        jl = jm;
      else
        ju = jm;
    }
    cor = abs(jl - jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
  }
}