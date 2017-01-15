#include <pricer/stdafx.h>
#include "UVolInverser.h"

#include <functional>
#include "CMinimazer.h"
#include <pricer/tools/UFormulas.h>

namespace Pricer
{
  namespace UVolInverser
  {

    double InversCallBach(double p_s, double p_t, double p_k, double p_price, double p_mat,
      double p_volMin, double p_volMax)
    {
       std::function<double(double)> l_function =
        [p_s, p_t, p_k, p_price, p_mat](double x) -> double
      {
        double l_diff = UFormulas::CallBachelier(p_s, p_t, p_k, x, p_mat) - p_price;
        return l_diff * l_diff;
      };

      Golden golden;
      golden.bracket(p_volMin, p_volMax, l_function);
      return golden.minimize(l_function);
    }
  }
}