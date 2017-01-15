#pragma once

namespace Pricer
{
  namespace UVolInverser
  {
    double const VolMin = 0.001;
    double const VolMax = 10.0;

    double InversCallBach(double p_s, double p_t, double p_k, 
      double p_price, double p_mat,
      double p_volMin, double p_volMax);
  }
}