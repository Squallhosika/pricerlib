#ifndef pricer_uformula
#define pricer_uformula

#include <algorithm>
#include <boost/math/distributions.hpp>

namespace Pricer
{
  double CallBlack(double p_s, double p_t, double p_k, double p_sig, double p_mat)
  {
    boost::math::normal_distribution<> l_normDistrib;
    double l_d1 = (std::log(p_s / p_k) + 0.5 * p_sig * p_sig * (p_mat - p_t)) / (p_sig * sqrt(p_mat - p_t));
    double l_d2 = (std::log(p_s / p_k) - 0.5 * p_sig * p_sig * (p_mat - p_t)) / (p_sig * sqrt(p_mat - p_t));
    return p_s * boost::math::cdf(l_normDistrib, l_d1) - p_k * boost::math::cdf(l_normDistrib, l_d2);
  }

  double CallDigitBlack(double p_s, double p_t, double p_k, double p_sig, double p_mat)
  {
    boost::math::normal_distribution<> l_normDistrib;
    double l_d2 = (std::log(p_s / p_k) - 0.5 * p_sig * p_sig * (p_mat - p_t)) / (p_sig * sqrt(p_mat - p_t));
    return boost::math::cdf(l_normDistrib, l_d2);
  }

  double CallBlackShifted(double p_s, double p_shift, double p_t, double p_k, double p_sig, double p_mat)
  {
    return CallBlack(p_s, p_t, p_k, p_sig, p_mat);
  }

  double CallDigitBlackShifted(double p_s, double p_shift, double p_t, double p_k, double p_sig, double p_mat)
  {
    return CallDigitBlack(p_s + p_shift, p_t, p_k + p_shift, p_sig, p_mat);
  }

  double CallBachelier(double p_s, double p_t, double p_k, double p_sig, double p_mat)
  {
    boost::math::normal_distribution<> l_normDistrib;
    double l_d1 = (p_s - p_k) / (p_sig * sqrt(p_mat - p_t));
    return (p_s - p_k)  * boost::math::cdf(l_normDistrib, l_d1) 
      + (p_sig * sqrt(p_mat - p_t)) * boost::math::pdf(l_normDistrib, l_d1);
  }


}


#endif
