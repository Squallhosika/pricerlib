#ifndef pricer_uformula
#define pricer_uformula

#include <algorithm>
#include <boost/math/distributions.hpp>

namespace Pricer
{
  double CallBlack(double p_s, double p_t, double p_k, double p_sig, double p_mat);
  double CallDigitBlack(double p_s, double p_t, double p_k, double p_sig, double p_mat);
  double CallBlackShifted(double p_s, double p_shift, double p_t, double p_k, double p_sig, double p_mat);
  double CallDigitBlackShifted(double p_s, double p_shift, double p_t, double p_k, double p_sig, double p_mat);
  double CallBachelier(double p_s, double p_t, double p_k, double p_sig, double p_mat);
}

#endif
