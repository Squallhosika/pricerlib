#include <pricer/stdafx.h>
#include "CPayoff.h"

#include <numeric>
#include <algorithm>

namespace Pricer
{
  CPayoff::~CPayoff()
  {
  }

  CPayoffEUropean::~CPayoffEUropean()
  {
  }

  double CPayoffId::operator()(double p_s) const
  {
    return p_s;
  }

  CPayoffConst::CPayoffConst(double p_const)
    : m_const(p_const)
  {
  }

  double CPayoffConst::operator()(double) const
  {
    return m_const;
  }

  CPayoffCall::CPayoffCall(double p_k)
    : m_k(p_k)
  {
  }

  double CPayoffCall::operator()(double p_s) const
  {
    return std::max(0.0, p_s - m_k);
  }

  CPayoffDigit::CPayoffDigit(double p_k)
    : m_k(p_k)
  {
  }

  double CPayoffDigit::operator()(double p_s) const
  {
    return p_s > m_k ? 1.0 : 0.0;
  }

  CPayoffPathDepPartial::CPayoffPathDepPartial(const timeSteps& p_oNeededTS)
    : m_spNeededTS(p_oNeededTS)
  {}

  CPayoffPathDepPartial::~CPayoffPathDepPartial()
  {
  }

  double CPayoffPathDepPartial::operator()(const std::vector<double>& p_s) const
  {
    if (m_spNeededTS.size() == p_s.size())
      std::exception("Not the good number of points gave.");

    return this->Value(p_s);
  }

  const timeSteps& CPayoffPathDepPartial::NeededTS() const
  {
    return m_spNeededTS;
  }

  CPayoffVolBach::CPayoffVolBach(const timeSteps& p_oNeededTS)
    : CPayoffPathDepPartial(p_oNeededTS)
  {}

  // externile this calcul because the two function are almost
  // similar
  double CPayoffVolBach::Value(const std::vector<double>& p_s) const
  {
    std::vector<double> l_oRnd(p_s.size());
    std::adjacent_difference(p_s.begin(), p_s.end(), l_oRnd.begin());
    l_oRnd.erase(l_oRnd.begin());
    std::transform(l_oRnd.begin(), l_oRnd.end(), l_oRnd.begin(),
      [](double x) {return x * x; });

    double l_dRes = 0.0;
    for (double value : l_oRnd)
      l_dRes += value;

    return std::sqrt(l_dRes / (m_spNeededTS.back() - m_spNeededTS.front()));
  }

  CPayoffVarBach::CPayoffVarBach(const timeSteps& p_oNeededTS)
    : CPayoffPathDepPartial(p_oNeededTS)
  {}

  double CPayoffVarBach::Value(const std::vector<double>& p_s) const
  {
    std::vector<double> l_oRnd(p_s.size());
    std::adjacent_difference(p_s.begin(), p_s.end(), l_oRnd.begin());
    l_oRnd.erase(l_oRnd.begin());
    std::transform(l_oRnd.begin(), l_oRnd.end(), l_oRnd.begin(),
      [](double x) {return x * x; });

    double l_dRes = 0.0;
    for (double value : l_oRnd)
      l_dRes += value;

    return l_dRes / (m_spNeededTS.back() - m_spNeededTS.front());
  }

}