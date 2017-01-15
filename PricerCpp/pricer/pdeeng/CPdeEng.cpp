#include <pricer/stdafx.h>
#include "CPdeEng.h"

#include <numeric>

namespace Pricer
{
  CPdeEng::CPdeEng(
    const ptr<CProcessVolLocBach>& p_spProcess,
    const ptr<CPayoffEUropean>& p_payoff,
    const std::vector<double>& p_timeSteps,
    double p_underMin, double p_underMax,
    double p_theta)
    : m_spProcess(p_spProcess), m_payoff(p_payoff),
    m_reversTS(p_timeSteps), m_underMin(p_underMin),
    m_underMax(p_underMax), m_theta(p_theta)
  {
    std::reverse(m_reversTS.begin(), m_reversTS.end());
  }

  std::vector<double> CPdeEng::UnderValues() const
  {
    return m_underValue;
  }

  size_t CPdeEng::SizeInT() const
  {
    // we only have the delta for Dt so we add + 1
    return m_reversTS.size() + 1;
  }

  size_t CPdeEng::SizeInH() const
  {
    return m_underValue.size();
  }

  void CPdeEng::OverloadBoundary(double p_mult, double p_s0, double p_mat, bool p_pos)
  {
    double l_window;

    // kind of a arbitrary choice for S
    // TODO Modifiy that !!!!
    l_window = p_mult * m_spProcess->StDev(0.0, p_s0, p_mat);

    if (!p_pos)
    {
      m_underMax = p_s0 + 0.5 * l_window;
      m_underMin = p_s0 - 0.5 * l_window;
    }
    else
    {
      m_underMax = p_s0 + l_window;
      m_underMin = 0.0;
    }
  }

  void CPdeEng::InitH(size_t p_nb)
  {
    m_h = (m_underMax - m_underMin) / (p_nb - 1);

    m_underValue.clear();
    m_underValue.reserve(p_nb);

    for (int l_i = p_nb - 1; l_i > -1; l_i--)
      m_underValue.push_back(l_i * m_h + m_underMin);
  }


  std::vector<double> CPdeEng::valuesAtMat() const
  {
    std::vector<double> l_res;
    l_res.reserve(SizeInH());

    for each (double l_h in m_underValue)
      l_res.push_back(m_payoff->operator()(l_h));

    return l_res;
  }

  std::vector<double> CPdeEng::vectorR(vciter p_iterTS, double p_time,
    const std::vector<double>& p_v) const
  {
    std::vector<double> l_res;
    l_res.reserve(SizeInH());

    double l_vKp1;
    double l_vK;
    double l_vKm1;
    double l_sigLoc;
    // theta * sig^2(l_vK) * delta t / h^2
    double l_coef;
    vciter l_iterV;
    vciter l_iterUnder;

    l_iterV = p_v.begin();
    // boundary value
    l_res.push_back(l_vK = *l_iterV);

    l_iterV++;
    l_vKm1 = *l_iterV;

    l_iterUnder = ++(m_underValue.begin());

    for (l_iterV++; l_iterV != p_v.end(); l_iterV++, l_iterUnder++)
    {
      l_vKp1 = l_vK;
      l_vK = l_vKm1;
      l_vKm1 = *l_iterV;
      l_sigLoc = m_spProcess->Sigma()->operator()(p_time, *l_iterUnder);

      l_coef = 0.5 * m_theta * l_sigLoc * l_sigLoc * (*p_iterTS) / (m_h * m_h);
      l_res.push_back(l_coef * (l_vKp1 + l_vKm1) + (1 - 2.0 * l_coef) * l_vK);
    }
    // boundary value
    l_res.push_back(l_vKm1);

    return l_res;
  }

  CTriDagMatrix<double> CPdeEng::matrixTriDag(vciter p_iterTS, double p_time,
    const std::vector<double>& p_v) const
  {
    std::vector<double> l_sup;
    std::vector<double> l_cent;
    std::vector<double> l_inf;
    double              l_sigLoc;
    // theta * sig^2(l_vK) * delta t / h^2
    double               l_coef;
    vciter               l_iterUnder;

    l_sup.reserve(SizeInH() - 1);
    l_cent.reserve(SizeInH());
    l_inf.reserve(SizeInH() - 1);

    // boundary value
    l_sup.push_back(0.0);
    l_cent.push_back(1.0);

    l_iterUnder = ++(m_underValue.begin());

    for (vciter l_iterV = ++(p_v.begin()); l_iterV != --(p_v.end()); l_iterV++, l_iterUnder++)
    {
      l_sigLoc = m_spProcess->Sigma()->Value(p_time, *l_iterUnder);
      l_coef = -0.5 * (1.0 - m_theta) * l_sigLoc * l_sigLoc * (*p_iterTS) / (m_h * m_h);

      l_sup.push_back(l_coef);
      l_cent.push_back(1.0 - 2.0 * l_coef);
      l_inf.push_back(l_coef);
    }
    // boundary value
    l_cent.push_back(1.0);
    l_inf.push_back(0.0);

    return CTriDagMatrix<double>(l_sup, l_cent, l_inf);
  }

  std::vector<double> CPdeEng::Compute() const
  {
    int l_counter = 0;
    std::vector<double> l_res = valuesAtMat();
    double l_time = std::accumulate(m_reversTS.begin(), m_reversTS.end(), 0.0);

    for (vciter l_iterTS = m_reversTS.begin(); l_iterTS != m_reversTS.end(); l_iterTS++)
    {
      //CTriDagMatrix<double> l_matrix = matrixTriDag(l_iterTS, l_res);
      //std::vector<double>   l_vector = vectorR(l_iterTS, l_res);
      //l_res = l_matrix.LUDemcop(l_vector);
      l_res = matrixTriDag(l_iterTS, l_time, l_res)
        .LUDemcop(vectorR(l_iterTS, l_time, l_res));

      l_time -= *l_iterTS;
    }

    return l_res;
  }

}