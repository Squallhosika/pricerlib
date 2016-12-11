#ifndef pricer_cpathgen
#define pricer_cpathgen

#include <pricer/stdafx.h>
#include "CEurStat.h"
#include <pricer/processes/IProcess.h>
#include <pricer/UType.h>
#include <vector>

namespace Pricer
{
  template<typename GEN>
  class CPathGen
  {
  public:
    CPathGen(const ptr<IProcess>& p_spProcess, const GEN& p_gen, 
      const ptr<timeStepsDiff>& p_spTs, size_t p_nbSimu);
    //~CPathGen();

    void GenSequence();
    void FillStat(double p_s0, CStat& p_stat);

  private:
    CPathGen();
    ptr<IProcess> m_spProcess;
    GEN           m_gen;
    size_t        m_nbSimu;
    // TODO coherence between CPathGen TimeStep and STAT TimeStep
    ptr<timeSteps> m_spTimeSteps;
    ptr<timeStepsDiff> m_spTimeStepsDiff;
    std::vector<ptr<path>> m_paths;
    bool m_isGen;
  };

  template<typename GEN>
  CPathGen<GEN>::CPathGen(const ptr<IProcess>& p_spProcess, const GEN& p_gen, 
    const ptr<timeStepsDiff>& p_spTs, size_t p_nbSimu)
    : m_spProcess(p_spProcess), m_gen(p_gen), m_spTimeStepsDiff(p_spTs),
    m_nbSimu(p_nbSimu), m_isGen(false), m_spTimeSteps(new timeSteps())
  {
    m_spTimeSteps->resize(m_spTimeStepsDiff->size() + 1);
    (*m_spTimeSteps)[0] = 0.0;
    std::partial_sum(m_spTimeStepsDiff->begin(), m_spTimeStepsDiff->end(), ++m_spTimeSteps->begin());
  }

  template<typename GEN>
  void CPathGen<GEN>::GenSequence()
  {
    if (m_isGen)
      return;

    m_paths.reserve(m_nbSimu);

    for (size_t l_i = 0; l_i < m_nbSimu; l_i++)
    {
      ptr<path> l_spPath = std::make_shared<path>();
      l_spPath->reserve(m_spTimeSteps->size());

      for (double l_step : *m_spTimeStepsDiff)
        l_spPath->push_back(m_gen.Next());

      m_paths.push_back(l_spPath);
    }

    m_isGen = true;
  }

  template<typename GEN>
  void CPathGen<GEN>::FillStat(double p_s0, CStat& p_stat)
  {
    path::const_iterator l_iterPathGen;
    double l_s;
    double l_t;

    for (ptr<path> l_spPathGen : m_paths)
    {
      l_s = p_s0;
      l_t = 0.0;

      ptr<path> l_spPath = std::make_shared<path>();
      l_spPath->reserve(m_spTimeSteps->size());
      l_spPath->push_back(l_s);

      l_iterPathGen = (*l_spPathGen).begin();

      for (double l_step : *m_spTimeStepsDiff)
      {
        l_s = m_spProcess->evol(l_t, l_s, l_step, *l_iterPathGen);
        l_t += l_step;
        l_spPath->push_back(l_s);
        l_iterPathGen++;
      }

      p_stat.Add(l_spPath);
    }
  }

}


#endif
