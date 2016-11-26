#ifndef pricer_cpathgen
#define pricer_cpathgen


#include "IProcess.h"
#include "UType.h"
#include <vector>

namespace Pricer
{
  template<typename GEN, typename STAT>
  class CPathGen
  {
  public:
    typedef typename STAT::timeSteps timeSteps;
    typedef typename STAT::path path;

    CPathGen(const ptr<IProcess>& p_spProcess, const GEN& p_gen, const ptr<timeSteps>& p_spTs, size_t p_nbSimu);
    //~CPathGen();

    void GenSequence();
    void FillStat(double p_s0, STAT& p_stat);

  private:
    CPathGen();
    ptr<IProcess> m_spProcess;
    GEN           m_gen;
    size_t        m_nbSimu;
    // TODO coherence between CPathGen TimeStep and STAT TimeStep
    ptr<timeSteps> m_spTimeSteps;
    std::vector<ptr<path>> m_paths;
    bool m_isGen;
  };

  template<typename GEN, typename STAT>
  CPathGen<GEN, STAT>::CPathGen(const ptr<IProcess>& p_spProcess, const GEN& p_gen, const ptr<timeSteps>& p_spTs, size_t p_nbSimu)
    : m_spProcess(p_spProcess), m_gen(p_gen), m_spTimeSteps(p_spTs), m_nbSimu(p_nbSimu), m_isGen(false)
  {
  }

  template<typename GEN, typename STAT>
  void CPathGen<GEN, STAT>::GenSequence()
  {
    if (m_isGen)
      return;

    m_paths.reserve(m_nbSimu);

    for (size_t l_i = 0; l_i < m_nbSimu; l_i++)
    {
      ptr<path> l_spPath = std::make_shared<path>();
      l_spPath->reserve(m_spTimeSteps->size());

      for (double l_step : *m_spTimeSteps)
        l_spPath->push_back(m_gen.Next());

      m_paths.push_back(l_spPath);
    }

    m_isGen = true;
  }

  template<typename GEN, typename STAT>
  void CPathGen<GEN, STAT>::FillStat(double p_s0, STAT& p_stat)
  {
    path::const_iterator l_iterPathGen;
    double l_s;

    for (ptr<path> l_spPathGen : m_paths)
    {
      l_s = p_s0;
      ptr<path> l_spPath = std::make_shared<path>();
      l_spPath->reserve(m_spTimeSteps->size());
      
      l_iterPathGen = (*l_spPathGen).begin();

      for (double l_step : *m_spTimeSteps)
      {
        l_s = m_spProcess->evol(l_s, l_step, *l_iterPathGen);
        l_spPath->push_back(l_s);
        l_iterPathGen++;
      }

      p_stat.Add(l_spPath);
    }
  }

}


#endif
