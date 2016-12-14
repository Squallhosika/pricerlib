#ifndef pricer_cgrid
#define pricer_cgrid

// TODO it is ok to add this file here
#include <pricer/stdafx.h>

namespace Pricer
{
  class CGrid
  {
    typedef std::vector<double> vectd;
  public:
    CGrid(const vectd& p_timeSteps, const vectd& p_oSpaceSteps);
    vectd TimeStepsDiff() const;
    vectd TimeSteps() const;
    vectd SpaceStepsDiff() const;
    vectd SpaceSteps() const;

  private:
    vectd m_oTimeSteps;
    vectd m_oSpaceSteps;
    
  };

  // Init of the grid
  struct SGener_iter
  {
    SGener_iter(double start, double end, unsigned int nbIter)
    : m_dCurrent(start), m_dstep((end - start) / (nbIter - 1))
    { }

    double operator()()
    {
      double res = m_dCurrent;
      m_dCurrent = res + m_dstep;
      return res;
    }

    double m_dCurrent;
    const double m_dstep;
  };
}

#endif
