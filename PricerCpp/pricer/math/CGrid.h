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
}

#endif
