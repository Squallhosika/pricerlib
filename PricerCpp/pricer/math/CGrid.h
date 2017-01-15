#ifndef pricer_cgrid
#define pricer_cgrid

// TODO it is ok to add this file here
#include <pricer/stdafx.h>
#include <pricer/UType.h>
#include <pricer/math/CMatrix.h>
#include <functional>

namespace Pricer
{
  class CGrid
  {
  public:
    CGrid(const vectd& p_timeSteps, const vectd& p_oSpaceSteps);
    CGrid(double p_MinTS, double p_MaxTS, size_t p_nTS,
      double p_MinSpaceStep, double p_MaxSpaceStep, size_t p_nSpaceStep);

    vectd TimeStepsDiff() const;
    vectd TimeSteps() const;
    vectd SpaceStepsDiff() const;
    vectd SpaceSteps() const;

  protected:
    vectd m_oTimeSteps;
    vectd m_oSpaceSteps;

    void sort();
  };

  class CGridWithValue : public CGrid
  {
  public:
    CGridWithValue(const vectd& p_timeSteps, const vectd& p_oSpaceSteps);
    const CMatrix& Values() const;
    void SetValue(size_t p_iTS, size_t p_iSS, double value);
    void Fill(const std::function<double(double, double)>& p_fctTimeSpace);
  private:
    // First by TimeStep then by Tenor
    CMatrix m_oValues;

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
