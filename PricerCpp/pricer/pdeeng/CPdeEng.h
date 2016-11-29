#ifndef pricer_cpdeeng
#define pricer_cpdeeng

#include <vector>
#include <algorithm>

#include "pricer/UType.h"
#include "pricer/processes/CProcessVolLoc.h"
#include "pricer/math/CTriDagMatrix.h"
#include "pricer/instrs/CPayoff.h"

namespace Pricer
{
  class CPdeEng
  {
    // implicit for theta  = 0.0
    // explicit for thetat = 1.0
  public:
    CPdeEng(
      const ptr<CProcessVolLoc>& p_spProcess,
      const ptr<CPayoffEUropean>& p_payoff,
      const std::vector<double>& p_timeSteps,
      double p_underMin, double p_underMax,
      double p_theta = 0.0);

    std::vector<double> UnderValues() const;
    size_t SizeInT() const;
    size_t SizeInH() const;

    void OverloadBoundary(double p_mult, double p_s0, double p_mat, bool p_pos);
    void InitH(size_t p_nb);
    std::vector<double> Compute() const;

  private:
    typedef std::vector<double>::const_iterator vciter;
    CPdeEng();

    ptr<CProcessVolLoc>  m_spProcess;
    ptr<CPayoffEUropean> m_payoff;
    std::vector<double>  m_reversTS;
    double               m_underMin;
    double               m_underMax;
    double               m_h;
    // TODO change h constant to variant:
    //std::vector<double> m_underStep;
    std::vector<double>  m_underValue;
    double               m_theta;

    
    std::vector<double>   valuesAtMat() const;
    std::vector<double>   vectorR(vciter p_iterTS, double p_time, 
      const std::vector<double>& p_v) const;
    CTriDagMatrix<double> matrixTriDag(vciter p_iterTS, double p_time, 
      const std::vector<double>& p_v) const;
  };

}
#endif
