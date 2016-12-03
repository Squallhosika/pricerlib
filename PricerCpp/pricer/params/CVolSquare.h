#ifndef pricer_params_cvolsquare
#define pricer_params_cvolsquare

#include <pricer/UType.h>
#include <pricer/math/Cinterp2d.h>
#include <vector>
#include "CVol.h"

namespace Pricer
{
  // TODO make the interpolatyion generic
  class CVolSquare : public CVol
  {
  public:
    CVolSquare(const std::vector<double>& p_oTenors, 
      const std::vector<double> p_oStrikes, 
      const DMatrix& p_oPoints);
    virtual ~CVolSquare() {};
    double virtual GetPoint(double tenor, double strike) const;

    

  private:
    std::vector<double> m_oTenors;
    std::vector<double> m_oStrikes;
    DMatrix m_oPoints;
    CVar_spline_interp2d m_oInter;
    
  };

}


#endif
