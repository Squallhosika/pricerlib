#ifndef pricer_params_cvol
#define pricer_params_cvol

#include <pricer/UType.h>
#include <pricer/math/Cinterp2d.h>
#include <vector>


namespace Pricer
{
  class CVolSquare
  {
  public:
    CVolSquare(const std::vector<double>& p_oTenors, 
      const std::vector<double> p_oStrikes, 
      const DMatrix& p_oPoints);
    virtual ~CVolSquare() {};
    double virtual GetPoint(double tenor, double strike);

  private:
    std::vector<double> m_oTenors;
    std::vector<double> m_oStrikes;
    DMatrix m_oPoints;
    CVar_spline_interp2d m_oInter;
  };

}


#endif
