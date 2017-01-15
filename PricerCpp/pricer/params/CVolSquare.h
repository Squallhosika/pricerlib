#ifndef pricer_params_cvolsquare
#define pricer_params_cvolsquare

#include <pricer/UType.h>
#include <pricer/math/Cinterp2d.h>
#include <pricer/math/CGrid.h>
#include <vector>
#include "CVol.h"
#include "CVolDerivSimple.h"

namespace Pricer
{
  // TODO make the interpolatyion generic
  class CVolSquare : public CVol
  {
  public:
    CVolSquare(double p_spot, const std::vector<double>& p_oTenors, 
      const std::vector<double> p_oStrikes, const CMatrix& p_oPoints);
    CVolSquare(double p_spot, const CGridWithValue& p_spGrid);
    virtual ~CVolSquare() {};
    double virtual GetPoint(double tenor, double strike) const;
    // Give a grid for the point you need to match
    ptr<CGridWithValue> GridToMath();

  private:
    std::vector<double> m_oTenors;
    std::vector<double> m_oStrikes;
    CMatrix m_oPoints;
    CSplineInY_interp2d m_oInter;
  };

}


#endif
