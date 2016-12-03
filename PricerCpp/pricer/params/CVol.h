#ifndef pricer_params_cvol
#define pricer_params_cvol

#include <pricer/UType.h>
#include "IVolDerivator.h"

namespace Pricer
{
  // TODO in the actual implementation all the volatilty are
  // volatility in strike for simplicity reason
  // but we should implement after different type of voltility
  // monyness, ... and vol type here bachelier but could ne BS,
  // BS shifted
  class IVolDerivator;
  class CVol
  {
  public:
    virtual ~CVol() {};
    double virtual GetPoint(double tenor, double strike) const = 0;
    double Spot() const { return m_spot; };

    // TODO check the influence of the return type betwwen
    // ptr<IVolDerivator>, IVolDerivator, IVolDerivator&, ...

    ptr<IVolDerivator> Derivator() const { return m_spDerivator; };

  private:
    // TODO simpification not all the volatility should include their
    // spot
    ptr<IVolDerivator> m_spDerivator;
    double m_spot;
  };

}


#endif
