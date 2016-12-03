#ifndef pricer_params_ivolderivator
#define pricer_params_ivolderivator

#include "CVol.h"
#include <pricer/UType.h>

namespace Pricer
{
  class IVolDerivator
  {
  public:
    virtual ~IVolDerivator() {};

    double DerivInT(double T, double K);
    double DerivInK(double T, double K);
    // here the 2 is for the square
    double DerivInK2(double T, double K);

    double DerivCallBachInT(double T, double K);
    double DerivCallBachInK(double T, double K);
    double DerivCallBachInK2(double T, double K);

  private:
    // TODO add a param Date member
  };

}


#endif
