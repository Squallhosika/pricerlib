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

    double virtual DerivInT(double T, double K) const = 0;
    double virtual DerivInK(double T, double K) const = 0;
    // here the 2 is for the square
    double virtual DerivInK2(double T, double K) const = 0;

    double virtual DerivCallBachInT(double T, double K) const = 0;
    double virtual DerivCallBachInK(double T, double K) const = 0;
    double virtual DerivCallBachInK2(double T, double K) const = 0;

    double virtual DerivCallBSInT(double T, double K) const = 0;
    double virtual DerivCallBSInK(double T, double K) const = 0;
    double virtual DerivCallBSInK2(double T, double K) const = 0;

  private:
    // TODO add a param Date member
  };

}


#endif
