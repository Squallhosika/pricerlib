#ifndef pricer_iprocess
#define pricer_iprocess

namespace Pricer
{
  class IProcess
  {
  public:
    virtual ~IProcess() {};
    double virtual evol(double p_t, double p_s, double p_step, double p_w) const = 0;
  };

}

#endif
