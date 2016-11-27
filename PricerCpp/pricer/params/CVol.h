#ifndef pricer_params_cvol
#define pricer_params_cvol


namespace Pricer
{
  class CVol
  {
  public:
    virtual ~CVol() {};
    double virtual GetPoint(double tenor, double strike) = 0;

  private:
    // TODO add a param Date member
  };

}


#endif
