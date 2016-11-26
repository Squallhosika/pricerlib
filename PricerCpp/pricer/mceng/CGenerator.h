#ifndef pricer_cgenerator
#define pricer_cgenerator

// To desable security on visual studio
// add -D_SCL_SECURE_NO_WARNINGS in the command
// line
#include <boost/random.hpp>

namespace Pricer
{
  class CNormGen
  {
  public:
    typedef boost::mt19937 gen;
    typedef boost::normal_distribution<> distrib;
    
    CNormGen(size_t p_seed = 0, double p_mean = 0.0, double p_std = 1.0);
    double Next();
  private:
    boost::variate_generator<gen, distrib> m_impl;
  };
}


#endif