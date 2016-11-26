#include "stdafx.h"
#include "CGenerator.h"
#include <cstdlib>



namespace Pricer
{
  CNormGen::CNormGen(size_t p_seed, double p_mean, double p_std)
    : m_impl(gen(p_seed), distrib(p_mean, p_std))
  {

  }

  double CNormGen::Next()
  {
    return m_impl();
  }
}