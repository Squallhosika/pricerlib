#pragma once

namespace Pricer
{
  enum class EPricingType
  {
    E_MC,
    E_PDE,
    E_MAX
  };

  struct SPricerOption
  {
    EPricingType PricingType = EPricingType::E_MAX;

    // MC Part
    bool MCisMilstein = false;
    size_t MCnbSimu = 10000;
    size_t MCnbTimeStep = 100;

    // PDE Part
    double PDEtheta = 0.0;
    size_t PDEsizeInT = 100;
    size_t PDEsizeInH = 500;

    // Dupire rep
    unsigned int DupNInTime = 101;
    unsigned int DupNInSpace = 1001;
  };
}

