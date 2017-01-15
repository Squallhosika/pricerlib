#pragma once

#include <pricer/UType.h>

namespace Pricer
{
  class CMatrix : public std::vector<vectd>
  {
  public:
    CMatrix();
    CMatrix(size_t p_iNRow, size_t p_iNCol);
  };
}