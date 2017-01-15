#include <pricer/stdafx.h>
#include "CMatrix.h"

namespace Pricer
{
  CMatrix::CMatrix()
  {}

  CMatrix::CMatrix(size_t p_iNRow, size_t p_iNCol)
    : std::vector<vectd>(p_iNRow, vectd(p_iNCol))
  {
  }
}