#ifndef pricer_utype
#define pricer_utype

#include <memory>

template<typename T>
using ptr = std::shared_ptr<T>;

namespace Pricer {

  // TODO change to ublast after or if more need
  // knowing that the class CTriDagMatrix.h
  // already exist
  typedef std::vector<std::vector<double>> DMatrix;

}

#endif
