#ifndef pricer_utype
#define pricer_utype

#include <memory>
#include <vector>

template<typename T>
using ptr = std::shared_ptr<T>;

namespace Pricer {

  typedef double Doub;

  // TODO change to ublast after or if more need
  // knowing that the class CTriDagMatrix.h
  // already exist
  typedef std::vector<std::vector<double>> DMatrix;
  typedef std::vector<double>    timeSteps;
  typedef std::vector<double>    timeStepsDiff;
  typedef std::vector<double>    path;
  typedef std::vector<ptr<path>> paths;
}

#endif
