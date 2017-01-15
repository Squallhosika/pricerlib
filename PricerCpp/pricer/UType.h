#pragma once

#include <memory>
#include <vector>
#include <exception>

template<typename T>
using ptr = std::shared_ptr<T>;

namespace Pricer {

  typedef double Doub;

  // TODO change to ublast after or if more need
  // knowing that the class CTriDagMatrix.h
  // already exist
  typedef std::vector<double>    vectd;
  typedef std::vector<double>    timeSteps;
  typedef std::vector<double>    timeStepsDiff;
  typedef std::vector<double>    path;
  typedef std::vector<ptr<path>> paths;

  class NotImplementedException : public std::exception
  {
  public:
    NotImplementedException()
    : std::exception("Not implemented exception")
    {}
    
  };
}

