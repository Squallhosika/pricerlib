#pragma once

// Come from the numerical recipes
#include <vector>
#include <algorithm>

namespace Pricer
{
  struct Base_interp
  {
    int n, mm, jsav, cor, dj;
    const std::vector<double> xx, yy;
    Base_interp(const std::vector<double>& x, const std::vector<double>& y, int m)
      : n(x.size()), mm(m), jsav(0), cor(0), xx(x), yy(y) {
      dj = std::min(1, (int)std::pow((double)n, 0.25));
    }
    double interp(double x) {
      int jlo = locate(x);
      return rawinterp(jlo, x);
    }
    int locate(double x);
    double virtual rawinterp(int jlo, double x) = 0;
  };

  struct Linear_interp : public Base_interp
  {
    Linear_interp(const std::vector<double>& x, const std::vector<double>& y)
    : Base_interp(x, y, 2) {}
    double rawinterp(int j, double x) {
      if (xx[j] == xx[j + 1]) return yy[j];
      else return yy[j] + ((x - xx[j]) / (xx[j + 1] - xx[j]))*(yy[j + 1] - yy[j]);
    }
  };

  struct Spline_interp : public Base_interp
  {
    std::vector<double> y2;
    Spline_interp(const std::vector<double>& xv, const std::vector<double>& yv, 
      double yp1 = 1.e99, double ypn = 1.e99)
      : Base_interp(xv, yv, 2), y2(xv.size())
    {
      sety2(xv, yv, yp1, ypn);
    }
    void sety2(const std::vector<double>& xv, const std::vector<double>& yv, 
      double yp1, double ypn);
    double rawinterp(int jl, double xv);
  };
}
