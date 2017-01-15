#pragma once

#include <pricer/UType.h>
#include <pricer/external/nr3.h>

namespace Pricer
{
  struct Bracketmethod {
    Doub ax, bx, cx, fa, fb, fc;
    template <class T>
    void bracket(const Doub a, const Doub b, T &func)
    {
      const Doub GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
      ax = a; bx = b;
      Doub fu;
      fa = func(ax);
      fb = func(bx);
      if (fb > fa) {

        SWAP(ax, bx);
        SWAP(fb, fa);
      }
      cx = bx + GOLD*(bx - ax);
      fc = func(cx);
      while (fb > fc) {
        Doub r = (bx - ax)*(fb - fc);
        Doub q = (bx - cx)*(fb - fa);
        Doub u = bx - ((bx - cx)*q - (bx - ax)*r) /
          (2.0*SIGN(MAX(abs(q - r), TINY), q - r));
        Doub ulim = bx + GLIMIT*(cx - bx);

        if ((bx - u)*(u - cx) > 0.0) {
          fu = func(u);
          if (fu < fc) {
            ax = bx;
            bx = u;
            fa = fb;
            fb = fu;
            return;
          }
          else if (fu > fb) {
            cx = u;
            fc = fu;
            return;
          }
          u = cx + GOLD*(cx - bx);
          fu = func(u);
        }
        else if ((cx - u)*(u - ulim) > 0.0) {
          fu = func(u);
          if (fu < fc) {
            shft3(bx, cx, u, u + GOLD*(u - cx));
            shft3(fb, fc, fu, func(u));
          }
        }
        else if ((u - ulim)*(ulim - cx) >= 0.0) {
          u = ulim;
          fu = func(u);
        }
        else {
          u = cx + GOLD*(cx - bx);
          fu = func(u);
        }
        shft3(ax, bx, cx, u);
        shft3(fa, fb, fc, fu);
      }
    }
    inline void shft2(Doub &a, Doub &b, const Doub c)
    {
      a = b;
      b = c;
    }

    inline void shft3(Doub &a, Doub &b, Doub &c, const Doub d)
    {
      a = b;
      b = c;
      c = d;
    }

    inline void mov3(Doub &a, Doub &b, Doub &c, const Doub d, const Doub e,
      const Doub f)
    {
      a = d; b = e; c = f;
    }
  };


  struct Golden : Bracketmethod {
    Doub xmin, fmin;
    const Doub tol;
    Golden(const Doub toll = 3.0e-8) : tol(toll) {}
    template <class T>
    Doub minimize(T &func)
    {
      const Doub R = 0.61803399, C = 1.0 - R;
      Doub x1, x2;
      Doub x0 = ax;
      Doub x3 = cx;
      if (abs(cx - bx) > abs(bx - ax)) {
        x1 = bx;
        x2 = bx + C*(cx - bx);
      }
      else {
        x2 = bx;
        x1 = bx - C*(bx - ax);
      }
      Doub f1 = func(x1);
      Doub f2 = func(x2);
      while (abs(x3 - x0) > tol*(abs(x1) + abs(x2))) {
        if (f2 < f1) {
          shft3(x0, x1, x2, R*x2 + C*x3);
          shft2(f1, f2, func(x2));
        }
        else {
          shft3(x3, x2, x1, R*x1 + C*x0);
          shft2(f2, f1, func(x1));
        }
      }
      if (f1 < f2) {
        xmin = x1;
        fmin = f1;
      }
      else {
        xmin = x2;
        fmin = f2;
      }
      return xmin;
    }
  };

}