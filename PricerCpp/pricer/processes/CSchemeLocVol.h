#ifndef pricer_cschemelocvol
#define pricer_cschemelocvol

#include "pricer/UType.h"
#include "CSigmaLoc.h"

namespace Pricer
{
  class CSchemeLocalVol
  {
  public:
    CSchemeLocalVol(const ptr<CSigmaLoc>& p_sigma);
    double virtual evol(double p_t, double p_s, double p_step, double p_w) const = 0;

  protected:
    ptr<CSigmaLoc> m_sigma;
  };

  class CSchemeLocVolBach : public CSchemeLocalVol
  {
  public:
    CSchemeLocVolBach(const ptr<CSigmaLoc>& p_sigma);
    double virtual evol(double p_t, double p_s, double p_step, double p_w) const = 0;

    static ptr<CSchemeLocVolBach> Factory(const ptr<CSigmaLoc>& p_sigma, 
      bool p_isMillstein = false);
  };

  class CSchemeLocVolBachEuler : public CSchemeLocVolBach
  {
  public:
    CSchemeLocVolBachEuler(const ptr<CSigmaLoc>& p_sigma)
      : CSchemeLocVolBach(p_sigma) {};

    double virtual evol(double p_t, double p_s, double p_step, double p_w) const;
  };

  class CSchemeLocVolBachMilstein : public CSchemeLocVolBach
  {
  public:
    CSchemeLocVolBachMilstein(const ptr<CSigmaLoc>& p_sigma)
      : CSchemeLocVolBach(p_sigma) {};

    double virtual evol(double p_t, double p_s, double p_step, double p_w) const;
  };

  class CSchemeLocVolBS : public CSchemeLocalVol
  {
  public:
    CSchemeLocVolBS(const ptr<CSigmaLoc>& p_sigma);
    double virtual evol(double p_t, double p_s, double p_step, double p_w) const = 0;

    static ptr<CSchemeLocVolBS> Factory(const ptr<CSigmaLoc>& p_sigma,
      bool p_isMillstein = false);
  };

  class CSchemeLocVolBSEuler : public CSchemeLocVolBS
  {
  public:
    CSchemeLocVolBSEuler(const ptr<CSigmaLoc>& p_sigma)
      : CSchemeLocVolBS(p_sigma) {};

    double virtual evol(double p_t, double p_s, double p_step, double p_w) const;
  };

  class CSchemeLocVolBSMilstein : public CSchemeLocVolBS
  {
  public:
    CSchemeLocVolBSMilstein(const ptr<CSigmaLoc>& p_sigma)
      : CSchemeLocVolBS(p_sigma) {};

    double virtual evol(double p_t, double p_s, double p_step, double p_w) const;
  };
}



// Old Scheme maybe to complex fo what we need

//
//template<typename SIGMA>
//class CSchemeLocVol
//{
//public:
//  // We do not overload the destructor because
//  // the scheme is containt by the object
//  // on which he point so we do not want to 
//  // destroy it when the scheme is destroy
//  CSchemeLocVol(CProcessVolLoc<SIGMA> *p_process);
//  double virtual evol(double p_s, double p_step, double p_w) const = 0;
//
//  static ptr<CSchemeLocVol> Factory(CProcessVolLoc<SIGMA> *p_process, bool p_isMillstein = false);
//protected:
//  CProcessVolLoc<SIGMA> *m_process;
//};
//
//template<typename SIGMA>
//CSchemeLocVol<SIGMA>::CSchemeLocVol(CProcessVolLoc<SIGMA> *p_process)
//:m_process(p_process)
//{}
//
//template<typename SIGMA>
//static ptr<CSchemeLocVol> CSchemeLocVol<SIGMA>::Factory(CProcessVolLoc<SIGMA> *p_process, bool p_isMillstein)
//{
//  if (p_isMillstein)
//    return std::make_shared<CSchemeLocVolMilstein>(p_process);
//
//  return std::make_shared<CSchemeLocVolEuler>(p_process);
//}
//
//template<typename SIGMA>
//class CSchemeLocVolEuler : public CSchemeLocVol
//{
//public:
//  CSchemeLocVolEuler(CProcessVolLoc<SIGMA> *p_process)
//    : CSchemeLocVol(p_process) {};
//
//  double evol(double p_s, double p_step, double p_w) const = 0;
//};
//
//template<typename SIGMA>
//double CSchemeLocVolEuler<SIGMA>::evol(double p_s, double p_step, double p_w) const
//{
//  return p_s + m_process->Sigma(p_s) * std::sqrt(p_step) * p_w;
//}
//
//template<typename SIGMA>
//class CSchemeLocVolMilstein : public CSchemeLocVol
//{
//public:
//  CSchemeLocVolEuler(CProcessVolLoc<SIGMA> *p_process)
//    : CSchemeLocVol(p_process) {};
//
//  double virtual evol(double p_s, double p_step, double p_w) const = 0;
//};
//
//template<typename SIGMA>
//double CSchemeLocVolEuler<SIGMA>::evolMill(double p_s, double p_step, double p_w) const
//{
//  return p_s + m_process->Sigma(p_s)* std::sqrt(p_step) * p_w
//    + 0.5 * m_process->Sigma(p_s)* m_process->Sigma().Deriv(p_s) * p_step * (p_w  * p_w - p_step);
//}

#endif