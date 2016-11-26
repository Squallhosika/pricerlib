#ifndef pricer_cschemelocvol
#define pricer_cschemelocvol

#include "UType.h"

namespace Pricer
{
  template<typename SIGMA>
  class CSchemeLocVol
  {
  public:
    // We do not overload the destructor because
    // the scheme is containt by the object
    // on which he point so we do not want to 
    // destroy it when the scheme is destroy
    CSchemeLocVol(SIGMA p_sigma);
    double virtual evol(double p_s, double p_step, double p_w) const = 0;

    static ptr<CSchemeLocVol> Factory(SIGMA p_sigma, bool p_isMillstein = false);
  protected:
    SIGMA m_sigma;
  };

  template<typename SIGMA>
  CSchemeLocVol<SIGMA>::CSchemeLocVol(SIGMA p_sigma)
    :m_sigma(p_sigma)
  {}

  template<typename SIGMA>
  ptr<CSchemeLocVol<SIGMA>> CSchemeLocVol<SIGMA>::Factory(SIGMA p_sigma, bool p_isMillstein)
  {
    if (p_isMillstein)
      return std::static_pointer_cast<CSchemeLocVol<SIGMA>>(std::make_shared<CSchemeLocVolMilstein<SIGMA>>(p_sigma));

    return std::static_pointer_cast<CSchemeLocVol<SIGMA>>(std::make_shared<CSchemeLocVolEuler<SIGMA>>(p_sigma));
  }

  template<typename SIGMA>
  class CSchemeLocVolEuler : public CSchemeLocVol<SIGMA>
  {
  public:
    CSchemeLocVolEuler(SIGMA p_sigma)
      : CSchemeLocVol(p_sigma) {};

    double evol(double p_s, double p_step, double p_w) const;
  };

  template<typename SIGMA>
  double CSchemeLocVolEuler<SIGMA>::evol(double p_s, double p_step, double p_w) const
  {
    return p_s + m_sigma(p_s) * std::sqrt(p_step) * p_w;
  }

  template<typename SIGMA>
  class CSchemeLocVolMilstein : public CSchemeLocVol<SIGMA>
  {
  public:
    CSchemeLocVolMilstein(SIGMA p_sigma)
      : CSchemeLocVol(p_sigma) {};

    double evol(double p_s, double p_step, double p_w) const;
  };

  template<typename SIGMA>
  double CSchemeLocVolMilstein<SIGMA>::evol(double p_s, double p_step, double p_w) const
  {
    double l_sigma = m_sigma(p_s);
    return p_s + l_sigma * std::sqrt(p_step) * p_w
      + 0.5 * l_sigma * m_sigma.Deriv(p_s) * p_step * (p_w  * p_w - 1.0);
  }


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