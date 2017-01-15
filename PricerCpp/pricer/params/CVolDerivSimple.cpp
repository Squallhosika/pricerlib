#include "pricer/stdafx.h"
#include "CVolDerivSimple.h"

#include <pricer/tools/UFormulas.h>

namespace Pricer
{

  CVolDerivSimple::CVolDerivSimple(CVol* p_vol)
    : m_vol(p_vol)
  {
  }


  CVolDerivSimple::~CVolDerivSimple()
  {
  }

  double CVolDerivSimple::DerivInT(double T, double K) const
  {
    return 0.5 * (m_vol->GetPoint(T + m_tenorEpsi, K)
      - m_vol->GetPoint(T - m_tenorEpsi, K)) / m_tenorEpsi;
  }

  double CVolDerivSimple::DerivInK(double T, double K) const
  {
    return 0.5 * (m_vol->GetPoint(T, K + m_strikeEpsi)
      - m_vol->GetPoint(T, K - m_strikeEpsi)) / m_strikeEpsi;
  }

  double  CVolDerivSimple::DerivInK2(double T, double K) const
  {
    return (m_vol->GetPoint(T, K + m_strikeEpsi) + m_vol->GetPoint(T, K - m_strikeEpsi)
      - 2.0 * m_vol->GetPoint(T, K)) / (m_strikeEpsi * m_strikeEpsi);
  }

  double CVolDerivSimple::DerivCallBachInT(double T, double K) const
  {
    // TODO no symetric to resolve issue about negative T
    // it is the good way to fo it ??
    return (UFormulas::CallBachelier(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T + m_tenorEpsi, K), T + m_tenorEpsi)
      - UFormulas::CallBachelier(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T, K), T))
      / m_tenorEpsi;
  }

  double CVolDerivSimple::DerivCallBachInK(double T, double K) const
  {
    return 0.5 * (UFormulas::CallBachelier(m_vol->Spot(), 0.0, K + m_strikeEpsi, m_vol->GetPoint(T, K + m_strikeEpsi), T)
      - UFormulas::CallBachelier(m_vol->Spot(), 0.0, K - m_strikeEpsi, m_vol->GetPoint(T, K - m_strikeEpsi), T))
      / m_strikeEpsi;
  }

  double CVolDerivSimple::DerivCallBachInK2(double T, double K) const
  {
    // DEBUGTODO the below line
    double l_vol = m_vol->GetPoint(T, K + m_strikeEpsi);
    return (UFormulas::CallBachelier(m_vol->Spot(), 0.0, K + m_strikeEpsi, m_vol->GetPoint(T, K + m_strikeEpsi), T)
      + UFormulas::CallBachelier(m_vol->Spot(), 0.0, K - m_strikeEpsi, m_vol->GetPoint(T, K - m_strikeEpsi), T)
      - 2.0 * UFormulas::CallBachelier(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T, K), T))
      / (m_strikeEpsi * m_strikeEpsi);
  }

  // WARNTODO Black instead of BS in the CallBlack Formula !!
  double CVolDerivSimple::DerivCallBSInT(double T, double K) const
  {
    // TODO no symetric to resolve issue about negative T
    // it is the good way to fo it ??
    return (UFormulas::CallBlack(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T + m_tenorEpsi, K), T + m_tenorEpsi)
      - UFormulas::CallBlack(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T, K), T))
      / m_tenorEpsi;
  }

  // WARNTODO Black instead of BS in the CallBlack Formula !!
  double CVolDerivSimple::DerivCallBSInK(double T, double K) const
  {
    return 0.5 * (UFormulas::CallBlack(m_vol->Spot(), 0.0, K + m_strikeEpsi, m_vol->GetPoint(T, K + m_strikeEpsi), T)
      - UFormulas::CallBlack(m_vol->Spot(), 0.0, K - m_strikeEpsi, m_vol->GetPoint(T, K - m_strikeEpsi), T))
      / m_strikeEpsi;
  }

  // WARNTODO Black instead of BS in the CallBlack Formula !!
  double CVolDerivSimple::DerivCallBSInK2(double T, double K) const
  {
    // DEBUGTODO the below line
    double l_vol = m_vol->GetPoint(T, K + m_strikeEpsi);
    return (UFormulas::CallBlack(m_vol->Spot(), 0.0, K + m_strikeEpsi, m_vol->GetPoint(T, K + m_strikeEpsi), T)
      + UFormulas::CallBlack(m_vol->Spot(), 0.0, K - m_strikeEpsi, m_vol->GetPoint(T, K - m_strikeEpsi), T)
      - 2.0 * UFormulas::CallBlack(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T, K), T))
      / (m_strikeEpsi * m_strikeEpsi);
  }

}