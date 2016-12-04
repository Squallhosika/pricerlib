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
    return (CallBachelier(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T + m_tenorEpsi, K), T + m_tenorEpsi)
      - CallBachelier(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T, K), T))
      / m_tenorEpsi;
  }

  double CVolDerivSimple::DerivCallBachInK(double T, double K) const
  {
    return 0.5 * (CallBachelier(m_vol->Spot(), 0.0, K + m_strikeEpsi, m_vol->GetPoint(T, K + m_strikeEpsi), T)
      - CallBachelier(m_vol->Spot(), 0.0, K - m_strikeEpsi, m_vol->GetPoint(T, K - m_strikeEpsi), T))
      / m_strikeEpsi;
  }

  double CVolDerivSimple::DerivCallBachInK2(double T, double K) const
  {
    return (CallBachelier(m_vol->Spot(), 0.0, K + m_strikeEpsi, m_vol->GetPoint(T, K + m_strikeEpsi), T)
      + CallBachelier(m_vol->Spot(), 0.0, K - m_strikeEpsi, m_vol->GetPoint(T, K - m_strikeEpsi), T)
      - 2.0 * CallBachelier(m_vol->Spot(), 0.0, K, m_vol->GetPoint(T, K), T))
      / (m_strikeEpsi * m_strikeEpsi);
  }

}