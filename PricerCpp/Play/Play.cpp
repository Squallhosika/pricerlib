#include "stdafx.h"
#include <pricer/Pricer.h>
#include <pricer/math/CVolInverser.h>

using namespace Pricer;

int _tmain(int argc, _TCHAR* argv[])
{
  // Test the pricing by MC and PDE throw
  typedef CStatEurPriceOnly stat;
  typedef CProcessVolLoc process;
  typedef CPathGen<CNormGen> pathGen;
  typedef CPdeEng pdeEng;

  double l_s0 = 0.0;
  double l_strike = 0.0;
  double l_mat = 1.0;

  // MC Part
  bool l_isMilstein = false;
  size_t l_nbSimu = 10000;
  size_t l_nbTimeStep = 100;

  // PDE Part
  double l_theta = 0.0;
  size_t l_sizeInT = 100;
  size_t l_sizeInH = 500;
  
  // Dupire rep
  unsigned int l_nInTime = 101;
  unsigned int l_nInSpace = 1001;

  // bound 
  double l_dVolMin = 0.001;
  double l_dVolMax = 2.0;


  // Creation of a volatility
  std::vector<double> l_oTenors = { 0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0 };
  std::vector<double> l_oStrikes = { -5.0, -0.2, -0.1 , 0.0, 0.1, 0.2, 5.0 };

  // Vol flat
  //double l_sigma = 0.15;
  //std::vector<double> l_oSmile(l_oStrikes.size(), l_sigma);
  //DMatrix l_oSurface(l_oTenors.size(), l_oSmile);

  // Vol with smile only
  // std::vector<double> l_oSmile = { 0.25, 0.23, 0.19, 0.15, 0.19, 0.23, 0.25 };
  std::vector<double> l_oSmile = { 0.10, 0.11, 0.13, 0.15, 0.17, 0.19, 0.20 };
  DMatrix l_oSurface(l_oTenors.size(), l_oSmile);

  ptr<CVol> l_spVolSquare(new CVolSquare(l_s0, l_oTenors, l_oStrikes, l_oSurface));

  std::vector<double> l_oTimeSteps(l_nInTime);
  std::vector<double> l_oSpaceSteps(l_nInSpace);
  std::generate(l_oTimeSteps.begin(), l_oTimeSteps.end(),
    SGener_iter(l_oTenors.front(), l_oTenors.back(), l_nInTime));
  std::generate(l_oSpaceSteps.begin(), l_oSpaceSteps.end(),
    SGener_iter(l_oStrikes.front(), l_oStrikes.back(), l_nInSpace));
  ptr<CGrid> l_spGrid(new CGrid(l_oTimeSteps, l_oSpaceSteps));

  // Generate the dupire local volatility
  ptr<CSigmaDupire> l_sigDupire(new CSigmaDupire(l_spVolSquare));
  l_sigDupire->Init(l_spGrid);
  ptr<CSigmaLoc> l_sig = l_sigDupire;

  double l_shift = 0.1;

  // ptr<CPayoffEUropean> l_payoff(new CPayoffId());
  ptr<CPayoffEUropean> l_payoff(new CPayoffCall(l_strike));
  ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

  // MC Part
  ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
  stat    l_stat(l_payoff, l_spMcTS);
  stat    l_stat_m10(l_payoff, l_spMcTS);
  stat    l_stat_p10(l_payoff, l_spMcTS);
  pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
  l_pathGen.GenSequence();
  l_pathGen.FillStat(l_s0, l_stat);
  l_pathGen.FillStat(l_s0 - l_shift, l_stat_m10);
  l_pathGen.FillStat(l_s0 + l_shift, l_stat_p10);

  // PDE part;
  ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
  pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS,
    l_s0 - 50.0, l_s0 + 50.0, l_theta);
  l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
  std::cout << "m_underMin: " << l_pdeEng.m_underMin << std::endl;
  std::cout << "m_underMax: " << l_pdeEng.m_underMax << std::endl;
  l_pdeEng.InitH(l_sizeInH);
  std::vector<double> l_res = l_pdeEng.Compute();
  Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

  double l_priceMC = l_stat.Price();
  double l_priceMC_m10 = l_stat_m10.Price();
  double l_priceMC_p10 = l_stat_p10.Price();
  double l_pricePDE = l_prices.interp(l_s0);
  double l_pricePDE_m10 = l_prices.interp(l_s0 - l_shift);
  double l_pricePDE_p10 = l_prices.interp(l_s0 + l_shift);

  std::cout << "Test 11" << std::endl;
  std::cout << "Theo Value: " << CallBachelier(l_s0, 0.0, l_strike,
    l_spVolSquare->GetPoint(l_mat, l_strike), l_mat) << std::endl;
  std::cout << "Theo Value -10: " << CallBachelier(l_s0 - l_shift, 0.0, l_strike,
    l_spVolSquare->GetPoint(l_mat, l_strike), l_mat) << std::endl;
  std::cout << "Theo Value +10: " << CallBachelier(l_s0 + l_shift, 0.0, l_strike,
    l_spVolSquare->GetPoint(l_mat, l_strike), l_mat) << std::endl;
  std::cout << "MC result: " << l_priceMC << std::endl;
  std::cout << "MC -10 result: " << l_priceMC_m10 << std::endl;
  std::cout << "MC +10 result: " << l_priceMC_p10 << std::endl;
  std::cout << "PDE result: " << l_pricePDE << std::endl;
  std::cout << "PDE -10 result: " << l_pricePDE_m10 << std::endl;
  std::cout << "PDE +10 result: " << l_pricePDE_p10 << std::endl;

  std::cout << "Theo Vol: " << UVolInverser::InversCallBach(l_s0, 0.0, l_strike,
    CallBachelier(l_s0, 0.0, l_strike, l_spVolSquare->GetPoint(l_mat, l_strike), l_mat),
    l_mat, l_dVolMin, l_dVolMax) << std::endl;
  std::cout << "Theo Vol -10: " << UVolInverser::InversCallBach(l_s0 - l_shift, 0.0, l_strike,
    CallBachelier(l_s0 - l_shift, 0.0, l_strike, l_spVolSquare->GetPoint(l_mat, l_strike - l_shift), l_mat),
    l_mat, l_dVolMin, l_dVolMax) << std::endl;
  std::cout << "Theo Vol +10: " << UVolInverser::InversCallBach(l_s0 + l_shift, 0.0, l_strike,
    CallBachelier(l_s0 + l_shift, 0.0, l_strike, l_spVolSquare->GetPoint(l_mat, l_strike + l_shift), l_mat),
    l_mat, l_dVolMin, l_dVolMax) << std::endl;

  std::cout << "MC result: " << UVolInverser::InversCallBach(l_s0, 0.0, l_strike,
    l_priceMC, l_mat, l_dVolMin, l_dVolMax)   << std::endl;
  std::cout << "MC -10 result: " << UVolInverser::InversCallBach(l_s0 - l_shift, 0.0, l_strike,
    l_priceMC_m10, l_mat, l_dVolMin, l_dVolMax) << std::endl;
  std::cout << "MC +10 result: " << UVolInverser::InversCallBach(l_s0 + l_shift, 0.0, l_strike,
    l_priceMC_p10, l_mat, l_dVolMin, l_dVolMax) << std::endl;
  std::cout << "PDE result: " << UVolInverser::InversCallBach(l_s0, 0.0, l_strike,
    l_pricePDE, l_mat, l_dVolMin, l_dVolMax) << std::endl;
  std::cout << "PDE -10 result: " << UVolInverser::InversCallBach(l_s0 - l_shift, 0.0, l_strike,
    l_pricePDE_m10, l_mat, l_dVolMin, l_dVolMax) << std::endl;
  std::cout << "PDE +10 result: " << UVolInverser::InversCallBach(l_s0 + l_shift, 0.0, l_strike,
    l_pricePDE_p10, l_mat, l_dVolMin, l_dVolMax) << std::endl;

  std::cout << std::endl;

  std::cout << "END NOW COME THE VOL AND VAR PAYOFF" << std::endl;

  ptr<timeSteps> l_oTrueTS(new timeSteps()) ;

  l_oTrueTS->resize(l_spMcTS->size() + 1);
  (*l_oTrueTS)[0] = 0.0;
  std::partial_sum(l_spMcTS->begin(), l_spMcTS->end(), ++l_oTrueTS->begin());


  ptr<CPayoffPathDepPartial> l_payoffVol(new CPayoffVolBach(*l_oTrueTS));
  ptr<CPayoffPathDepPartial> l_payoffVar(new CPayoffVarBach(*l_oTrueTS));
  
  CStatPathDepend l_statVol(l_payoffVol, l_oTrueTS);
  CStatPathDepend l_statVar(l_payoffVar, l_oTrueTS);

  l_pathGen.FillStat(l_s0, l_statVol);
  l_pathGen.FillStat(l_s0, l_statVar);

  double l_priceVol = l_statVol.Price();
  double l_priceVar = l_statVar.Price();

  std::cout << "Theo Value: " << CallBachelier(l_s0, 0.0, l_strike,
    l_spVolSquare->GetPoint(l_mat, l_strike), l_mat) << std::endl;
  std::cout << "MC  l_priceVol: " << l_priceVol << std::endl;
  std::cout << "MC l_priceVar: " << l_priceVar << std::endl;
  std::cout << "MC l_priceVar normalize: " << std::sqrt(l_priceVar) << std::endl;
  std::cout << std::endl;

  // std::cin.get();

  return 0;

}