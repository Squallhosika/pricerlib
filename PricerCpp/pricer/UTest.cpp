#include <pricer/stdafx.h>
#include "UTest.h"

#include "pricer/Pricer.h"

namespace Pricer
{
  namespace UTest
  {
    void Test1()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = false;
      size_t l_nbSimu = 1000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 100;
      size_t l_sizeInT = 100;

      typedef CPayoffId PAYOFF;
      ptr<CSigmaConst> l_sig(new CSigmaConst(0.2));
      PAYOFF l_payoff;
      //PAYOFF   l_payoff;

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc  process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 1" << std::endl;
      std::cout << "Theo Value: " << l_s0 << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test2()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = false;
      size_t l_nbSimu = 10000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 1000;
      size_t l_sizeInT = 100;

      typedef CPayoffCall PAYOFF;
      ptr<CSigmaConst> l_sig(new CSigmaConst(0.2));
      PAYOFF l_payoff(l_s0);
      //PAYOFF   l_payoff;

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc  process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 2" << std::endl;
      std::cout << "Theo Value: " << CallBachelier(l_s0, 0.0, l_s0, 0.2, l_mat) << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test3()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = true;
      size_t l_nbSimu = 10000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 1.0;
      size_t l_sizeInH = 1000;
      size_t l_sizeInT = 1000;

      typedef CPayoffCall PAYOFF;
      ptr<CSigmaConst> l_sig(new CSigmaConst(0.2));
      PAYOFF l_payoff(l_s0);
      //PAYOFF   l_payoff;

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc  process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 3" << std::endl;
      std::cout << "Theo Value: " << CallBachelier(l_s0, 0.0, l_s0, 0.2, l_mat) << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test4()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = true;
      size_t l_nbSimu = 1000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 100;
      size_t l_sizeInT = 100;

      typedef CPayoffDigit PAYOFF;
      ptr<CSigmaConst> l_sig(new CSigmaConst(0.2));
      PAYOFF l_payoff(l_s0);
      //PAYOFF   l_payoff;

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc  process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 4" << std::endl;
      std::cout << "Theo Value: " << 0.5 << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test5()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = true;
      size_t l_nbSimu = 1000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 100;
      size_t l_sizeInT = 100;

      typedef CPayoffId PAYOFF;
      ptr<CSigmaLNShifted> l_sig(new CSigmaLNShifted(0.2, 0.5));
      PAYOFF l_payoff;

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc  process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 5" << std::endl;
      std::cout << "Theo Value: " << 10.0 << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test6()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = true;
      size_t l_nbSimu = 10000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 100;
      size_t l_sizeInT = 100;

      typedef CPayoffCall        PAYOFF;
      ptr<CSigmaLNShifted> l_sig(new CSigmaLNShifted(0.2, 0.5));
      PAYOFF l_payoff(l_s0 - 1.0);

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 6" << std::endl;
      std::cout << "Theo Value: " << CallBlackShifted(l_s0, 0.5, 0.0, l_s0 - 1.0, 0.8, l_mat) << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test7()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = true;
      size_t l_nbSimu = 1000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 1000;
      size_t l_sizeInT = 100;

      typedef CPayoffDigit      PAYOFF;
      ptr<CSigmaLNShifted> l_sig(new CSigmaLNShifted(0.2, 0.3));
      PAYOFF l_payoff(l_s0 + 1.0);

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc  process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 7" << std::endl;
      std::cout << "Theo Value: " << CallDigitBlackShifted(l_s0, 0.3, 0.0, l_s0 + 1.0, 0.2, l_mat) << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test8()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = true;
      size_t l_nbSimu = 10000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 100;
      size_t l_sizeInT = 100;

      typedef CPayoffId PAYOFF;
      ptr<CSigmaCEV> l_sig(new CSigmaCEV(0.2, 0.5));
      PAYOFF l_payoff;

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 8" << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test9()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = true;
      size_t l_nbSimu = 10000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 100;
      size_t l_sizeInT = 100;

      typedef CSigmaCEV  SIGMA;
      typedef CPayoffCall  PAYOFF;
      ptr<CSigmaCEV> l_sig(new CSigmaCEV(0.8, 0.5));
      PAYOFF l_payoff(l_s0 - 1.0);

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 9" << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test10()
    {
      double l_s0 = 10.0;
      double l_mat = 1.0;

      // MC Part
      bool l_isMilstein = true;
      size_t l_nbSimu = 10000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta = 0.0;
      size_t l_sizeInH = 1000;
      size_t l_sizeInT = 1000;

      typedef CSigmaCEV   SIGMA;
      typedef CPayoffDigit PAYOFF;
      ptr<CSigmaCEV> l_sig(new CSigmaCEV(0.8, 1.3));
      PAYOFF l_payoff(l_s0 - 1.0);

      typedef CStatEurPriceOnly<PAYOFF> stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen, stat> pathGen;
      typedef CPdeEng<PAYOFF> pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<stat::timeSteps> l_spMcTS = std::make_shared<stat::timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<stat::timeSteps> l_spPdeTS = std::make_shared<stat::timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      std::cout << "Test 10" << std::endl;
      std::cout << "MC  result: " << l_stat.Price() << std::endl;
      std::cout << "PDE result: " << l_prices.interp(l_s0) << std::endl;
      std::cout << std::endl;
    }

    void Test11()
    {
      std::vector<double> l_oTenor = { 0.0, 0.1, 0.25, 0.5, 0.1 };
      std::vector<double> l_oStrike = { 0.1, 0.5, 1.0, 2.0, 10.0 };

      double l_sigma = 0.15;
      std::vector<double> l_oSmile(l_oStrike.size(), l_sigma);
      DMatrix l_oSurface(l_oTenor.size(), l_oSmile);


    }
  }
}