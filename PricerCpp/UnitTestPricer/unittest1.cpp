#include "stdafx.h"
#include "CppUnitTest.h"

#include <pricer/Pricer.h>

#define _SCL_SECURE_NO_WARNINGS

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestPricer
{		
  using namespace Pricer;

	TEST_CLASS(UnitTest1)
	{
  public:

    // Test 1 : Sigma constant Payoff Identity 
    // Euler Scheme, implicit scheme
		TEST_METHOD(TestMethod1)
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

      ptr<CSigmaLoc> l_sig(new CSigmaConst(0.2));
      ptr<CPayoffEUropean> l_payoff(new CPayoffId());

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "Test 1" << std::endl;
      std::cout << "Theo Value: " << l_s0 << std::endl;
      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(10.003746636298395, l_priceMC);
      Assert::AreEqual(9.9999999999999769, l_pricePDE);
			// TODO: Your test code here
		}

    // Test 2 : Sigma constant Payoff Call
    // Euler Scheme, implicit scheme
    TEST_METHOD(TestMethod2)
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

      ptr<CSigmaLoc> l_sig(new CSigmaConst(0.2));
      ptr<CPayoffEUropean> l_payoff(new CPayoffCall(l_s0));

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "Theo Value: " << CallBachelier(l_s0, 0.0, l_s0, 0.2, l_mat) << std::endl;
      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(0.079536843193572515, l_priceMC);
      Assert::AreEqual(0.079713865239046167, l_pricePDE);
    }

    // Test 3 : Sigma constant Payoff Call
    // Milstein Scheme, explicit scheme 
    // rmq : for sigma constant Milstein == Euler
    TEST_METHOD(TestMethod3)
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

      ptr<CSigmaLoc> l_sig(new CSigmaConst(0.2));
      ptr<CPayoffEUropean> l_payoff(new CPayoffCall(l_s0));

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "Theo Value: " << CallBachelier(l_s0, 0.0, l_s0, 0.2, l_mat) << std::endl;
      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(0.079536843193572515, l_priceMC);
      Assert::AreEqual(0.079823408639291368, l_pricePDE);
    }

    // Test 4 : Sigma constant Payoff Digit
    // Milstein Scheme, implicit scheme
    TEST_METHOD(TestMethod4)
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

      ptr<CSigmaLoc> l_sig(new CSigmaConst(0.2));
      ptr<CPayoffEUropean> l_payoff(new CPayoffDigit(l_s0));

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(50.0, l_s0, l_mat, false);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "Theo Value: " << 0.5 << std::endl;
      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(0.52000000000000002, l_priceMC);
      Assert::AreEqual(0.49999999999999611, l_pricePDE);
    }

    // Test 5 : Sigma LN shifted Payoff Identity 
    // Milstein Scheme, implicit scheme
    TEST_METHOD(TestMethod5)
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

      ptr<CSigmaLoc> l_sig(new CSigmaLNShifted(0.2, 0.5));
      ptr<CPayoffEUropean> l_payoff(new CPayoffId());

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "Theo Value: " << 10.0 << std::endl;
      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(10.043414891676351, l_priceMC);
      Assert::AreEqual(10.000000000000000, l_pricePDE);
    }

    // Test 6 : Sigma LN shifted Payoff Call
    // Milstein Scheme, implicit scheme
    TEST_METHOD(TestMethod6)
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
      ptr<CSigmaLoc> l_sig(new CSigmaLNShifted(0.2, 0.5));
      ptr<CPayoffEUropean> l_payoff(new CPayoffCall(l_s0 - 1.0));

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "Theo Value: " << CallBlackShifted(l_s0, 0.5, 0.0, l_s0 - 1.0, 0.8, l_mat) << std::endl;
      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(1.3854367986049263, l_priceMC);
      Assert::AreEqual(1.3971814031761520, l_pricePDE);
    }

    // Test 7 : Sigma LN shifted Payoff Digit
    // Milstein Scheme, implicit scheme
    TEST_METHOD(TestMethod7)
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

      ptr<CSigmaLoc> l_sig(new CSigmaLNShifted(0.2, 0.3));
      ptr<CPayoffEUropean> l_payoff(new CPayoffDigit(l_s0 + 1.0));

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "Theo Value: " << CallDigitBlackShifted(l_s0, 0.3, 0.0, l_s0 + 1.0, 0.2, l_mat) << std::endl;
      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(0.28999999999999998, l_priceMC);
      Assert::AreEqual(0.28779485688267609, l_pricePDE);
    }

    // Test 8 : Sigma CEV (beta < 1) Payoff Identity
    // Milstein Scheme, implicit scheme
    TEST_METHOD(TestMethod8)
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

      ptr<CSigmaLoc> l_sig(new CSigmaCEV(0.2, 0.5));
      ptr<CPayoffEUropean> l_payoff(new CPayoffId());

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(9.9931294854125650, l_priceMC);
      Assert::AreEqual(10.000000000000000, l_pricePDE);
    }

    // Test 9 : Sigma CEV (beta < 1)  Payoff Call
    // Milstein Scheme, implicit scheme
    TEST_METHOD(TestMethod9)
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

      ptr<CSigmaLoc> l_sig(new CSigmaCEV(0.8, 0.5));
      ptr<CPayoffEUropean> l_payoff(new CPayoffCall(l_s0 - 1.0));

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(1.5528401627364483, l_priceMC);
      Assert::AreEqual(1.5663962946100614, l_pricePDE);
    }

    // Test 10 : Sigma CEV (beta > 1)  Payoff Digit
    // Milstein Scheme, implicit scheme
    TEST_METHOD(TestMethod10)
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

      ptr<CSigmaLoc> l_sig(new CSigmaCEV(0.8, 1.3));
      ptr<CPayoffEUropean> l_payoff(new CPayoffDigit(l_s0 - 1.0));

      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      ptr<process> l_spProcess = std::make_shared<process>(l_sig, l_isMilstein);

      // MC Part
      ptr<timeSteps> l_spMcTS = std::make_shared<timeSteps>(l_nbTimeStep, l_mat / l_nbTimeStep);
      stat    l_stat(l_payoff, l_spMcTS);
      pathGen l_pathGen(std::static_pointer_cast<IProcess>(l_spProcess), CNormGen(), l_spMcTS, l_nbSimu);
      l_pathGen.GenSequence();
      l_pathGen.FillStat(l_s0, l_stat);

      // PDE part;
      ptr<timeSteps> l_spPdeTS = std::make_shared<timeSteps>(l_sizeInT, l_mat / l_sizeInT);
      pdeEng l_pdeEng(std::make_shared<process>(l_sig), l_payoff, *l_spPdeTS, 0.0, 100.0 + l_s0, l_theta);
      l_pdeEng.OverloadBoundary(25.0, l_s0, l_mat, true);
      l_pdeEng.InitH(l_sizeInH);
      std::vector<double> l_res = l_pdeEng.Compute();
      Linear_interp l_prices(l_pdeEng.UnderValues(), l_res);

      double l_priceMC = l_stat.Price();
      double l_pricePDE = l_prices.interp(l_s0);

      std::cout << "MC  result: " << l_priceMC << std::endl;
      std::cout << "PDE result: " << l_pricePDE << std::endl;
      std::cout << std::endl;

      Assert::AreEqual(0.16339999999999999, l_priceMC);
      Assert::AreEqual(0.16876916735316458, l_pricePDE);
    }

    TEST_METHOD(TestMethod11)
    {
      typedef CStatEurPriceOnly stat;
      typedef CProcessVolLoc process;
      typedef CPathGen<CNormGen> pathGen;
      typedef CPdeEng pdeEng;

      double l_s0     = 0.0;
      double l_strike = 0.2;
      double l_mat    = 1.0;

      // MC Part
      bool l_isMilstein   = false;
      size_t l_nbSimu     = 1000;
      size_t l_nbTimeStep = 100;

      // PDE Part
      double l_theta   = 0.0;
      size_t l_sizeInH = 100;
      size_t l_sizeInT = 100;

      // Dupire rep
      unsigned int l_nInTime = 101;
      unsigned int l_nInSpace = 101;

      // Creation of a volatility
      std::vector<double> l_oTenors = { 0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0 };
      std::vector<double> l_oStrikes = { -50.0, 0.0, 50.0 };

      // Vol flat
      double l_sigma = 0.15;
      std::vector<double> l_oSmile(l_oStrikes.size(), l_sigma);
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

      double l_priceMC      = l_stat.Price();
      double l_priceMC_m10  = l_stat_m10.Price();
      double l_priceMC_p10  = l_stat_p10.Price();
      double l_pricePDE     = l_prices.interp(l_s0);
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
      std::cout << std::endl;
    }


	};
}