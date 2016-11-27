// Pricer.cpp : Defines the entry point for the console application.
//

#include <pricer/stdafx.h>
#include "UTest.h"
#include <iostream>

using namespace Pricer;

int _tmain(int argc, _TCHAR* argv[])
{
  // Test the pricing by MC and PDE throw

  std::cout << ((int)std::pow(2, 16) >> 1) << std::endl;
  std::cout << ((int)std::pow(2, 16) >> 2) << std::endl;
  std::cout << ((int)std::pow(2, 16) >> 3) << std::endl;
  std::cout << ((int)std::pow(2, 16) >> 4) << std::endl;
  std::cout << ((int)std::pow(2, 16) >> 5) << std::endl;
  std::cout << ((int)std::pow(2, 16) >> 6) << std::endl;

  UTest::Test1();
  UTest::Test2();
  UTest::Test3();
  UTest::Test4();
  UTest::Test5();
  UTest::Test6();
  UTest::Test7();
  UTest::Test8();
  UTest::Test9();
  UTest::Test10();

  std::cin.get();

	return 0;

}

