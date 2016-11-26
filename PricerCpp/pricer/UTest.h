#ifndef pricer_ctest
#define pricer_ctest


namespace Pricer
{
  // Dimmy class to do test
  // kind of dirty should create an other project
  // to do those test
  namespace UTest
  {
    // Test 1 : Sigma constant Payoff Identity 
    // Euler Scheme, implicit scheme
    void Test1();
    // Test 2 : Sigma constant Payoff Call
    // Euler Scheme, implicit scheme
    void Test2();
    // Test 3 : Sigma constant Payoff Call
    // Milstein Scheme, explicit scheme 
    // rmq : for sigma constant Milstein == Euler
    void Test3();
    // Test 4 : Sigma constant Payoff Digit
    // Milstein Scheme, implicit scheme
    void Test4();
    // Test 5 : Sigma LN shifted Payoff Identity 
    // Milstein Scheme, implicit scheme
    void Test5();
    // Test 6 : Sigma LN shifted Payoff Call
    // Milstein Scheme, implicit scheme
    void Test6();
    // Test 7 : Sigma LN shifted Payoff Digit
    // Milstein Scheme, implicit scheme
    void Test7();
    // Test 8 : Sigma CEV (beta < 1) Payoff Identity
    // Milstein Scheme, implicit scheme
    void Test8();
    // Test 9 : Sigma CEV (beta < 1)  Payoff Call
    // Milstein Scheme, implicit scheme
    void Test9();
    // Test 10 : Sigma CEV (beta > 1)  Payoff Digit
    // Milstein Scheme, implicit scheme
    void Test10();

  };

}


#endif
