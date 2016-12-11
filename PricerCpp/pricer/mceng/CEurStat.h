#ifndef pricer_ceurstat
#define pricer_ceurstat

#include <vector>

#include <pricer/UType.h>
#include <pricer/instrs/CPayoff.h>

namespace Pricer
{
  // TODO add interval error price
  // info before full pricing

  class CStat
  {
  public:
    CStat(const ptr<timeSteps>& p_spTimeSteps, size_t p_iNbSimu);
    virtual ~CStat() {};

    virtual double Price() const = 0;
    void Add(ptr<path> p_spPath);

  protected:
    ptr<timeSteps> m_spTimeSteps;
    paths          m_paths;
  };

  class CStatEurPriceOnly : public CStat
  {
  public:
    CStatEurPriceOnly(const ptr<CPayoffEUropean>& p_payoff, 
      const ptr<timeSteps>& p_spTimeSteps, size_t p_iNbSimu = 1000);
    //~CStatEurPriceOnly();

    double Price() const;

  private:
    ptr<CPayoffEUropean> m_payoff;
  };

  // TODO question: Where to put the intelligence of
  // selected the point in the path that are used for
  // the pricing. In first it'seems logical to give this 
  // logic only to the CPayoff. We choose to put part of 
  // this intelligence in the CStatPathDepend in order
  // to determine which point in the path to select only
  // once. 
  class CStatPathDepend : public CStat
  {
  public:
    CStatPathDepend(const ptr<CPayoffPathDepPartial>& p_payoff,
      const ptr<timeSteps>& p_spTimeSteps, size_t p_iNbSimu = 1000);
    //~CStatEurPriceOnly();

    double Price() const;

  private:
    ptr<CPayoffPathDepPartial> m_payoff;

    std::vector<size_t> pointToSelect() const;
  };
}

#endif