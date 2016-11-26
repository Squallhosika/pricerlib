#ifndef pricer_ctridagmatrix
#define pricer_ctridagmatrix

#include <vector>

namespace Pricer
{
  // //LU decomposition for the tridiagonal matrix extract from the numerical recipes
  //template<typename T>
  //void tridag(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c, 
  //  const std::vector<T> &r, std::vector<T> &u)
  //{
  //  int j, n = a.size();
  //  T bet;
  //  std::vector<T> gam(n);
  //  if (b[0] == 0.0) throw("Error 1 in tridag");
  //  u[0] = r[0] / (bet = b[0]);
  //  for (j = 1; j<n; j++) {
  //    gam[j] = c[j - 1] / bet;
  //    bet = b[j] - a[j] * gam[j];
  //    if (bet == 0.0) throw("Error 2 in tridag");
  //    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  //  }
  //  for (j = (n - 2); j >= 0; j--)
  //    u[j] -= gam[j + 1] * u[j + 1];
  //}

  //LU decomposition for the tridiagonal matrix extract from the numerical recipes
  template<typename T>
  void tridag(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c,
    const std::vector<T> &r, std::vector<T> &u)
  {
    int j, n = b.size();
    T bet;
    std::vector<T> gam(n);
    if (b[0] == 0.0) throw("Error 1 in tridag");
    u[0] = r[0] / (bet = b[0]);
    for (j = 1; j<n; j++) {
      gam[j] = c[j - 1] / bet;
      bet = b[j] - a[j - 1] * gam[j];
      if (bet == 0.0) throw("Error 2 in tridag");
      u[j] = (r[j] - a[j - 1] * u[j - 1]) / bet;
    }
    for (j = (n - 2); j >= 0; j--)
      u[j] -= gam[j + 1] * u[j + 1];
  }

  template<typename T>
  class CTriDagMatrix
  {
  public:
    CTriDagMatrix(size_t p_size);
    //~CProcessVolLoc();
    CTriDagMatrix(const std::vector<T>& p_sup, const std::vector<T>& p_cent, const std::vector<T>& p_inf);

    inline size_t Size() const;
    std::vector<T> LUDemcop(const std::vector<T>& p_r) const;
    //TODO copy of the vector or not ???
    void Fill(std::vector<T>& p_sup, std::vector<T>& p_cent, std::vector<T>& p_inf);

  private:
    std::vector<T> m_sup;
    std::vector<T> m_cent;
    std::vector<T> m_inf;

  };

  template<typename T>
  CTriDagMatrix<T>::CTriDagMatrix(size_t p_size)
    :m_sup(p_size - 1), m_cent(p_size), m_inf(p_size - 1)
  {
  }

  template<typename T>
  CTriDagMatrix<T>::CTriDagMatrix(const std::vector<T>& p_sup, const std::vector<T>& p_cent, const std::vector<T>& p_inf)
    : m_sup(p_sup), m_cent(p_cent), m_inf(p_inf)
  {
  }

  template<typename T>
  inline size_t CTriDagMatrix<T>::Size() const
  {
    return m_cent.size();
  }

  template<typename T>
  std::vector<T> CTriDagMatrix<T>::LUDemcop(const std::vector<T>& p_r) const
  {
    // TODO add check size exception
    std::vector<T> l_res(Size());
    tridag(m_inf, m_cent, m_sup, p_r, l_res);

    return l_res;
  }

  template<typename T>
  void CTriDagMatrix<T>::Fill(std::vector<T>& p_sup, std::vector<T>& p_cent, std::vector<T>& p_inf)
  {
    m_sup  = p_sup;
    m_cent = p_cent;
    m_inf  = p_inf;
  }

}
#endif
