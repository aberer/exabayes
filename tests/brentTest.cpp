#include "math/brent.hpp"


class OptFunc : public brent::func_base
{
public:   
  OptFunc(std::function< double(double)> fun)
    : _fun(fun)
  {
  }

  virtual double operator() (double val)
  {
    return _fun(val); 
  }

private: 
  std::function< double(double) > _fun; 

}; 



TEST(BRENT_TEST, general)
{

  auto square =  [](double x ) -> double 
  {
    return 3 * pow(x -2 ,2) + 2 ; 
  }; 

  auto optFun = OptFunc(square); 
  
  auto result = brent::zero(-123,345, 1e-3, optFun); 
  tout << SHOW(result) << std::endl; 

}
