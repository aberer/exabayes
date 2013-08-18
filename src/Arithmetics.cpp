
#include "Arithmetics.hpp"
#include <cmath>

namespace Arithmetics
{
  std::pair<double,double> getMeanAndVar (const std::vector<double> &data )
  {
    double mean = 0; 
    for(auto d: data)
      mean += d; 
    mean /= data.size(); 

    double var = 0; 
    for(auto d : data)
      var += pow(d - mean, 2); 
    var /= data.size(); 

    return std::pair<double,double>(mean,var);  
  }
}
