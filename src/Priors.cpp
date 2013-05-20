
#include "Priors.hpp"


 double DirichletPrior::getLogProb(vector<double> values) const 
{       
  // cout << "values are "; 
  // for(auto v : values)
  //   cout <<  v << "," ; 
  // cout << endl; 

  double result = densityDirichletWrapper(values, alphas, true ); 

  // result = log(result);
  // cout << "RESULT is " <<  result << endl; 
    
  return result; 
}
