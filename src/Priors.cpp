#include "Priors.hpp"

double DirichletPrior::getLogProb(vector<double> values) const 
{       
  assert(values.size() == alphas.size() ); 

  // cout << "values are "; 
  // for(auto v : values)
  //   cout <<  v << "," ; 
  // cout << endl; 

  double result = densityDirichletWrapper(values, alphas, true ); 

  // cout << "RESULT is " <<  result << endl; 
    
  return result; 
}
