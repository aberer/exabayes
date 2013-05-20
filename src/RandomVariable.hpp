/**
   @file Parameter.hpp

   @brief models (free) parameter we want to integrate   
 */ 


#ifndef _RANDOM_VARIABLE_H
#define _RANDOM_VARIABLE_H

#include <memory> 
#include <vector>

#include "Priors.hpp"
#include "GlobalVariables.hpp"

using namespace std; 

class RandomVariable
{
public: 
  RandomVariable(category_t cat, nat id)
    : cat(cat) , id(id)  {} 

  void setPrior (shared_ptr<AbstractPrior> _prior)
  {
    prior = _prior; 
  }
  
  void addPartition(nat id) { partitions.push_back(id); }
  const vector<nat>& getPartitions() const {return partitions; }
  
private: 
  category_t cat; 
  vector<nat> partitions; 	// the partitions, this parameter
				// applies to
  nat id ; 			//  which id does the  parameter have among all of its kind  

  shared_ptr<AbstractPrior> prior;   

  friend ostream& operator<<(ostream &out, const RandomVariable& rhs)
  {
    bool  isFirst = true; 
    out << "{" ; 
    for(auto v : rhs.partitions)
      {
	if (not isFirst) 
	  out << "," ; 
	else 
	  isFirst = false; 
	out << v; 
      }
    out << "}"; 
    return out;
  }
};

#endif
