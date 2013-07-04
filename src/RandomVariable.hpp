/**
   @file Parameter.hpp

   @brief models (free) parameter we want to integrate   
 */ 


#ifndef _RANDOM_VARIABLE_H
#define _RANDOM_VARIABLE_H

#include <memory> 
#include <vector>
#include <map>

#include "Priors.hpp"
#include "GlobalVariables.hpp"
#include "Category.hpp"


class RandomVariable
{
public: 
  RandomVariable(Category cat, nat id)
    : cat(cat) , id(id)  {} 

  void setPrior (std::shared_ptr<AbstractPrior> _prior) { prior = _prior; }  
  void addPartition(nat id) { partitions.push_back(id); }
  const std::vector<nat>& getPartitions() const {return partitions; }
  Category getCategory() const { return cat; }
  nat getId() const { return id; }
  void setId(nat _id ) {id = _id; }
  std::shared_ptr<AbstractPrior> getPrior() const {return prior; }
  std::ostream&  printShort(std::ostream& out); 

  // make cloneable, if there is information in it   
  // RandomVariablePtr clone() const ; 


  RandomVariable(RandomVariable& rhs) = delete ; 
  RandomVariable& operator=(RandomVariable &rhs) = delete; 

private: 
  Category cat; 
  std::vector<nat> partitions; 	// the partitions, this parameter
				// applies to
  nat id ; 			//  which id does the  parameter have among all of its kind  
  std::shared_ptr<AbstractPrior> prior; 

  friend std::ostream& operator<<(std::ostream &out, const RandomVariable& rhs); 
};




#endif
