#ifndef _ARRAYPOLICY_HPP
#define _ARRAYPOLICY_HPP

#include <vector>
#include <memory>

#include "Branch.hpp"
#include "ArrayOrientation.hpp"

class TreeAln; 

class ArrayPolicy
{
public: 
  virtual void imprintPolicy(const TreeAln &traln, ArrayOrientation &arrayOrient)  = 0; 

  void imprint(const TreeAln &traln, ArrayOrientation &arrayOrient) 
  {
    prevLnl = traln.getTrHandle().likelihood; 
    partitionLnls = traln.getPartitionLnls(); 
    imprintPolicy(traln, arrayOrient); 
  }

  virtual void freeMemory() = 0;  


  /**
     @brief deal with rejection
   */ 
  void accountForRejection(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat>& invalidNodes, ArrayOrientation &arrayOrient) 
  {
    accountForRejectionPolicy(traln, partitions, invalidNodes, arrayOrient); 
   
    auto lnls = traln.getPartitionLnls(); 
    auto sum = 0.; 
    for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i)
      {
	if(partitions[i])
	  lnls[i] =  partitionLnls[i] ; 
	sum += lnls[i]; 
      }
    traln.setPartitionLnls(lnls); 
    traln.getTrHandle().likelihood = sum;  
  }

  virtual void accountForRejectionPolicy(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat>& invalidNodes, ArrayOrientation &arrayOrient)  = 0; 


  virtual std::unique_ptr<ArrayPolicy> clone() const = 0; 
  virtual void prepareForEvaluation(TreeAln &traln, BranchPlain virtualRoot, nat models, ArrayOrientation& arrayOrientation ) = 0; 

  virtual void enableRestoreGapVector(){}

protected: 
  double prevLnl; 
  std::vector<double> partitionLnls;
}; 

#endif
