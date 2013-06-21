#ifndef _INSERTION_SCORE 
#define _INSERTION_SCORE 


#include "Branch.hpp"

class  InsertionScore
{
public: 
  InsertionScore(Branch _b, vector<nat> _tmp) : b(_b), partitionParsimony(_tmp){}  
  Branch getBranch() const  {return b; }

  double getWeight() const {return  logProb; }
  void setWeight(double w) { logProb = w; }

  nat getScore() const
  {
    nat result = 0; 
    for(auto b : partitionParsimony)
      result += b; 
    return result; 
  }

  nat getPartitionScore(int model) const{return partitionParsimony[model] ; }


private: 
  Branch b; 
  vector<nat> partitionParsimony; 
  double logProb;

  friend ostream& operator<< (ostream &out, const InsertionScore &rhs) { 
    out <<  "(" << rhs.b << "=" ; 
    for(auto elem : rhs.partitionParsimony)
      out << elem << "," ; 
    return out; 
  }
} ; 


#endif
