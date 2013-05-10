
#ifndef __CATEGORY_H
#define __CATEGORY_H

#include <vector>

typedef struct _pfun proposalFunction; 

#include "categoryType.h"
#include "Randomness.hpp"

class AbstractProposal; 

class Category
{
public: 
  Category(string name, category_t _type, double catFreq, vector<AbstractProposal*> _proposals ); 
  ~Category(); 

  void copyDeep(Category& rhs); 
  
  /** @brief adds a proposal */ 
  void addProposal(AbstractProposal *proposal){proposals.push_back(proposal); }

  /** @brief samples according to relative frequency of the proposal */ 
  AbstractProposal* drawProposal(Randomness &rand);

  void setCatFreq(double freq){categoryFrequency = freq; }
  double getCatFreq( ) const {return categoryFrequency ; }
  string getName()  const {return name; }

  vector<AbstractProposal*> getProposals(){return proposals; }


private: 
  category_t type;
  vector<AbstractProposal*> proposals; 
  double categoryFrequency;   
  string name; 
};  


#endif
