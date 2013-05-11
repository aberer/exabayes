
// #include "bayes.h"

#include "Category.hpp"
#include "AbstractProposal.hpp"



/**
   @brief initializes a category of proposals. 

   @param catFreq -- must be the relative frequency (that sums up to one across all frequencies)
 */
Category::Category(string _name, category_t _type, double catFreq, vector<AbstractProposal*> _proposals )
  : type(_type)  
  , proposals(_proposals)
  , categoryFrequency(catFreq)    
  , name(_name) 
{
  // normalize relative proposal weights  
  double sum = 0;   
  for (auto proposal: proposals)
    sum += proposal->getRelativeProbability(); 

  for(auto proposal : proposals) 
    proposal->setRelativeProbability(proposal->getRelativeProbability() / sum); 
}

Category::Category(const Category &rhs)
  : type(rhs.type)
  , categoryFrequency(rhs.categoryFrequency)   
  , name(rhs.name)
{
  assert(proposals.size() == 0);  
  for(auto p : rhs.proposals)
    proposals.push_back(p->clone());     
}

Category& Category::operator=(Category rhs)
{
  swap(*this, rhs); 
  return *this; 
}


Category::~Category()
{
  for(auto p :  getProposals())
    delete p; 
}



AbstractProposal* Category::drawProposal(Randomness &rand)
{
  double r = rand.drawRandDouble01();

  for(auto proposal : proposals)
    {
      double prob = proposal->getRelativeProbability(); 
      if(r <= prob )
	return proposal; 
      else 
	r -= prob; 	
    }

  assert(0); 
  return NULL;   
}




void swap(Category &c1, Category &c2)
{
  std::swap(c1.type, c2.type); 
  std::swap(c1.name, c2.name); 
  std::swap(c1.categoryFrequency, c2.categoryFrequency); 
  std::swap(c1.proposals, c2.proposals); 
} 
