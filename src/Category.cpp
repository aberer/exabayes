
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


Category::~Category()
{
}
