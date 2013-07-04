#include "RevMatParameter.hpp"
#include "axml.h"

void RevMatParameter::applyParameter(TreeAln& traln, std::vector<nat> model, ParameterContent &content) const
{
  for(auto &m : model)
    {
      pInfo *partition = traln.getPartition(m);
      assert(content.values.size() == numStateToNumInTriangleMatrix(partition->states)); 

      for(nat i = 0; i < numStateToNumInTriangleMatrix(partition->states); ++i)
	{
	  // partition->substRates
	}
    }  
} 


ParameterContent RevMatParameter::extractParameter(const TreeAln &traln, vector<nat> model)  const
{
  return ParameterContent();
}   

