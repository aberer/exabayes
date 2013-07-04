#include "Category.hpp"
#include <cassert>
#include <algorithm>

#include "parameters/AbstractParameter.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "parameters/FrequencyParameter.hpp"
#include "parameters/ParameterContent.hpp"
#include "parameters/RateHetParameter.hpp"
#include "parameters/RevMatParameter.hpp"
#include "parameters/TopologyParameter.hpp"


namespace CategoryFuns
{
  std::string getLongName(Category cat)
  {
    switch(cat)
      {
      case Category::TOPOLOGY :
	return "Topology";
      case Category::BRANCH_LENGTHS:
	return "BranchLen" ;
      case Category::FREQUENCIES :
	return "Frequencies" ;
      case Category::SUBSTITUTION_RATES :
	return "ReversiMatr" ;
      case Category::RATE_HETEROGENEITY:
	return "RateHetero" ;
      case Category::AA_MODEL :
	return "aaModel" ;
      default : 
	{
	  assert(0); 
	  return "NOTHING"; 
	}
      }
  }


  std::string getShortName(Category cat)
  {
    switch(cat)
      {
      case Category::TOPOLOGY:
	return "topo" ;
      case Category::BRANCH_LENGTHS:
	return "bl" ;
      case Category::FREQUENCIES :
	return "pi" ;
      case Category::SUBSTITUTION_RATES:
	return "revMat";
      case Category::RATE_HETEROGENEITY :
	return "shape" ;
      case Category::AA_MODEL:
	return "aaModel" ;
      default: 
	assert(0); 
	return "NOTHING"; 
      }
  }
 
  std::string getPriorName(Category cat)
  {
    switch(cat)
      {
      case Category::TOPOLOGY: 
	return "TOPOPR"; 
      case Category::BRANCH_LENGTHS: 
	return "BRLENPR";
      case Category::FREQUENCIES: 
	return "STATEFREQPR";
      case Category::SUBSTITUTION_RATES: 
	return "REVMATPR";
      case Category::RATE_HETEROGENEITY: 
	return "SHAPEPR"; 
      default: 
	assert(0); 
      }
  } 


  std::vector<Category> getAllCategories()
  {
    return {  Category::TOPOLOGY, Category::BRANCH_LENGTHS, Category::FREQUENCIES, Category::SUBSTITUTION_RATES, Category::RATE_HETEROGENEITY} ; 
  }


  Category getCategoryByPriorName(std::string name)
  {
    // bad but this is a bit exhausting... 
    std::vector<Category> cats = getAllCategories();

    for(auto c : cats)
      {
	if(getPriorName(c).compare(name) == 0)
	  return c; 
      }

    assert(0); 
    return cats[0]; 
  }


  Category getCategoryFromLinkLabel(std::string name)
  {
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    if(name.compare("statefreq") == 0)
      return Category::FREQUENCIES; 
    else if(name.compare("ratehet") == 0)
      return Category::RATE_HETEROGENEITY; 
    else if(name.compare("revmat") == 0)
      return Category::SUBSTITUTION_RATES; 
    else if(name.compare("aamodel") == 0)
      return Category::AA_MODEL; 
    else if(name.compare("branchlength") == 0)
      return Category::BRANCH_LENGTHS; 
    else
      {
	assert(0); 
	return Category::TOPOLOGY; 
      }
  }



  std::shared_ptr<AbstractParameter> getParameterFromCategory(Category cat, nat id)
  {
    switch(cat)
      {
      case Category::TOPOLOGY :
	return std::make_shared<TopologyParameter>(id);
      case Category::BRANCH_LENGTHS:
	return std::make_shared<BranchLengthsParameter>(id); 
      case Category::FREQUENCIES :
	return std::make_shared<FrequencyParameter>(id); 
      case Category::SUBSTITUTION_RATES :
	return std::make_shared<RevMatParameter>(id); 
      case Category::RATE_HETEROGENEITY:
	return std::make_shared<RateHetParameter>(id); 
      case Category::AA_MODEL :
      default : 
	{
	  assert(0); 
	  return std::make_shared<TopologyParameter>(id);
	}
      }    
  } 
 
}
