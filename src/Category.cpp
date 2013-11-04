#include "Category.hpp"
#include <cassert>
#include <memory>
#include <algorithm>

#include "parameters/ProtModelParameter.hpp"
#include "parameters/AbstractParameter.hpp"
#include "parameters/BranchLengthsParameter.hpp"
#include "parameters/FrequencyParameter.hpp" 
#include "parameters/ParameterContent.hpp"
#include "parameters/RateHetParameter.hpp"
#include "parameters/RevMatParameter.hpp"
#include "parameters/TopologyParameter.hpp"


std::ostream&  operator<<(std::ostream& out, const Category &rhs)
{
  return out << CategoryFuns::getLongName(rhs) ; 
}


namespace CategoryFuns
{
  std::string getLongName(Category cat)
  {
    switch(cat)
      {
      case Category::TOPOLOGY :
	return "Topology";
      case Category::BRANCH_LENGTHS:
	return "BranchLengths" ;
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
	return "v" ;
      case Category::FREQUENCIES :
	return "pi" ;
      case Category::SUBSTITUTION_RATES:
	return "r";
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
      case Category::AA_MODEL:
	return "AAPR"; 
      case Category::RATE_HETEROGENEITY: 
	return "SHAPEPR"; 
      default: 
	assert(0); 
      }
  } 


  std::vector<Category> getAllCategories()
  {
    return {  
      Category::TOPOLOGY, 
	Category::BRANCH_LENGTHS, 
	Category::FREQUENCIES, 
	Category::AA_MODEL,
	Category::SUBSTITUTION_RATES, 
	Category::RATE_HETEROGENEITY} ; 
  }


  Category getCategoryByPriorName(std::string name)
  {
    // bad but this is a bit exhausting... 
    auto cats = getAllCategories();
    
    std::transform(begin(name), end(name), begin(name), ::toupper);

    for(auto c : cats)
      {
	if(getPriorName(c).compare(name) == 0)
	  return c; 
      }

    tout << "Error in config file: did not find >" << name << "<" << std::endl; 

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
    else if(name.compare("brlens") == 0)
      return Category::BRANCH_LENGTHS; 
    else
      {
	assert(0); 
	return Category::TOPOLOGY; 
      }
  }



  std::unique_ptr<AbstractParameter> getParameterFromCategory(Category cat, nat id, nat idOfMyKind, std::vector<nat> partitions)
  {
    switch(cat)
      {
      case Category::TOPOLOGY :
	return  std::unique_ptr<AbstractParameter>( new TopologyParameter(id, idOfMyKind,partitions ));
      case Category::BRANCH_LENGTHS:
	return  std::unique_ptr<AbstractParameter>( new BranchLengthsParameter(id, idOfMyKind, partitions));
      case Category::FREQUENCIES :
	return  std::unique_ptr<AbstractParameter>( new FrequencyParameter(id, idOfMyKind,partitions));
      case Category::SUBSTITUTION_RATES :
	return  std::unique_ptr<AbstractParameter>( new RevMatParameter(id, idOfMyKind,partitions));
      case Category::RATE_HETEROGENEITY:
	return  std::unique_ptr<AbstractParameter>( new RateHetParameter(id, idOfMyKind,partitions));
      case Category::AA_MODEL :
	return std::unique_ptr<AbstractParameter>( new ProtModelParameter(id, idOfMyKind,partitions));
      default : 
	{
	  assert(0); 
	  return std::unique_ptr<AbstractParameter>( new RateHetParameter(id, idOfMyKind, partitions));
	}
      }    
  }  
}
