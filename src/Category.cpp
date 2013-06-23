#include "Category.hpp"
#include <cassert>

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
