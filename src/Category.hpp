#ifndef _CATEGORIES_H
#define _CATEGORIES_H

#include <vector>
#include <string>
#include <memory>
#include "parameters/AbstractParameter.hpp"

enum class Category
{  
  TOPOLOGY = 0, 
  BRANCH_LENGTHS = 1, 
  FREQUENCIES = 2,
  SUBSTITUTION_RATES = 3,
  RATE_HETEROGENEITY = 4,	
  AA_MODEL= 5  

} ; 


namespace CategoryFuns 
{
  std::vector<Category> getAllCategories(); 
  std::string getLongName(Category cat); 
  std::string getShortName(Category cat); 
  std::string getPriorName(Category cat); 
  Category getCategoryByPriorName(std::string name); 
  Category getCategoryFromLinkLabel(std::string name); 

  std::shared_ptr<AbstractParameter> getParameterFromCategory(Category cat, nat id ); 
} 

#endif
 
