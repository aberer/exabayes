#ifndef REV_MAT_PARAMETER
#define REV_MAT_PARAMETER


#include "Category.hpp"
#include "AbstractParameter.hpp"
  
class RevMatParameter : public AbstractParameter
{
public: 
  RevMatParameter(nat id// , nat numComponents
		  ) 
    : AbstractParameter(Category::SUBSTITUTION_RATES, id// , numComponents
			)
  {
    // modifiesBL = true; 
  }

  virtual void applyParameter(TreeAln& traln, const ParameterContent &content) const; 
  virtual ParameterContent extractParameter(const TreeAln &traln )  const;   
  virtual AbstractParameter* clone () const {return new RevMatParameter(*this); } 

  virtual void printSample(std::ostream& fileHandle, const TreeAln &traln) const ; 
  virtual void printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const ; 
}; 

#endif
