#ifndef _LENGTH_PART_H
#define _LENGTH_PART_H

#include <iostream>

#include "BranchFwd.hpp"



class TreeAln; 
class AbstractParameter; 


template<typename TYPE = void> class LengthPart
{
public: 
  friend std::ostream& operator<<(std::ostream &out, const LengthPart& rhs)
  {
    return out; 
  }

  // this could be a solution to the parameter dilemma. 
  // but let's do more hard-coding for now 
  static typename std::conditional< std::is_same<TYPE,double>::value, AbstractParameter* , std::vector<AbstractParameter*> >::type PARAM; 

  void extractLength(const TreeAln &traln ){} 

  void lengthToString(std::ostream &out)  const
  {
  }

}; 


template<> class LengthPart<double>
{
public: 
  /** 
      @brief gets the absolute (true) length of the branch
   */ 
  double getInterpretedLength(const TreeAln &traln, const AbstractParameter* param) const; 
  void setConvertedInternalLength(const TreeAln& traln,  const AbstractParameter* param, double length) ; 
  /**
     @brief sets the branch length (internal representation)
   */ 
  void setLength(double intLength){ length = intLength; }
  /**
     @brief gets the (internal) branch length
   */ 
  double getLength () const {return length; }


  friend std::ostream& operator<<(std::ostream &out, const LengthPart &rhs)
  {
    return  out << ":" << rhs.length; 
  }

  void lengthToString(std::ostream &out) const 
  {
    out << ":" << length ; 
  }

  void extractLength(const TreeAln &traln, const BranchPlain& branch, const AbstractParameter*  param); 

protected: 
  double length; 
}; 


template<> class LengthPart<std::vector<double>>
{
 public: 
  double getLength(const AbstractParameter* param) const ; 
  const std::vector<double>& getLengths() const {return lengths; }
  void setLengths(std::vector<double> _lengths) { lengths = _lengths; } 

  void extractLength(const TreeAln &traln, const BranchPlain &branch, const std::vector<AbstractParameter*> &params); 

  friend std::ostream& operator<<(std::ostream &out, const LengthPart &rhs)
  {
    out << ":"; 
    bool isFirst = true; 
    for(auto &v : rhs.lengths)
      {
	if(isFirst)
	  isFirst = false; 
	else 
	  out << ","; 
	out << v ; 
      }
    return out; 
  }

  
  void lengthToString(std::ostream &out) const 
  {
    out << ":[" ; 
    bool isFirst = false; 
    for(auto &v : lengths)
      {
	if(isFirst)
	  isFirst = false; 
	else 
	  out << "," ; 
	out << v ; 
      }
    out << "]"; 
  }

 protected: 
  std::vector<double> lengths;
}; 



#endif
