#ifndef _LENGTH_PART_H
#define _LENGTH_PART_H

#include <vector>
#include <iostream>


template<typename TYPE> class Branch; 
typedef Branch<double> BranchLength; 
typedef Branch<std::vector<double>> BranchLengths; 
typedef Branch<void> BranchPlain; 

class TreeAln; 
class AbstractParameter; 


template<typename TYPE = void> class LengthPart
{
public: 
  LengthPart(){}

  // trait-like (could be better...)
  virtual ~LengthPart(){}
  
  friend std::ostream& operator<<(std::ostream &out, const LengthPart& rhs)
  {
    return out; 
  }

  void extractLength(const TreeAln &traln ){} 

  void lengthToString(std::ostream &out)  const
  {
  }

  bool lenEqual(const LengthPart &rhs) const 
  {
    return true; 
  }
}; 




template<> class LengthPart<double>
{
public: 
  LengthPart(double val = 0 )
    : length{val}  { }

  virtual ~LengthPart(){}
  /** 
      @brief gets the absolute (true) length of the branch
   */ 
  double getInterpretedLength( const AbstractParameter* param) const; 
  void setConvertedInternalLength(  const AbstractParameter* param, double length) ; 
  /**
     @brief sets the branch length (internal representation)
   */ 
  void setLength(double intLength){ length = intLength; }
  /**
     @brief gets the (internal) branch length
   */ 
  double getLength () const {return length; }

  bool lenEqual(const LengthPart<double>& rhs) const ; 

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
  LengthPart()
    : lengths{}    
  { } 

  virtual ~LengthPart(){}
  double getLength(const AbstractParameter* param) const ; 
  const std::vector<double>& getLengths() const {return lengths; }
  void setLengths(std::vector<double> _lengths) { lengths = _lengths; } 

  void extractLength(const TreeAln &traln, const BranchPlain &branch, const std::vector<AbstractParameter*> &params); 

  bool lenEqual(const LengthPart<std::vector<double> > &rhs) const ; 

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
