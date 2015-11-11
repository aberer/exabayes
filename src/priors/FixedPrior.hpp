#ifndef _FIXE_PRIOR
#define  _FIXE_PRIOR


class FixedPrior : public AbstractPrior
{
public: 
  FixedPrior(std::vector<double> fixedValues) : fixedValues(fixedValues) 
  {
  } 
  
  virtual double getLogProb(std::vector<double> values)  const
  {    
    assert(values.size() == fixedValues.size()); 
    for(nat i = 0; i < fixedValues.size() ; ++i)
      assert(fixedValues[i] == values[i]);
    return 0; 
  }

  virtual std::vector<double> drawFromPrior(Randomness &rand)  const
  {
    return fixedValues; 
  }

  virtual void print(std::ostream &out) const 
  {
    out << "Fixed(" ;     
    bool first = true; 
    for(auto v : fixedValues)
      {
	out << (first ? "" : ",") << v ; 
	if(first) first = false; 
      }
    out << ")"; 
  }

private: 
  std::vector<double> fixedValues; 
}; 

#endif
