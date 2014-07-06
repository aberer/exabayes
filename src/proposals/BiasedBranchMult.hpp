#ifndef BIASEDBRANCHMULT_H
#define BIASEDBRANCHMULT_H



class BiasedBranchMult : public BranchLengthMultiplier
{
public: 
  BiasedBranchMult(double m)
    : BranchLengthMultiplier(m)
  {
    this->_relativeWeight = 1.;
    this->_name = "biasBLMult"; 
  }
  
  BiasedBranchMult(const BiasedBranchMult &rhs)  = default; 
  BiasedBranchMult(BiasedBranchMult&& rhs) = default; 
  BiasedBranchMult& operator=(const BiasedBranchMult &rhs)  = default; 
  BiasedBranchMult& operator=( BiasedBranchMult &&rhs)  = default; 
 
  virtual AbstractProposal* clone() const { return new BiasedBranchMult(*this); }


  BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const 
  {
    auto result = BranchPlain{}; 
    auto myEps = 1e-6; 
    
    auto param = getPrimaryParameterView()[0]; 
    
    auto b2len = std::unordered_map<BranchPlain,double> {}; 
    
    auto sum = 0.; 
    for(auto bl : traln.extractBranches(param))
      {
	auto len = bl.getInterpretedLength( param );
	auto prob = std::exp(- std::log(len ) * 0.25 )  + myEps; 
	sum += prob; 
	b2len.insert(std::make_pair(bl.toPlain(), prob)); 
      }
    
    for(auto elem :  b2len)
      b2len[elem.first] /= sum; 
    
    auto r =  rand.drawRandDouble01();
    for(auto elem : b2len)
      {
	r -= elem.second; 
	if(r <= 0 )
	  {
	    result = elem.first; 
	    break; 
	  }
      }

    return result; 
  }
}; 

#endif /* BIASEDBRANCHMULT_H */
