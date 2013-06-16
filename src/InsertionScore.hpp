#ifndef _INSERTION_SCORE 
#define _INSERTION_SCORE 


class  InsertionScore
{
public: 
  InsertionScore(branch _b, vector<nat> _tmp) : b(_b), partitionParsimony(_tmp){}  
  branch getBranch() const  {return b; }

  Branch getNewBranch() const {Branch br; br.initFromLegacy(b); return br; }

  double getWeight() const {return  logProb; }
  void setWeight(double w) { logProb = w; }

  nat getScore() const
  {
    nat result = 0; 
    for(auto b : partitionParsimony)
      result += b; 
    return result; 
  }

  nat getPartitionScore(int model) const{return partitionParsimony[model] ; }


private: 
  branch b; 
  vector<nat> partitionParsimony; 
  double logProb;

  friend ostream& operator<< (ostream &out, const InsertionScore &rhs) { 
    out <<  "(" << rhs.b.thisNode << "," << rhs.b.thatNode << ")" << "=" ; 
    for(auto elem : rhs.partitionParsimony)
      out << elem << "," ; 
    return out; 
  }
} ; 


#endif
