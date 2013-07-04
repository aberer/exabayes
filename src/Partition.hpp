#ifndef _PARTITION_H
#define _PARTITION_H

#include <vector>

#include "TreeAln.hpp"


class Partition 
{
public: 
  Partition(int model, const TreeAln &traln )
    : alpha(traln.getAlpha(model))      
  {
    pInfo *partition = traln.getPartition(model);
    for(int i = 0; i < partition->states; ++i)
      stateFreqs.push_back(partition->frequencies[i]);     

    for(nat i = 0; i < numStateToNumInTriangleMatrix(partition->states); ++i)
      revMat.push_back(partition->substRates[i]); 
  }

  void setAlpha(double _alpha){ alpha = _alpha; }
  void setRevMat(vector<double> _revMat ) {revMat = _revMat; }
  void setStateFreqs(vector<double> _stateFreqs) {stateFreqs = _stateFreqs; }

  double getAlpha() const {return alpha; }
  vector<double> getRevMat() const {return revMat; }
  vector<double> getStateFreqs() const {return stateFreqs; }

private: 
  double alpha; 
  vector<double> revMat; 
  vector<double> stateFreqs;   

}; 



#endif
