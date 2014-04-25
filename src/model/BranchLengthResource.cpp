#include "BranchLengthResource.hpp"

#include "model/TreeAln.hpp"


void BranchLengthResource::initialize(nat numTax, nat numPart )
{
  _numTax = numTax; 
  _numPart = numPart; 

  // meh meh meh 
  // PLL_NUM_BRANCHES = _numPart;

  _zqr.resize(_numPart); 
  _currentZQR.resize(_numPart); 
  _currentLZR.resize(_numPart); 
  _currentLZQ.resize(_numPart); 
  _currentLZS.resize(_numPart); 
  _currentLZI.resize(_numPart); 
  _lzs.resize(_numPart); 
  _lzq.resize(_numPart); 
  _lzr.resize(_numPart); 
  _lzi.resize(_numPart); 
  _coreLZ.resize(_numPart); 
  _curvatOK.resize(_numPart); 
  _partitionSmoothed.resize(_numPart); 
  _partitionConverged.resize(_numPart); 
  
  for(nat i = 0; i < _numTax; ++i)
    {
      _qz.push_back(std::vector<double>(_numPart,0));
      _rz.push_back(std::vector<double>(_numPart,0));
    }
  
  for(nat i = 0; i < _numTax + 3 * (_numTax-1); ++i)
    _z.push_back(std::vector<double>(_numPart,0));
} 


void BranchLengthResource::assign(TreeAln &traln) 
{
  auto &tr = traln.getTrHandle();

  tr.zqr = _zqr.data();
  tr.currentZQR = _currentZQR.data(); 
  tr.currentLZR = _currentLZR.data(); 
  tr.currentLZQ = _currentLZQ.data(); 
  tr.currentLZS = _currentLZS.data(); 
  tr.currentLZI = _currentLZI.data(); 
  tr.lzs = _lzs.data(); 
  tr.lzq = _lzq.data(); 
  tr.lzr = _lzr.data(); 
  tr.lzi = _lzi.data(); 
  tr.coreLZ = _coreLZ.data(); 
  tr.curvatOK = _curvatOK.data(); 

  tr.partitionSmoothed = _partitionSmoothed.data(); 
  tr.partitionConverged = _partitionConverged.data(); 

  for(nat i = 0; i < _numTax ; ++i)
    {
      tr.td[0].ti[i].qz = _qz.at(i).data(); 
      tr.td[0].ti[i].rz = _rz.at(i).data(); 
    }

  nat ctr = 0; 
  for(auto &n : traln._nodes)
    {
      n.z = _z.at(ctr).data(); 
      ++ctr; 
    }
  assert(ctr == traln._nodes.size()); 
}  
