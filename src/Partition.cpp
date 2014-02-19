#include <algorithm>
#include <cstring>
#include <cassert>

#include "Partition.hpp"
#include "GlobalVariables.hpp"

#include "time.hpp"

#include "BitMask.hpp"

extern const char inverseMeaningDNA[16]; // TODO remove 

int Partition::maxCategories = 4; 	// TODO!!!

Partition::Partition(nat numTax, std::string name, int dataType, int states , int maxTipStates, bool saveMemory)
  : _numTax{numTax}
  , _name{name}
  , _wgtPtr{nullptr}
  , _y{nullptr}
  , _parsVect(nullptr)
  , _saveMemory{saveMemory}
  , _parsimonyInformative(0) 
{
  memset(&_partition, 0, sizeof(pInfo)); 
  defaultInit();
 
  // TODO the gap vector! 
  assert(not saveMemory); 

  _partition.dataType = dataType; 
  _partition.maxTipStates = maxTipStates; 
  _partition.states = states; 

  const partitionLengths 
    *pl = getPartitionLengths(&_partition); 
  
  // allocate the memory 
  _left.resize( pl->leftLength * (maxCategories + 1));  
  _right.resize( pl->rightLength * (maxCategories + 1)); 
  _EV.resize(pl->evLength);
  _tipVector.resize(pl->tipVectorLength);
  _EIGN.resize(pl->eignLength ); 
  _EI.resize(pl->eiLength);
  _substRates.resize(pl->substRatesLength); 
  _frequencies.resize(pl->frequenciesLength); 
  _empiricalFrequencies.resize(pl->frequenciesLength); 
  _globalScaler.resize(2 * _numTax); 
  _yPtrs.resize(_numTax + 1);
  
  _parsimonyScore.resize( 2 * numTax); 

  _gammaRates.resize(Partition::maxCategories); // =(

  _xVector.resize(_numTax, nullptr);
  _xSpaceVector.resize(_numTax); 
  
  setPtrs();
}


Partition::Partition(const Partition &rhs)
  : _numTax{rhs._numTax}
  , _partition(rhs._partition)
  , _name {rhs._name }
  , _wgtPtr {rhs._wgtPtr }
  , _y {rhs._y }
  , _left {rhs._left }
  , _right {rhs._right }
  , _EV {rhs._EV }
  , _tipVector {rhs._tipVector }
  , _EIGN {rhs._EIGN }
  , _EI {rhs._EI }
  , _substRates {rhs._substRates }
  , _frequencies {rhs._frequencies }
  , _empiricalFrequencies{rhs._empiricalFrequencies} 
  , _gammaRates  {rhs._gammaRates  }
  , _globalScaler {rhs._globalScaler }
  , _yPtrs {rhs._yPtrs }
  , _xVector {rhs._xVector }
  , _xSpaceVector {rhs._xSpaceVector }
  , _parsVect {rhs._parsVect }
  , _parsimonyScore {rhs._parsimonyScore }
  , _saveMemory {rhs._saveMemory }
  , _parsimonyInformative{rhs._parsimonyInformative}
{
  defaultInit();
  setPtrs();
}


void Partition::defaultInit()
{
  _partition.optimizeBaseFrequencies = PLL_FALSE; 
  _partition.partitionLH = 1; 
  _partition.executeModel = PLL_TRUE; 
}


void Partition::setPtrs()
{
  _partition.left = _left.data(); 
  _partition.right = _right.data();
  _partition.EV = _EV.data();
  _partition.tipVector = _tipVector.data(); // ***TODO*** actually needed? 
  _partition.EIGN = _EIGN.data();
  _partition.EI = _EI.data();
  _partition.substRates = _substRates.data();
  _partition.frequencies = _frequencies.data();
  _partition.empiricalFrequencies = _empiricalFrequencies.data();
  _partition.gammaRates = _gammaRates.data();
  _partition.globalScaler = _globalScaler.data(); 
  _partition.xVector = _xVector.data(); 
  _partition.xSpaceVector = _xSpaceVector.data();
  _partition.yVector = _yPtrs.data();
  _partition.parsVect = _parsVect.get();
  _partition.parsimonyScore = _parsimonyScore.data();		
}


// TODO should be movable, non-copyable 


Partition& Partition::operator=( Partition rhs)
{
  swap(rhs,*this); 
  return *this; 
} 


void swap(Partition& lhs, Partition& rhs)
{
  using std::swap; 

  swap(lhs._numTax, rhs._numTax); 
  swap(lhs._wgtPtr, rhs._wgtPtr); 
  swap(lhs._y, rhs._y); 
  swap(lhs._name, rhs._name); 
  swap(lhs._partition, rhs._partition); 
  swap(lhs._tipVector, rhs._tipVector); 
  swap(lhs._yPtrs, rhs._yPtrs); 
  swap(lhs._saveMemory, rhs._saveMemory); 
  swap(lhs._left, rhs._left);
  swap(lhs._right,rhs._right);
  swap(lhs._EV,rhs._EV);
  swap(lhs._EIGN,rhs._EIGN);
  swap(lhs._EV,rhs._EV);
  swap(lhs._EI,rhs._EI);
  swap(lhs._substRates,rhs._substRates);
  swap(lhs._frequencies,rhs._frequencies);
  swap(lhs._empiricalFrequencies, rhs._empiricalFrequencies); 
  swap(lhs._gammaRates,rhs._gammaRates);  
  swap(lhs._globalScaler, rhs._globalScaler); 
  swap(lhs._parsVect, rhs._parsVect); 
  swap(lhs._parsimonyScore, rhs._parsimonyScore); 
  swap(lhs._xVector, rhs._xVector); 
  swap(lhs._xSpaceVector, rhs._xSpaceVector); 
  swap(lhs._wgtPtr, rhs._wgtPtr); 
  swap(lhs._parsimonyInformative, rhs._parsimonyInformative); 
} 


void Partition::setAlignment( shared_pod_ptr<unsigned char> aln, int width)
{
  _y = aln; 
  _partition.width = width; 

  _yPtrs.at(0) = nullptr; 
  nat pos = 0; 
  for(nat i = 1 ; i < _numTax+1; ++i)
    {
      _yPtrs.at(i) = _y.get() + (pos * width) ; 
      ++pos; 
    }

  _partition.yVector = _yPtrs.data(); 

  prepareParsimony();
}
 

void Partition::setWeights(shared_pod_ptr<int> wgts, int width)
{
  _wgtPtr = wgts; 
  _partition.width = width; 
  _partition.wgt = _wgtPtr.get();
} 


void Partition::compressDNA(const std::vector<bool> &informative)
{
  auto totalNodes = 2 * _numTax;

  auto states = getStates(); 
  auto width = getWidth(); 

  auto entries = 0; 
  for(int i = 0; i < width; i++)    
    if(informative[i])
      entries += _partition.wgt[i];  

  nat compressedEntries = entries / PLL_PCF;

  if(entries % PLL_PCF != 0)
    compressedEntries++;

  nat compressedEntriesPadded = compressedEntries; 
#if (defined(__SSE3) || defined(__AVX))
  if(compressedEntries % INTS_PER_VECTOR != 0)
    compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
  else
    compressedEntriesPadded = compressedEntries;
#else
  compressedEntriesPadded = compressedEntries;
#endif     

  auto mask32 = BitMask<nat>();

  _partition.parsimonyLength = compressedEntriesPadded; 

  auto ptr = aligned_malloc<nat, size_t(EXA_ALIGN)>(compressedEntriesPadded * states * totalNodes); 
  std::fill(ptr, ptr  + compressedEntriesPadded * states * totalNodes, 0 );
  _parsVect = shared_pod_ptr<nat>(ptr);
  _partition.parsVect = _parsVect.get(); 

  const nat 
    *bitValue = getBitVector(_partition.dataType);

  for(nat taxIter = 1; taxIter <= _numTax; taxIter++)
    {			
      size_t
	compressedIndex = 0,
	compressedCounter = 0; 

      auto alnLine = _partition.yVector[taxIter]; 

      auto compressedValues = std::vector<nat>(states , 0 );
      auto compressedTips = std::vector<nat*>(states, nullptr);
      for(int k = 0 ; k < states; k++)
	compressedTips[k] = _partition.parsVect + ( (compressedEntriesPadded * states * taxIter) + (compressedEntriesPadded * k) ) ; 

      for(int index = 0; index < width; index++)
	{
	  if(informative[index])
	    {
	      auto 
		value = bitValue[alnLine[index]]; 
              
	      for(int w = 0; w < _partition.wgt[index]; w++)
		{      
		  for(int k = 0; k < states; k++)
		    {
		      if(value & mask32[k])
			compressedValues[k] |= mask32[compressedCounter];
		    }
                     
		  compressedCounter++;
                  
		  if(compressedCounter == PLL_PCF)
		    {
		      for(int k = 0; k < states; k++)
			{
			  compressedTips[k][compressedIndex] = compressedValues[k];
			  compressedValues[k] = 0;
			}                    
                          
		      compressedCounter = 0;
		      compressedIndex++;
		    }
		}
	    }
	}
                           
      for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	{   
	  for(;compressedCounter < PLL_PCF; compressedCounter++)              
	    for(int k = 0; k < states; k++)
	      compressedValues[k] |= mask32[compressedCounter];               
          
	  for(int k = 0; k < states; k++)
	    {
	      compressedTips[k][compressedIndex] = compressedValues[k];
	      compressedValues[k] = 0;
	    }                     
              
	  compressedCounter = 0;
	}           
    }               
}




bool Partition::isInformative(int site) const 
{
  auto dataType = getDataType();

  int
    check[256],   
    undetermined = getUndetermined(dataType);
  memset(check, 0, sizeof(int) * 256); 

  const nat
    *bitVector = getBitVector(dataType);
  
  for(nat i = 1 ; i <= _numTax; ++i)
    {      
      auto nucleotide = _partition.yVector[i][site];            
      ++check[nucleotide];
      assert(bitVector[nucleotide] > 0);                   
    }
  
  nat infoCtr = 0; 
  for( int i = 0; i < undetermined; ++i)
    if(check[i] > 0)
      ++infoCtr ; 

  if(infoCtr <= 1 )
    return false;    
  else
    {        
      for(int i = 0; i < undetermined; i++)
        {
          if(check[i] > 1)
            return true;
        } 
    }
     
  return false; 
}


std::vector<bool> Partition::determineUninformativeSites() const 
{
  auto result = std::vector<bool>(getWidth(), false);
  for(int i = 0; i <  getWidth() ; ++i)
    if(isInformative(i))
      result[i]= true; 
  return result; 
}


void Partition::prepareParsimony()
{
  if( not _parsimonyInformative.empty() )
    {
      auto tp = getTimePoint(); 
      compressDNA(_parsimonyInformative);
      tout << "compress: "<< getDuration(tp) << " sec"  << std::endl; 
    }
  else 
    {
      auto tp = getTimePoint(); 
      auto informative = determineUninformativeSites();
      tout << "determine uninfo: "<< getDuration(tp)  << " sec"<< std::endl; 

      tp = getTimePoint(); 
      compressDNA(informative);
      tout << "compress: "<< getDuration(tp) << " sec"  << std::endl; 
    }
}


void Partition::printAlignment()
{
  assert(0); 
  tout << "alignment " << std::endl; 
  for(nat i = 1; i < _numTax+1; ++i)
    {
      for(int j = 0; j < getWidth(); ++j)
	tout << inverseMeaningDNA[_partition.yVector[i][j]] ; 
      tout << "\n"; 
    }
  tout << std::endl; 
}
