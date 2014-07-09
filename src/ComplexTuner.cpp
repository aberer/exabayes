#include "ComplexTuner.hpp"
#include <cmath>
#include "system/GlobalVariables.hpp"


double ComplexTuner::tuneParameter(int batch, double parameter, bool increase ) const 
{
  auto delta = 1. / (sqrt(   (batch + 1))); 

  if(delta > _maxDelta )
    delta = _maxDelta; 

  double val = _logScale ? log(parameter) : parameter; 
  if(increase)
    val += delta;
  else 
    val -= delta; 
  
  if(_logScale)
    val = exp(val); 

  val = fmax(_minBound, val); 
  val = fmin(_maxBound, val); 

  return val; 
} 

void ComplexTuner::tune() 
{
  auto curRatio = _sctr.getRatioInLastInterval();
  auto batch = _sctr.getBatch();

  if( curRatio < _prevSuccess )	
    _tuneUp = not _tuneUp; 
  _prevSuccess = curRatio; 

  auto newParam = tuneParameter(batch, _parameter, _tuneUp);

  _parameter = newParam; 
  _sctr.nextBatch();
}  
