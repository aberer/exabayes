#ifndef _PARAMETER_FILE
#define _PARAMETER_FILE

#include <sstream>
#include <string>
#include <iostream> 
#include <limits>
#include <vector>

#include "parameters/AbstractParameter.hpp"
#include "Branch.hpp"


/**
   @notice opening and closing the stream all the time may appear
   weird. Since properly the copy/move constructor for streams is not
   implemented properly in g++, it appears as a good choice to do it
   as implemented below.
 */

class  ParameterFile 
{
public: 
  ParameterFile(std::string workdir, std::string runname, nat runid, nat couplingId); 
  void initialize(const TreeAln& traln, std::vector<AbstractParameter*> parameters,  nat someId ) const ; 
  void sample(const TreeAln &traln, const std::vector<AbstractParameter*> parameters, nat gen, double lnPr) const ; 

  void finalize() const  { }	// NO IMPLEMENT
  
  void regenerate(std::string prevId, nat gen) ; 
  
private: 
  std::string fullFilename; 
  nat runid; 
  nat couplingId; 
}; 


#endif
