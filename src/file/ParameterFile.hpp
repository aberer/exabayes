#ifndef _PARAMETER_FILE
#define _PARAMETER_FILE


#include <sstream>
#include <string>
#include <iostream> 
#include <limits>
#include <vector>

#include "parameters/AbstractParameter.hpp"
#include "Branch.hpp"
#include "OutputFile.hpp"

/**
   @notice opening and closing the stream all the time may appear
   weird. Since properly the copy/move constructor for streams is not
   implemented properly in g++, it appears as a good choice to do it
   as implemented below.
 */

class  ParameterFile  : public OutputFile
{
public: 
  ParameterFile(std::string workdir, std::string runname, nat runid, nat couplingId); 
  void initialize(const TreeAln& traln, std::vector<AbstractParameter*> parameters,  nat someId )  ; 
  void sample(const TreeAln &traln, const std::vector<AbstractParameter*> parameters, nat gen, double lnPr)  ; 

  void finalize() const  { }	// NO IMPLEMENT  
  void regenerate(std::string workdir, std::string prevId, nat gen) ; 
  
private: 
  nat runid; 
  nat couplingId; 
}; 


#endif
