#ifndef _PARAMETER_FILE
#define _PARAMETER_FILE

#include <sstream>
#include <string>
#include <iostream> 
#include <limits>

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

  ParameterFile(std::string workdir, std::string runname, nat runid, nat couplingId)
  {
    std::stringstream ss ; 
    // TODO portability 
    ss << workdir << (workdir.compare("") == 0 ? "" : "/")  << "ExaBayes_parameters." << runname << "." << runid ; 
    if(couplingId != 0 )
      ss << ".hot-"<<  couplingId; 
    fullFilename = ss.str();

    std::ofstream fh(fullFilename); 
    fh << "" ; 
    fh.close(); 
  }


  void initialize(const TreeAln& traln, std::vector<AbstractParameter*> parameters,  nat someId ) const 
  {    
    std::ofstream fh(fullFilename,std::fstream::app|std::fstream::out); 
    fh << "[ID: " << someId << "]" << std::endl; 

    fh << "Gen\t";
    fh << "LnPr\t"; 
    fh << "LnL\t" ; 
    fh << "TL\t" ; 

    bool isFirst = true; 
    for(auto &p : parameters)
      {
	if(p->isPrintToParamFile())
	  {
	    if(isFirst) 
	      isFirst = false; 
	    else 
	      fh << "\t" ; 
	    p->printAllComponentNames(fh, traln); 
	  }
      }

    fh << std::endl; 

    fh.close(); 
  }


  void sample(const TreeAln &traln, const std::vector<AbstractParameter*> parameters, nat gen, double lnPr) const 
  {
    std::ofstream fh(fullFilename, std::fstream::app|std::fstream::out); 
    
    fh << gen << "\t"; 
    fh << std::setprecision(std::numeric_limits<double>::digits10)  << std::scientific; 
    fh << lnPr << "\t"; 
    fh << traln.getTr()->likelihood << "\t" ; 
    fh << Branch(0,0,traln.getTreeLengthExpensive()).getInterpretedLength(traln) << "\t"; 

    bool isFirst = true; 
    for(auto &p : parameters)
      {
	if(p->isPrintToParamFile())
	  {
	    if(isFirst)
	      isFirst = false; 
	    else 
	      fh << "\t"; 
	    p->printSample(fh,traln);
	  }
      }
    fh << std::endl; 
    
    fh.close(); 
  }

  void finalize() const 
  {
  }
  
private: 
  std::string fullFilename; 
}; 


#endif
