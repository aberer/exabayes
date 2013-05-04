#ifndef _PATH_H
#define _PATH_H

#include <list>

#include "axml.h"
#include "bayes.h"
#include "branch.h"

class Path
{
public: 
  Path();
  ~Path();   
  void pushStack(branch value); 


private: 
  list<branch> stack; 

  
  

  
}; 


#endif
