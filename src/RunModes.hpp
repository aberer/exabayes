#ifndef _RUNMODES_HPP
#define _RUNMODES_HPP

/** 
    @brief various modes relevant to likelihood computation.  
    Mostly, the TreeAln shouldo be aware of this  
 */ 

enum class RunModes : int
{
NOTHING = 0, 
  PARTITION_DISTRIBUTION = 1, 
  MEMORY_SEV = 2
}; 

#endif
