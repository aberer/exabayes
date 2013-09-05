#ifndef _BRANCH_FWD_HPP
#define _BRANCH_FWD_HPP

/**
   @file BranchFwd.hpp
   @brief forward declarations for Branches. 

   These are needed sometimes because of the mutual inclusion of
   TreeAln and Branch. But it is challenging.
 */


#include <vector>

template<typename TYPE>
class Branch; 
typedef Branch<double> BranchLength; 
typedef Branch<std::vector<double>> BranchLengths; 
typedef Branch<void> BranchPlain; 

#endif
