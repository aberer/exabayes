#ifndef PROT_MODEL_HPP_
#define PROT_MODEL_HPP_

#include <vector>
#include <iostream>
#include <string>
#include <cassert>

#include "axml.h"

enum class ProtModel : int 
{
  DAYHOFF_T =    DAYHOFF,
    DCMUT_T =      DCMUT,
    JTT_T =        JTT,
    MTREV_T =      MTREV,
    WAG_T =        WAG,
    RTREV_T =      RTREV,
    CPREV_T =      CPREV,
    VT_T =         VT,
    BLOSUM62_T =   BLOSUM62,
    MTMAM_T =      MTMAM,
    LG_T =         LG,
    MTART_T =      MTART,
    MTZOA_T =      MTZOA,
    PMB_T =        PMB,
    HIVB_T =       HIVB,
    HIVW_T =       HIVW,
    JTTDCMUT_T =   JTTDCMUT,
    FLU_T =        FLU
}; 


namespace ProtModelFun
{
  std::string getName(ProtModel mod); 
  std::vector<ProtModel> getAllModels();
  std::tuple<bool,ProtModel> getModelFromStringIfPossible(const std::string & modelString); 
}

std::ostream& operator<<(std::ostream &out, const ProtModel &rhs); 

namespace std
{
  template<> struct less<ProtModel> {
    bool operator()(const ProtModel& a, const ProtModel& b)  const 
    {
      return int(a) < int(b);
    } 
  }; 


  template<> struct hash<ProtModel>
  {
    unsigned int operator() (const ProtModel& rhs) const 
    {
      return std::hash<int>()(int(rhs));
    }
  }; 
}


#endif
