#include "LikeArrayManager.hpp"
#include <cassert>


void LikeArrayManager::flushBackup( std::unordered_map<bitvector,CondLikeSubtree> & map  )
{
  for(auto iter = _backupArrays.begin(); iter != _backupArrays.end(); ++iter)
    {
      iter->second.releaseArray(*_arrayReservoir);
      _backupArrays.erase(iter);
    }
}



static inline void moveFromMapToMap( std::unordered_map<bitvector,CondLikeSubtree>::iterator &it,
                                     std::unordered_map<bitvector,CondLikeSubtree> &toBeRemoved,
                                     std::unordered_map<bitvector,CondLikeSubtree> &toBeInserted )
{
  auto key = it->first;
  auto val = it->second;
  toBeRemoved.erase(it);
  toBeInserted.emplace( key, val);


  //that was the original idea ... 
  
  // auto elem = it; 
  // auto&& elem = std::move(*it);
  // toBeRemoved.erase(it);
  // toBeInserted.emplace(elem); 
}


void LikeArrayManager::imprint()
{
  flushBackup(_backupArrays);
}


void LikeArrayManager::makeVanish(bitvector bv)
{
  auto iter = _curArrays.find(bv);
  if(  iter != _curArrays.end() ) 
    {
      // computed in previous generation, must create backup 
      moveFromMapToMap( iter, _curArrays, _backupArrays);
      return; 
    }

  iter = _newArrays.find(bv);
  if( iter != _newArrays.end() )
    {
      moveFromMapToMap( iter, _newArrays, _unusedArrays); 
      return; 
    }
}


void LikeArrayManager::makeExist(bitvector bv)
{
  auto iter = _unusedArrays.find(bv);
  if( iter != _unusedArrays.end() ) 
    {
      moveFromMapToMap( iter, _unusedArrays, _newArrays);
      return; 
    }

  // that does not quite work out ; not likely to occur anyway
  
  // iter = _backupArrays.find(bv);
  // if(iter != _backupArrays.end() )
  //   {
  //     moveFromMapToMap(iter, _backupArrays, _curArrays);
  //     return; 
  //   }
  
  // _newArrays
  _newArrays.emplace(bv,bv); 
}



std::ostream& operator<<(std::ostream& s,  LikeArrayManager const& c)
{
  s << "BACKUP " << std::endl;
  for(auto & v : c._backupArrays)
    s << v.first << std::endl;
  
  s << "CUR" << std::endl;
  for(auto &v : c._curArrays)
    s << v.first << std::endl; 

  s << "NEW"<< std::endl;
  for(auto &v : c._newArrays)
    s << v.first << std::endl;

  s << "UNUSED" << std::endl;
  for(auto &v : c._unusedArrays)
    s << v.first << std::endl;
  
  return s;
}
