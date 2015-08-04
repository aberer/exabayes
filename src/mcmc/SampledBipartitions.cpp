#include "SampledBipartitions.hpp"
#include "Topology.hpp"
#include "TreePrinter.hpp"
#include "NewickParser.hpp"
#include "ParallelSetup.hpp"

using std::make_pair; 

void SampledBipartitions::deserialize( std::istream &in )
{
  auto num = cRead<size_t>(in);
  for(auto i = 0u; i < num; ++i )
    {
      auto bipKey = bitvector();
      auto bipVal =  bitvector(); 
      bipKey.deserialize(in);
      bipVal.deserialize(in);
      _bip2trees.emplace(bipKey, bipVal);
    }
}

void SampledBipartitions::serialize( std::ostream &out) const
{
  cWrite<size_t>(out,_bip2trees.size());
  for(auto const &v : _bip2trees)
    {
      v.first.serialize(out);
      v.second.serialize(out); 
    }
}


void SampledBipartitions::addTree(Topology const &t, uint64_t trIndex)
{
  for(auto tIter = t.begin(); tIter != t.end(); ++tIter)
    {
      if( tIter->isOuterBranch() )
        continue; 
      
      auto bip = tIter.getBipOrDummy();
      
      auto bipIter = _bip2trees.find(bip);
      if( bipIter == _bip2trees.end())
        bipIter = _bip2trees.find( ~ bip );

      // either add it or update the mapping
      if( bipIter == _bip2trees.end() )
        {
          auto bipVal = bitvector();
          bipVal.set( trIndex);
          _bip2trees.emplace(bip, bipVal);
        }
      else
        {
          bipIter->second.set( trIndex ); 
        }
    }
}


void SampledBipartitions::addTree(TreeAln const& traln, uint64_t trIndex)
{
  auto tp = TreePrinter(false, false, false);
  auto str = tp.printTree(traln);

  auto nw = NewickParser<Topology>(str);
  nw.initalize(); 
  auto t = nw.getTree();

  addTree(t, trIndex);
}


void SampledBipartitions::reduce( SampledBipartitions const& others)
{
  for(auto &v : others._bip2trees)
    {
      auto &bip = v.first;
      auto &occ = v.second;

      auto iter = _bip2trees.find( bip) ;
      if( iter == _bip2trees.end())
        iter = _bip2trees.find( ~ bip ) ;

      if( iter == _bip2trees.end() ) // not there yet
        {
          _bip2trees.insert(make_pair(bip, occ)); 
        }
      else
        {
          iter->second |=  occ; 
        }
    }
}


std::ostream& operator<<(std::ostream& s,  SampledBipartitions const& rhs)
{
  for(auto const& v : rhs._bip2trees)
    s << v.first << "\t->\t" << v.second << std::endl; 
  
  return s; 
}

