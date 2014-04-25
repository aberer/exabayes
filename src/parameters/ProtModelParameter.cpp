#include "ProtModelParameter.hpp" 
#include "model/ProtModel.hpp"


void ProtModelParameter::applyParameter(TreeAln& traln,  const ParameterContent &content) const 
{
  for(auto &m : _partitions)
    {
      assert(content.protModel.size() == 1); 
      traln.setProteinModel(m, content.protModel[0]);
    }
}
 
ParameterContent ProtModelParameter::extractParameter(const TreeAln &traln)  const  
{

  auto& partition = traln.getPartition(_partitions.at(0)); 
  auto model = ProtModel(partition.getProtModels()); 
  auto result =   ParameterContent{} ; 
  result.protModel = {model}; 
  return result; 
}
   
void ProtModelParameter::printSample(std::ostream& fileHandle, const TreeAln &traln ) const 
{
  auto& partition = traln.getPartition(_partitions.at(0)); 
  auto name = ProtModelFun::getName(ProtModel(partition.getProtModels() )); 
  fileHandle << name; 
}
 
void ProtModelParameter::printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  
{
  fileHandle << "aaModel{" ; 
  bool isFirst = true; 
  for(auto &p : _partitions)
    {
      fileHandle << (isFirst ? "" : ",") << p ; 
      isFirst = false; 
    }
  fileHandle << "}"; 
}

 
void ProtModelParameter::verifyContent(const TreeAln &traln, const ParameterContent &content) const  
{
  if(content.protModel.size() != 1)
    {
      tout << "incorrect number of models in parameter content. This is a programming error." << std::endl; 
      assert(0); 
    }
} 


void ProtModelParameter::checkSanityPartitionsAndPrior(const TreeAln &traln) const 
{
  checkSanityPartitionsAndPrior_FreqRevMat(traln);
  
  if(traln.getPartition(_partitions.at(0)).getStates() != 20)
    {
      std::cerr << "Error: in the config file you specified partition " << _partitions.at(0) << " to have an amino acid model. However, previously this partition was declared to have a different data type." << std::endl; 
      exitFunction(-1, true); 
    }
}
