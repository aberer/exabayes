#define _TOPOLOGY_PARAMETER
#ifndef _TOPOLOGY_PARAMETER


class TopologyParameter : public AbstractParameter
{
public: 
  virtual void applyParameter(TreeAln& traln, vector<nat> model, ParameterContent &content) const; 
  virtual ParameterContent& extractParameter(const TreeAln &traln, vector<nat> model)  const;   
  
}; 

#endif
