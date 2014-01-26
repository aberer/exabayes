#include "extensions.hpp"


void formatRange(std::ostream &out, const std::vector<nat> &values)  
{
  bool inRange = false; 
  for(auto iter = begin(values); iter != end(values) ; ++iter )
    {
      if(iter == begin(values))
	out << *iter; 
      else if(iter == end(values) -1 )
	{
	  if(inRange)
	    out << "-" << *iter ;  
	  else 
	    out << "," << *iter; 
	}
      else 
	{ 
	  bool haveBeenInRange = ((iter-1) ==  begin(values)) || inRange; 
	  inRange = *iter == *(iter-1) + 1;  

	  if(not haveBeenInRange )
	    {
	      out << "," << *iter ; 
	    }
	  else
	    {
	      if(not inRange)
		out << "-" << *iter; 
	    }
	}
    }
  
}
