#ifndef UNSUPPORTEDTREEACTION_H
#define UNSUPPORTEDTREEACTION_H

#include <exception>

class UnsupportedTreeAction : public std::exception
{
public:
  UnsupportedTreeAction()
  {
  }
  
  virtual ~UnsupportedTreeAction()
  {
  }
  
};




#endif /* UNSUPPORTEDTREEACTION_H */
