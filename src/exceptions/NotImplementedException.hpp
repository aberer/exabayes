#ifndef NOTIMPLEMENTEDEXCEPTION_H
#define NOTIMPLEMENTEDEXCEPTION_H

#include <exception>

class NotImplementedException : public std::exception
{
public:
  NotImplementedException()
  {
  }
  
  virtual ~NotImplementedException()
  {
  }
};




#endif /* NOTIMPLEMENTEDEXCEPTION_H */
