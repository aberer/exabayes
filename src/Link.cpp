#include "Link.hpp"

#include <fstream>

std::ostream& operator<<(std::ostream& s, const Link& c)
{
  s << c._priNode << "," << c._secNode; 
  return s;
}
