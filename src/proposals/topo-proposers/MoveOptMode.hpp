#ifndef MOVEOPTMODE_H
#define MOVEOPTMODE_H

#include <iosfwd>

///////////////////////////////////////////////////////////////////////////////
//                               MOVE OPT MODE                               //
///////////////////////////////////////////////////////////////////////////////
enum class MoveOptMode : int
{
    NONE = 0,
    ONLY_SWITCHING = 1,
    ALL_INTERNAL = 2,
    ALL_IN_MOVE  = 3,
    ALL_SURROUNDING = 4
};


// ____________________________________________________________________________
std::ostream&                   operator<<(
    std::ostream&      s,
    const MoveOptMode& c);

#endif /* MOVEOPTMODE_H */

