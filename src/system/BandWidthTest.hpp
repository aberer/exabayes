#ifndef _BAND_WIDTH_TEST
#define _BAND_WIDTH_TEST

#include "common.h"
#include "Communicator.hpp"

///////////////////////////////////////////////////////////////////////////////
//                              BAND WIDTH TEST                              //
///////////////////////////////////////////////////////////////////////////////
class BandWidthTest
{
public:
    // ________________________________________________________________________
    nat                    determineOptimum(
        Communicator& comm,
        std::string   fileName);
};


#endif
