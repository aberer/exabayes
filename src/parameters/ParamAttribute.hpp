#ifndef _PARAM_ATTRIBUTE
#define _PARAM_ATTRIBUTE

#include "ComplexTuner.hpp"

// another one of these generic container classes

// not very much happy with it...

///////////////////////////////////////////////////////////////////////////////
//                              PARAM ATTRIBUTE                              //
///////////////////////////////////////////////////////////////////////////////
class ParamAttribute
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
public:
    ComplexTuner _convTuner;
    ComplexTuner _nonConvTuner;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    ParamAttribute()
        :  _convTuner{0., 0.01, 10, 0.1, false}
        , _nonConvTuner{0., 0.1, 10, 0.1, false}
    {}

    ParamAttribute(
        ComplexTuner c,
        ComplexTuner n)
        : _convTuner{c}
        , _nonConvTuner{n}
    {}
};


#endif
