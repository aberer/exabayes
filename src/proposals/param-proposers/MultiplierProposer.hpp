#ifndef MULTIPLIERPROPOSAL_H
#define MULTIPLIERPROPOSAL_H

#include "AbstractProposer.hpp"

///////////////////////////////////////////////////////////////////////////////
//                            MULTIPLIER PROPOSER                            //
///////////////////////////////////////////////////////////////////////////////
class MultiplierProposer : public AbstractProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    MultiplierProposer(
        double minVal,
        double maxVal);
    // ________________________________________________________________________
    MultiplierProposer(
        const MultiplierProposer& rhs)
        : AbstractProposer(rhs){}
    // ________________________________________________________________________
    virtual AbstractProposer*                      clone() const
    {return new MultiplierProposer(*this);}
    // ________________________________________________________________________
    virtual std::vector<double>                    proposeValues(
        std::vector<double>oldValues,
        double             parameter,
        Randomness&        rand,
        log_double&        hastings);
};


#endif
