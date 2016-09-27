#ifndef RATES_LIDINGP_ROPOSER_HPP
#define RATES_LIDINGP_ROPOSER_HPP

#include "AbstractProposer.hpp"


///////////////////////////////////////////////////////////////////////////////
//                           RATE SLIDING PROPOSER                           //
///////////////////////////////////////////////////////////////////////////////
class RateSlidingProposer : public AbstractProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    RateSlidingProposer(
        double minValI,
        double maxValI);
    // ________________________________________________________________________
    RateSlidingProposer(
        const RateSlidingProposer&rhs);
    // ________________________________________________________________________
    virtual std::vector<double>                    proposeValues(
        std::vector<double>oldValues,
        double             parameter,
        Randomness&        rand,
        log_double&        hastings);
    // ________________________________________________________________________
    virtual AbstractProposer*                      clone() const
    {return new RateSlidingProposer(*this);}
};


#endif
