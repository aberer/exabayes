#ifndef SLIDINGPROPOSAL
#define SLIDINGPROPOSAL

#include "AbstractProposer.hpp"

///////////////////////////////////////////////////////////////////////////////
//                              SLIDING PROPOSER                             //
///////////////////////////////////////////////////////////////////////////////
class SlidingProposer : public AbstractProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    SlidingProposer(
        double minVal,
        double maxVal,
        bool   minMaxIsRelative);
    // ________________________________________________________________________
    virtual ~SlidingProposer(){}
    // ________________________________________________________________________
    SlidingProposer(
        const SlidingProposer&rhs)
        : AbstractProposer(rhs)
        , minMaxIsRelative(rhs.minMaxIsRelative)
    {}
    // ________________________________________________________________________
    double                                         proposeOneValue(
        double     oldVal,
        double     parameter,
        Randomness&rand,
        log_double&hastings);
    // ________________________________________________________________________
    std::vector<double>                            proposeRelativeMany(
        std::vector<double>oldValues,
        double             parameter,
        Randomness&        rand,
        log_double&        hastings);
    // ________________________________________________________________________
    virtual std::vector<double>                    proposeValues(
        std::vector<double>oldValues,
        double             parameter,
        Randomness&        rand,
        log_double&        hastings);
    // ________________________________________________________________________
    virtual AbstractProposer*                      clone() const
    {return new SlidingProposer(*this);}

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    bool minMaxIsRelative;
};


#endif
