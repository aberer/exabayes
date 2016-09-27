#ifndef DIRICHLETPROPOSAL_H
#define DIRICHLETPROPOSAL_H

#include "AbstractProposer.hpp"
#include "RateHelper.hpp"

///////////////////////////////////////////////////////////////////////////////
//                             DIRICHLET PROPOSER                            //
///////////////////////////////////////////////////////////////////////////////
class DirichletProposer : public AbstractProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    /**
     *  @brief constructs a dirichlet proposal
     *
     *  @param minMaxIsRelative
     *  indicates whether the previous two boundary arguments are
     *  relative to a value (e.g., revmat) or whether they sum up to 1.
     */
    DirichletProposer(
        double minVal,
        double maxVal,
        bool   minMaxIsRelative);
    // ________________________________________________________________________
    virtual ~DirichletProposer(){}
    // ________________________________________________________________________
    virtual std::vector<double>                    proposeValues(
        std::vector<double>oldValues,
        double             parameter,
        Randomness&        rand,
        log_double&        hastings);
    // ________________________________________________________________________
    virtual AbstractProposer*                      clone() const
    {return new DirichletProposer(*this);}

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    bool minMaxIsRelative;
    RateHelper _rateHelper;
};


#endif
