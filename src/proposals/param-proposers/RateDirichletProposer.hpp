#ifndef _PARTIAL_DIRICHLET_HPP
#define _PARTIAL_DIRICHLET_HPP

#include "AbstractProposer.hpp"
#include "RateHelper.hpp"

///////////////////////////////////////////////////////////////////////////////
//                          RATE DIRICHLET PROPOSER                          //
///////////////////////////////////////////////////////////////////////////////
class RateDirichletProposer :  public AbstractProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    RateDirichletProposer(
        double minValI,
        double maxValI);
    // ________________________________________________________________________
    RateDirichletProposer(
        const RateDirichletProposer& rhs) = default;
    // ________________________________________________________________________
    virtual std::vector<double>                    proposeValues(
        std::vector<double>oldValues,
        double             parameter,
        Randomness&        rand,
        log_double&        hastings);
    // ________________________________________________________________________
    virtual AbstractProposer*                      clone() const;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    auto                                           correctValues(
        nat                whichFreq,
        nat                numFreq,
        std::vector<double>newRates,
        std::vector<double>allNewValues)
        ->std::tuple<std::vector<double>, std::vector<double> >;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    RateHelper _rateHelper;
};


#endif
