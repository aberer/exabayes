#ifndef _GAMMA_PROPOSER_HPP
#define _GAMMA_PROPOSER_HPP

#include "BranchLength.hpp"
#include "extensions.hpp"

class AbstractParameter;
class Randomness;
class TreeAln;

///////////////////////////////////////////////////////////////////////////////
//                               GAMMA PROPOSER                              //
///////////////////////////////////////////////////////////////////////////////
class GammaProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    GammaProposer(
        double nrOpt,
        double nrd1,
        double nrd2,
        double convParameter,
        double nonConvParameter);
    // ________________________________________________________________________
    BranchLength                                     proposeBranch(
        BranchPlain        b,
        TreeAln&           traln,
        AbstractParameter* param,
        Randomness&        rand) const;
    // ________________________________________________________________________
    log_double                                       getLogProbability(
        double val) const;
    // ________________________________________________________________________
    bool                                             isConverged() const
    {return _nrD2 < 0  && _nrD1 < 1; }
    // ________________________________________________________________________
    static std::string                               getName()
    {return std::string("Gamma"); }
    // ________________________________________________________________________
    friend std::ostream& operator                    <<(
        std::ostream&       out,
        const GammaProposer&rhs)
    {
        out << SHOW(rhs._nrOpt) << SHOW(rhs._nrD1) << SHOW(rhs._nrD2) << SHOW(
            rhs._alpha) << SHOW(rhs._beta);
        return out;
    }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double _nrOpt;      // branch length optimum
    double _nrD1;       // first derivative
    double _nrD2;       // second derivative

    double _alpha;
    double _beta;

    double _convParameter;
    double _nonConvParameter;
};


#endif
