#ifndef _CHOSEN_BRANCH_INTEGRATOR_HPP
#define _CHOSEN_BRANCH_INTEGRATOR_HPP

#include "BranchLengthMultiplier.hpp"


///////////////////////////////////////////////////////////////////////////////
//                          CHOSEN BRANCH INTEGRATOR                         //
///////////////////////////////////////////////////////////////////////////////
class ChosenBranchIntegrator : public BranchLengthMultiplier
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ChosenBranchIntegrator(
        double multiplier)
        : BranchLengthMultiplier(multiplier)
        , _chosenBranches{}
    {}
    // ________________________________________________________________________
    virtual BranchPlain                    proposeBranch(
        const TreeAln&traln,
        Randomness&   rand) const
    {
        auto elem = rand.drawIntegerOpen(_chosenBranches.size());
        return _chosenBranches[elem];
    }
    // ________________________________________________________________________
    void                                   setChosenBranches(
        std::vector<BranchPlain>branches)
    {_chosenBranches = branches; }
    // ________________________________________________________________________
    AbstractProposal*                      clone() const
    {return new ChosenBranchIntegrator(*this); }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<BranchPlain> _chosenBranches;
};


#endif
