/**
 *  @file ProposalSet.hpp
 *
 *  @brief represents a set of proposals
 *
 *  This class should be used to execute component-wise
 *  Metropolis-Hastings (a.k.a. component-in-sequence MH).
 *
 *  The main motivation to do so is efficiency.
 *
 *  The justification to do is e.g.,
 *
 * http://theclevermachine.wordpress.com/2012/11/04/mcmc-multivariate-distributions-block-wise-component-wise-updates/
 *
 *  and more specifically
 *  The Metropolitan-Hastings Algorithm and Extensions.S. Sawyer — Washington
 * University — Vs. August 22, 2010
 */


#ifndef PROPOSAL_SET
#define PROPOSAL_SET

#include "AbstractProposal.hpp"

#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//                                PROPOSAL SET                               //
///////////////////////////////////////////////////////////////////////////////
class ProposalSet : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ProposalSet(
        double                                         relWeight,
        std::vector<std::unique_ptr<AbstractProposal> >proposals);
    // ________________________________________________________________________
    ProposalSet(
        const ProposalSet&rhs);
    // ________________________________________________________________________
    ProposalSet&                                      operator=(
        const ProposalSet&rhs);
    // ________________________________________________________________________
    /**
     *  @brief prints the proposal set
     */
    void                                              printVerboseAbbreviated(
        std::ostream&out,
        double       sum) const;
    // ________________________________________________________________________
    /**
     *  @brief indicates whether for this proposal set a full tree
     *  traversal is necessary
     */
    bool                                              needsFullTraversal();
    // ________________________________________________________________________
    nat                                               numerateProposals(
        nat ctr);
    // ________________________________________________________________________
    void                                              setParameterListPtr(
        ParameterList* pPtr);
    // ________________________________________________________________________
    double                                            getRelativeWeight() const
    {return relativeWeight;}
    // ________________________________________________________________________
    auto                                              getProposalView() const
        ->std::vector<AbstractProposal*>;
    // ________________________________________________________________________
    virtual void                                      deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                                      serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    friend std::ostream&                              operator<<(
        std::ostream&     out,
        const ProposalSet&rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double relativeWeight;
    std::vector<std::unique_ptr<AbstractProposal> >proposals;
};


#endif
