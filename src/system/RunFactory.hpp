/**
 * @file RunFactory
 *
 * @brief an essential class that uses all available information to
 * create the exact setup of the run in terms of priors, proposals and
 * so on
 *
 *
 */


#ifndef _PROPOSALFACTORY
#define _PROPOSALFACTORY

#include "ProposalSet.hpp"
#include "BlockProposalConfig.hpp"
#include "BlockPrior.hpp"
#include "BlockParams.hpp"
#include "AbstractProposal.hpp"
#include "TreeAln.hpp"
#include "GlobalVariables.hpp"
#include "AbstractParameter.hpp"
#include "ParameterList.hpp"

#include <vector>


///////////////////////////////////////////////////////////////////////////////
//                                RUN FACTORY                                //
///////////////////////////////////////////////////////////////////////////////
class RunFactory
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    /**
     *  @brief configures the proposals
     */
    auto                    produceProposals(
        const BlockProposalConfig&propConfig,
        const BlockPrior&         priorInfo,
        ParameterList&            params,
        const TreeAln&            traln,
        bool                      componentWiseMH,
        ParallelSetup&            pl)
        ->std::tuple<std::vector<std::unique_ptr<AbstractProposal> >,
                     std::vector<ProposalSet> >;
    // ________________________________________________________________________
    /**
     *  @brief get a copy of the random variables to be integrated
     */
    void                    addStandardParameters(
        std::vector<std::unique_ptr<AbstractParameter> >&vars,
        const TreeAln&                                   traln) const;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                    addStandardPrior(
        AbstractParameter* var,
        const TreeAln&     traln);
    // ________________________________________________________________________
    void                    addPriorsToParameters(
        const TreeAln&   traln,
        const BlockPrior&priorInfo,
        ParameterList&   variables);
    // ________________________________________________________________________
    /**
     *  @brief adds secondary parameters to proposals, if necessary (currently
     * only branch lengths)
     */
    void                    addSecondaryParameters(
        AbstractProposal* proposal,
        ParameterList&    allParameters,
        nat               numTaxa);
};


#endif
