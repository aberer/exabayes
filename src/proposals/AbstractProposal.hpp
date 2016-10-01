#ifndef _ABSTRACTPROPOSAL_H
#define _ABSTRACTPROPOSAL_H

#include "Category.hpp"
#include "SuccessCounter.hpp"
#include "Randomness.hpp"
#include "PriorBelief.hpp"
#include "GlobalVariables.hpp"
#include "LikelihoodEvaluator.hpp"
#include "TreeRandomizer.hpp"
#include "Serializable.hpp"
#include "DistributionProposer.hpp"
#include "GammaProposer.hpp"
#include "ParameterList.hpp"

#include <string>
#include <vector>

class OptimizedParameter;

///////////////////////////////////////////////////////////////////////////////
//                             ABSTRACT PROPOSAL                             //
///////////////////////////////////////////////////////////////////////////////
class AbstractProposal : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    using UPtr = std::unique_ptr<AbstractProposal>;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    // INHERITED METHODS
    /**
     * @brief determines the proposal, applies it to the tree / model, updates
     * prior and hastings ratio
     */
    virtual void
                        applyToState(
        TreeAln&             traln,
        PriorBelief&         prior,
        log_double&          hastings,
        Randomness&          rand,
        LikelihoodEvaluator& eval) = 0;
    // ________________________________________________________________________
    /**
     * @brief evaluates the proposal
     * @todo remove the prior, we should not need it here
     */
    virtual void
                        evaluateProposal(
        LikelihoodEvaluator&evaluator,
        TreeAln&            traln,
        const BranchPlain&  branchSuggestion) = 0;
    // ________________________________________________________________________
    /**
     *  @brief resets the tree to its previous state; corrects the prior, if
     * necessary (@todo is this the case?)
     */
    virtual void
                        resetState(
        TreeAln&traln) = 0;
    // ________________________________________________________________________
    /**
     *  @brief tunes proposal parameters, if available
     */
    virtual void
                        autotune() =
        0;
    // ________________________________________________________________________
    virtual auto
                        clone() const
        ->AbstractProposal *  = 0;

    // ________________________________________________________________________
    virtual BranchPlain
                        determinePrimeBranch(
        const TreeAln&traln,
        Randomness&   rand) const = 0;
    // ________________________________________________________________________
    /**
     *  @brief gets nodes that are invalid by executed the proposal
     */
    virtual std::vector<nat>
                        getInvalidatedNodes(
        const TreeAln&traln) const = 0;
    // ________________________________________________________________________
    /**
     *  @brief inform proposal about acceptance
     */
    virtual void
                        accept()
    {_sctr.accept();}
    // ________________________________________________________________________
    /**
     * @brief inform proposal about rejection
     */
    virtual void
                        reject()
    {_sctr.reject();}
    // ________________________________________________________________________
    virtual void
                        extractProposer(
        TreeAln&                  traln,
        const OptimizedParameter& param){}
    // ________________________________________________________________________
    /**
     *  @brief prepare for set execution (only relevant for branch length +
     * node slider)
     */
    virtual std::pair<BranchPlain,
                      BranchPlain>
                        prepareForSetExecution(
        TreeAln&   traln,
        Randomness&rand)  = 0;
    // ________________________________________________________________________
    // we need to implement these
    virtual void
                        serialize(
        std::ostream&out)  const;
    // ________________________________________________________________________
    virtual void
                        deserialize(
        std::istream&in);
    // ________________________________________________________________________
    /**
     *  @brief writes proposal specific (tuned) parameters
     */
    virtual void
                        writeToCheckpointCore(
        std::ostream&out) const  = 0;
    // ________________________________________________________________________
    /**
     *  @brief reads proposal specific (tuned) parameters
     */
    virtual void
                        readFromCheckpointCore(
        std::istream&in) = 0;
    // ________________________________________________________________________
    virtual void
                        prepareForSetEvaluation(
        TreeAln&             traln,
        LikelihoodEvaluator& eval) const {}
    // ________________________________________________________________________
    virtual void
                        printParams(
        std::ostream&out)  const {}
    // ________________________________________________________________________
    AbstractProposal(
        Category    cat,
        std::string _name,
        double      weight,
        double      minTuning,
        double      maxTuning,
        bool        needsFullTraversal);
    // ________________________________________________________________________
    virtual ~AbstractProposal(){}
    // ________________________________________________________________________
    bool
                        isUsingOptimizedBranches()
    const {return _usingOptimizedBranches;}
    // ________________________________________________________________________
    std::array<bool,
               3>
                        getBranchProposalMode()
    const;
    // ________________________________________________________________________
    /**
     *  @brief gets the relative weight of this proposal
     */
    double
                        getRelativeWeight()
    const {return _relativeWeight;}
    // ________________________________________________________________________
    /**
     * @brief sets the relative weight of this proposal
     */
    void
                        setRelativeWeight(
        double tmp){_relativeWeight = tmp;}
    // ________________________________________________________________________
    /**
     *  @brief gets the category
     */
    Category
                        getCategory()
    const {return _category;}
    // ________________________________________________________________________
    /**
     *  @brief gets the name of the proposal
     */
    std::string
                        getName()
    const {return _name;}
    // ________________________________________________________________________
    /**
     *  @brief gets the success counter
     */
    const SuccessCounter&
                                            getSCtr()
    const {return _sctr;}
    // ________________________________________________________________________
    /**
     *  @brief gets the number of proposal invocation, since it has been tune
     * the last time
     */
    int
                                            getNumCallSinceTuning()
    const
    {return _sctr.getRecentlySeen();}
    // ________________________________________________________________________
    /**
     *  @brief add a parameter to be integrated over to the proposal
     */
    void
                                            addPrimaryParameter(
        nat id)
    {_primParamIds.push_back(id);}
    // ________________________________________________________________________
    /**
     *  @brief add a parameter  that is integrated over as a by-product of this
     * proposal
     */
    void
                                            addSecondaryParameter(
        nat id){_secParamIds.push_back(id);}
    // ________________________________________________________________________
    /**
     *  @brief indicates whether this proposal needs a full traversal
     */
    bool
                                            isNeedsFullTraversal()
    const
    {return _needsFullTraversal;}
    // ________________________________________________________________________
    std::vector<AbstractParameter*>
                                            getBranchLengthsParameterView()
    const;
    // ________________________________________________________________________
    /**
     *  @brief update the hatsings
     *  @param valToAdd is the absolute proposal ratio that shall be added
     */
    std::vector<nat>
                                            getAffectedPartitions()
    const;
    // ________________________________________________________________________
    std::ostream&
                                            printNamePartitions(
        std::ostream&out);
    // ________________________________________________________________________
    std::ostream&
                                            printShort(
        std::ostream&out)  const;
    // ________________________________________________________________________
    std::vector<AbstractParameter*>
                                            getPrimaryParameterView()
    const;
    // ________________________________________________________________________
    std::vector<AbstractParameter*>
                                            getSecondaryParameterView()
    const;
    // ________________________________________________________________________
    void
                                            setPreparedBranch(
        BranchPlain b)
    {_preparedBranch = b;}
    // ________________________________________________________________________
    void
                                            setOtherPreparedBranch(
        BranchPlain b)
    {_preparedOtherBranch = b;}
    // ________________________________________________________________________
    void
                                            setInSetExecution(
        bool exec)
    {_inSetExecution = exec;}
    // ________________________________________________________________________
    void
                                            setId(
        nat id)
    {_id = id;}
    // ________________________________________________________________________
    nat
                                            getId() const
    {return _id;}
    // ________________________________________________________________________
    friend std::ostream&
                                            operator<<(
        std::ostream&           out,
        const AbstractProposal& rhs);
    // ________________________________________________________________________
    void
                                            setForProteinsOnly(
        bool val){_forProteinOnly = val;}
    // ________________________________________________________________________
    bool
                                            isForProteinOnly()
    const {return _forProteinOnly;}
    // ________________________________________________________________________
    double
                                            tuneParameter(
        int    batch,
        double accRatio,
        double parameter,
        bool   inverse);
    // ________________________________________________________________________
    nat
                                            getNumTaxNeeded()
    const {return _numTaxNeeded;}

    // ________________________________________________________________________
    /*
     * this method should not be....rather unsafe. Ideally, we'd like to
     * have a factory that produces proposals and parameters alike...
     */
    void
                                            setParams(
        ParameterList*params){_allParams = params;}

    ///////////////////////////////////////////////////////////////////////////
    //                             PROTECTED DATA                            //
    ///////////////////////////////////////////////////////////////////////////
protected:
    // ATTRIBUTES
    std::string _name;
    SuccessCounter _sctr;
    Category _category;

    std::vector<nat>_primParamIds;  // it is the  primary purpose of this
                                    // proposal to integrate over these
                                    // parameters (in most cases only 1)
    std::vector<nat>_secParamIds;  // as a by-product also these random
                                   // variables are changed

    double _relativeWeight;
    bool _needsFullTraversal;
    bool _inSetExecution;

    // meh
    BranchPlain _preparedBranch;
    BranchPlain _preparedOtherBranch;

    nat _id;

    bool _forProteinOnly;

    double _minTuning;
    double _maxTuning;

    nat _numTaxNeeded;
    bool _usingOptimizedBranches;

    ParameterList*   _allParams; // call back to chain->params ; not owned here
};


#endif

