#ifndef _ABSTRACT_PARAMETER
#define _ABSTRACT_PARAMETER

#include "TreeAln.hpp"
#include "ParameterContent.hpp"
#include "AbstractPrior.hpp"
#include "ParamAttribute.hpp"

#include "Serializable.hpp"

enum class Category;

///////////////////////////////////////////////////////////////////////////////
//                             ABSTRACT PARAMETER                            //
///////////////////////////////////////////////////////////////////////////////
class AbstractParameter : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    // generally, there is
    // nothing to do ...
    virtual auto
                                                                                                             deserialize(
        std::istream& in)
        ->void
    {}
    // ________________________________________________________________________
    virtual auto
    serialize(
        std::ostream&out) const
        ->void
    {}
    // ________________________________________________________________________
    AbstractParameter(
        Category        cat,
        nat             id,
        nat             idOfMyKind,
        std::vector<nat>partitions,
        nat             paramPrio);
    // ________________________________________________________________________
    AbstractParameter(
        const AbstractParameter& rhs);
    // ________________________________________________________________________
    virtual auto
                                                                                                             getPriorValue(
        const TreeAln& traln) const
        ->log_double  = 0;
    // ________________________________________________________________________
    // TODO we should just cast to BranchLengthParameter and make this specific
    // to branch length parameter ...
    virtual double
                                                                                                             getMeanSubstitutionRate()
    const
    {
        assert(0);
        return 0;
    }
    // ________________________________________________________________________
    virtual void
                                                                                                             updateMeanSubstRate(
        const TreeAln& traln){assert(0); }
    // ________________________________________________________________________
    /**
     *  @brief applies the parameter content to the tree
     */
    virtual void
                                                                                                             applyParameter(
        TreeAln&               traln,
        const ParameterContent&content) = 0;
    // ________________________________________________________________________
    virtual void
                                                                                                             applyParameterRaw(
        TreeAln&                traln,
        const ParameterContent& content) const {}
    // ________________________________________________________________________
    /**
     *  @brief extracts the parameter
     */
    virtual ParameterContent
                                                                                                             extractParameter(
        const TreeAln&traln)  const  = 0;
    // ________________________________________________________________________
    virtual ParameterContent
                                                                                                             extractParameterRaw(
        const TreeAln& traln) const {return ParameterContent{}; }
    // ________________________________________________________________________
    /**
     *  @brief print a sample for this parameter
     */
    virtual void
                                                                                                             printSample(
        std::ostream& fileHandle,
        const TreeAln&traln) const = 0;
    // ________________________________________________________________________
    /**
     *  @brief print the names of all components of this parameter (e.g., the
     * meaning of the various rates )
     */
    virtual void
                                                                                                             printAllComponentNames(
        std::ostream& fileHandle,
        const TreeAln&traln) const  = 0;
    // ________________________________________________________________________
    /**
     *  @brief adds a partition to the parameter (during setup)
     */
    void
                                                                                                             addPartition(
        nat id){_partitions.push_back(id); }
    // ________________________________________________________________________
    /**
     *  @brief sets the prior for this parameter
     */
    virtual void
    setPrior(
        const std::unique_ptr<AbstractPrior>&prior)
    {_prior = std::unique_ptr<AbstractPrior>(prior->clone()); }
    // ________________________________________________________________________
    nat
                                                                                                             getIdOfMyKind()
    const
    {return _idOfMyKind; }
    // ________________________________________________________________________
    /**
     *  @brief veriffies that content is compatible to this parameter (e.g.,
     * not too many rates).
     *
     *  This is a crude method merely for initialization (user input
     * validation)
     */
    virtual void
                                                                                                             verifyContent(
        const TreeAln&         traln,
        const ParameterContent&content) const  =  0;

    // ________________________________________________________________________
    Category
                                                                                                             getCategory()
    const
    {return _cat; }
    // ________________________________________________________________________
    nat
    getId()
    const
    {return _id; }
    // ________________________________________________________________________
    std::vector<nat>
                                                                                                             getPartitions()
    const
    {return _partitions; }
    // ________________________________________________________________________
    auto
    getPrior()
    const
        ->AbstractPrior
    * {return _prior.get(); }
    // ________________________________________________________________________
    bool
                                                                                                             isPrintToParamFile()
    const
    {return _printToParamFile; }
    // ________________________________________________________________________
    virtual auto
    printShort(
        std::ostream& out) const
        ->std::ostream &;
    // ________________________________________________________________________
    friend auto
    operator<<(
        std::ostream&            out,
        const AbstractParameter* rhs)
        ->std::ostream &;
    // ________________________________________________________________________
    virtual AbstractParameter*
    clone()
    const = 0;
    // ________________________________________________________________________
    virtual bool
                                                                                                             priorIsFitting(
        const AbstractPrior&prior,
        const TreeAln&      traln) const;
    // ________________________________________________________________________
    virtual bool
                                                                                                             fitsToPartition(
        Partition& p) const
    {return true; }
    // ________________________________________________________________________
    virtual void
                                                                                                             checkSanityPartitionsAndPrior(
        const TreeAln&traln) const;
    // ________________________________________________________________________
    nat
                                                                                                             getParamPriority()
    const
    {return _paramPriority; }
    // ________________________________________________________________________
    virtual ParamAttribute
                                                                                                             getAttributes()
    const
    {
        assert(0);
        return ParamAttribute();
    }

    // ________________________________________________________________________
    // currently only used for tuning of branch length parameter
    virtual void
                                                                                                             setAttributes(
        ParamAttribute attr){assert(0); }                          // same here


    ///////////////////////////////////////////////////////////////////////////
    //                          PROTECTED INTERFACE                          //
    ///////////////////////////////////////////////////////////////////////////
protected:
    void
                                                                                                             checkSanityPartitionsAndPrior_FreqRevMat(
        const TreeAln&traln) const;

    ///////////////////////////////////////////////////////////////////////////
    //                             PROTECTED DATA                            //
    ///////////////////////////////////////////////////////////////////////////
protected:
    // ATTRIBUTES
    nat                            _id;
    nat                            _idOfMyKind;
    Category                       _cat;
    std::unique_ptr<AbstractPrior> _prior;
    bool                           _printToParamFile;
    std::vector<nat>               _partitions;
    nat                            _paramPriority;
};


#endif
