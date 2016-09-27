#ifndef REV_MAT_PARAMETER
#define REV_MAT_PARAMETER

#include "RateHelper.hpp"
#include "Category.hpp"
#include "AbstractParameter.hpp"

///////////////////////////////////////////////////////////////////////////////
//                             REV MAT PARAMETER                             //
///////////////////////////////////////////////////////////////////////////////
class RevMatParameter : public AbstractParameter
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    RevMatParameter(
        nat             id,
        nat             idOfMyKind,
        std::vector<nat>partitions)
        : AbstractParameter(Category::SUBSTITUTION_RATES, id, idOfMyKind,
            partitions, 1)
    {}
    // ________________________________________________________________________
    RevMatParameter(
        const RevMatParameter&rhs)
        : AbstractParameter(rhs)
    {}
    // ________________________________________________________________________
    virtual void                                 applyParameter(
        TreeAln&               traln,
        const ParameterContent&content);
    // ________________________________________________________________________
    virtual void                                 applyParameterRaw(
        TreeAln&                traln,
        const ParameterContent& content) const;
    // ________________________________________________________________________
    virtual ParameterContent                     extractParameter(
        const TreeAln&traln)  const;
    // ________________________________________________________________________
    virtual ParameterContent                     extractParameterRaw(
        const TreeAln& traln) const;
    // ________________________________________________________________________
    virtual AbstractParameter*                   clone() const
    {return new RevMatParameter(*this); }
    // ________________________________________________________________________
    virtual void                                 printSample(
        std::ostream& fileHandle,
        const TreeAln&traln) const;
    // ________________________________________________________________________
    virtual void                                 printAllComponentNames(
        std::ostream& fileHandle,
        const TreeAln&traln) const;
    // ________________________________________________________________________
    virtual void                                 verifyContent(
        const TreeAln&         traln,
        const ParameterContent&content) const;
    // ________________________________________________________________________
    virtual void                                 checkSanityPartitionsAndPrior(
        const TreeAln&traln) const;
    // ________________________________________________________________________
    virtual bool                                 priorIsFitting(
        const AbstractPrior&prior,
        const TreeAln&      traln) const
    {
        auto  content = prior.getInitialValue();
        auto& partition = traln.getPartition(_partitions.at(0));
        return content.values.size()  ==
               RateHelper::numStateToNumInTriangleMatrix(
            partition.getStates());
    }
    // ________________________________________________________________________
    virtual log_double                           getPriorValue(
        const TreeAln& traln) const
    {
        assert(0);
        return log_double::fromAbs(1);
    }
};


#endif
