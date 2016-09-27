#ifndef FREQ_PARAMETER
#define FREQ_PARAMETER

#include "TreeAln.hpp"
#include "AbstractParameter.hpp"
#include "Category.hpp"

///////////////////////////////////////////////////////////////////////////////
//                            FREQUENCY PARAMETER                            //
///////////////////////////////////////////////////////////////////////////////
class FrequencyParameter : public AbstractParameter
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    FrequencyParameter(
        nat             id,
        nat             idOfMyKind,
        std::vector<nat>partitions)
        : AbstractParameter(Category::FREQUENCIES, id, idOfMyKind, partitions,
            1)
    {}
    // ________________________________________________________________________
    FrequencyParameter(
        const FrequencyParameter&rhs)
        : AbstractParameter(rhs)
    {}
    // ________________________________________________________________________
    virtual void                                 applyParameter(
        TreeAln&               traln,
        const ParameterContent&content);
    // ________________________________________________________________________
    virtual ParameterContent                     extractParameter(
        const TreeAln&traln)  const;
    // ________________________________________________________________________
    virtual AbstractParameter*                   clone() const
    {return new FrequencyParameter(*this); }
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
        const TreeAln& traln) const;
    // ________________________________________________________________________
    virtual bool                                 priorIsFitting(
        const AbstractPrior&prior,
        const TreeAln&      traln) const
    {
        auto  content = prior.getInitialValue();
        auto& partition = traln.getPartition(_partitions.at(0));
        return content.values.size() == nat(partition.getStates());
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
