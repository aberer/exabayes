#ifndef PROT_MODEL_PARAMETER
#define PROT_MODEL_PARAMETER

#include "Category.hpp"
#include "AbstractParameter.hpp"

///////////////////////////////////////////////////////////////////////////////
//                            PROT MODEL PARAMETER                           //
///////////////////////////////////////////////////////////////////////////////
class ProtModelParameter : public AbstractParameter
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ProtModelParameter(
        nat             id,
        nat             idOfMyKind,
        std::vector<nat>partitions)
        : AbstractParameter(Category::AA_MODEL, id, idOfMyKind, partitions, 1)
    {}
    // ________________________________________________________________________
    virtual void                                 applyParameter(
        TreeAln&               traln,
        const ParameterContent&content);
    // ________________________________________________________________________
    virtual ParameterContent                     extractParameter(
        const TreeAln&traln)  const;
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
    virtual AbstractParameter*                   clone() const
    {return new ProtModelParameter(*this); }
    // ________________________________________________________________________
    virtual void                                 checkSanityPartitionsAndPrior(
        const TreeAln&traln) const;
    // ________________________________________________________________________
    virtual bool                                 fitsToPartition(
        Partition& p) const;
    // ________________________________________________________________________
    virtual log_double                           getPriorValue(
        const TreeAln& traln) const
    {
        assert(0);
        return log_double::fromAbs(1);
    }
};


#endif
