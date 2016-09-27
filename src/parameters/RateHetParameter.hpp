#ifndef RATE_HET_PARAMETER
#define RATE_HET_PARAMETER

#include "AbstractParameter.hpp"
#include "Category.hpp"

///////////////////////////////////////////////////////////////////////////////
//                             RATE HET PARAMETER                            //
///////////////////////////////////////////////////////////////////////////////
class RateHetParameter : public AbstractParameter
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    RateHetParameter(
        nat             id,
        nat             idOfMyKind,
        std::vector<nat>partitions)
        : AbstractParameter(Category::RATE_HETEROGENEITY, id, idOfMyKind,
            partitions, 1)
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
    {return new RateHetParameter(*this); }
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
    virtual log_double                           getPriorValue(
        const TreeAln& traln) const
    {
        assert(0);
        return log_double::fromAbs(1);
    }
};


#endif
