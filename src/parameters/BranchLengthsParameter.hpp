#ifndef _BRANCH_LENGTHS_PARAMETER
#define _BRANCH_LENGTHS_PARAMETER

#include "AbstractParameter.hpp"
#include "ComplexTuner.hpp"

#include "Category.hpp"

class BranchLengthsParameter : public AbstractParameter
{
public:
    // inherited from SERIALIZABLE
    // ________________________________________________________________________
    virtual void                                 deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                                 serialize(
        std::ostream&out) const;
public:
    // INHERITED METHODS
    // ________________________________________________________________________
    virtual void                                 applyParameter(
        TreeAln&               traln,
        const ParameterContent&content);
    // ________________________________________________________________________
    virtual ParameterContent                     extractParameter(
        const TreeAln&traln)  const;
    // ________________________________________________________________________
    virtual AbstractParameter*                   clone() const
    {
        return new BranchLengthsParameter(*this);
    }
    // ________________________________________________________________________
    virtual void                                 printSample(
        std::ostream& fileHandle,
        const TreeAln&traln) const {}
    // ________________________________________________________________________
    virtual void                                 printAllComponentNames(
        std::ostream& fileHandle,
        const TreeAln&traln) const {}
    // ________________________________________________________________________
    virtual void                                 verifyContent(
        const TreeAln&         traln,
        const ParameterContent&content) const;
    // ________________________________________________________________________
    virtual bool                                 priorIsFitting(
        const AbstractPrior&prior,
        const TreeAln&      traln) const;
    // ________________________________________________________________________
    virtual ParamAttribute                       getAttributes() const
    {
        return {_convTuner, _nonConvTuner};
    }
    // ________________________________________________________________________
    virtual void                                 setAttributes(
        ParamAttribute attr)
    {
        _convTuner = attr._convTuner;
        _nonConvTuner = attr._nonConvTuner;
    }
    // ________________________________________________________________________
    virtual double                               getMeanSubstitutionRate()
    const;
    // ________________________________________________________________________
    virtual void                                 updateMeanSubstRate(
        const TreeAln& traln);
    // ________________________________________________________________________
    virtual void                                 setMeanSubstitutionRate(
        double fac){_fracChange = fac; }
    // ________________________________________________________________________
    virtual log_double                           getPriorValue(
        const TreeAln& traln) const
    {
        assert(0);
        return log_double::fromAbs(1);
    }
public:
    // METHODS
    // ________________________________________________________________________
    BranchLengthsParameter(
        nat             id,
        nat             idOfMyKind,
        std::vector<nat>partitions);
private:
    // ATTRIBUTES
    ComplexTuner _convTuner; // convergence parameter for distribution
                             // proposals
    ComplexTuner _nonConvTuner; // non-convergence parameter for the
                                // distribution proposals
    double       _fracChange;
};


#endif
