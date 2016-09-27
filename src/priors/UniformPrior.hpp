#ifndef _UNIFORM_PRIOR
#define _UNIFORM_PRIOR

#include "AbstractPrior.hpp"

///////////////////////////////////////////////////////////////////////////////
//                               UNIFORM PRIOR                               //
///////////////////////////////////////////////////////////////////////////////
class UniformPrior : public AbstractPrior
{
public:
    // ________________________________________________________________________
    UniformPrior(
        double minVal,
        double maxVal);
    // ________________________________________________________________________
    virtual log_double                          getLogProb(
        const ParameterContent& content)  const;
    // ________________________________________________________________________
    virtual bool                                needsIntegration() const
    {return true; }
    // ________________________________________________________________________
    virtual void                                print(
        std::ostream& out) const;
    // ________________________________________________________________________
    virtual ParameterContent                    getInitialValue() const;
    // ________________________________________________________________________
    virtual ParameterContent                    drawFromPrior(
        Randomness&rand)  const {assert(0); return ParameterContent{}; };
    // ________________________________________________________________________
    virtual log_double                          accountForMeanSubstChange(
        TreeAln&                 traln,
        const AbstractParameter* param,
        double                   myOld,
        double                   myNew) const;
    // ________________________________________________________________________
    virtual AbstractPrior*                      clone() const
    {return new  UniformPrior(*this); }
    // ________________________________________________________________________
    virtual log_double                          getUpdatedValue(
        double                   oldRawVal,
        double                   newRawVal,
        const AbstractParameter* param) const
    {assert(0); return log_double(); }
    // ________________________________________________________________________
    double                                      getMin() const
    {return minVal; }
    // ________________________________________________________________________
    double                                      getMax() const
    {return maxVal; }
    // ________________________________________________________________________
    double                                      getFirstDerivative(
        const AbstractParameter& param) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double minVal;
    double maxVal;
};


#endif
