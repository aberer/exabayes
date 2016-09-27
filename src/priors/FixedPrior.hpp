#ifndef _FIXE_PRIOR
#define  _FIXE_PRIOR

#include "AbstractPrior.hpp"


///////////////////////////////////////////////////////////////////////////////
//                                FIXED PRIOR                                //
///////////////////////////////////////////////////////////////////////////////
class FixedPrior : public AbstractPrior
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    FixedPrior(
        std::vector<double>fixedValues);
    // ________________________________________________________________________
    virtual bool                                needsIntegration() const
    {return false;    }
    // ________________________________________________________________________
    virtual log_double                          getLogProb(
        const ParameterContent&content)  const;
    // ________________________________________________________________________
    virtual ParameterContent                    drawFromPrior(
        Randomness&rand)  const {assert(0); return ParameterContent{}; }
    // ________________________________________________________________________
    virtual void                                print(
        std::ostream&out) const;
    // ________________________________________________________________________
    virtual ParameterContent                    getInitialValue() const;
    // ________________________________________________________________________
    virtual AbstractPrior*                      clone() const
    {return new  FixedPrior(*this); }
    // ________________________________________________________________________
    virtual log_double                          accountForMeanSubstChange(
        TreeAln&                 traln,
        const AbstractParameter* param,
        double                   myOld,
        double                   myNew) const;
    // ________________________________________________________________________
    double                                      getFirstDerivative(
        const AbstractParameter& param) const {assert(0); return 0; }
    // ________________________________________________________________________
    virtual log_double                          getUpdatedValue(
        double                   oldRawVal,
        double                   newRawVal,
        const AbstractParameter* param) const
    {assert(0); return log_double(); }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<double> _fixedValues;
};


#endif
