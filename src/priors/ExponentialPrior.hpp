#ifndef _EXPONENTIAL_PRIOR
#define _EXPONENTIAL_PRIOR

#include "AbstractPrior.hpp"
#include "Density.hpp"

///////////////////////////////////////////////////////////////////////////////
//                             EXPONENTIAL PRIOR                             //
///////////////////////////////////////////////////////////////////////////////
class ExponentialPrior : public AbstractPrior
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ExponentialPrior(
        double lambda);
    // ________________________________________________________________________
    virtual bool                                needsIntegration() const
    {return true; }
    // ________________________________________________________________________
    virtual log_double                          getLogProb(
        const ParameterContent& content) const;
    // ________________________________________________________________________
    virtual log_double                          getUpdatedValue(
        double                   oldRawVal,
        double                   newRawVal,
        const AbstractParameter* param) const;
    // ________________________________________________________________________
    virtual ParameterContent                    drawFromPrior(
        Randomness&rand)  const;
    // ________________________________________________________________________
    virtual void                                print(
        std::ostream& out) const;
    // ________________________________________________________________________
    virtual double                              getLamda()  const
    {return _lambda; }
    // ________________________________________________________________________
    virtual ParameterContent                    getInitialValue() const;
    // ________________________________________________________________________
    virtual log_double                          accountForMeanSubstChange(
        TreeAln&                 traln,
        const AbstractParameter* param,
        double                   myOld,
        double                   myNew)  const;
    // ________________________________________________________________________
    virtual double                              getFirstDerivative(
        const AbstractParameter& param) const;
    // ________________________________________________________________________
    virtual AbstractPrior*                      clone() const
    {return new ExponentialPrior(*this); }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double _lambda;
};


#endif
