#ifndef DISCRETE_MODEL_PRIOR
#define DISCRETE_MODEL_PRIOR

#include "ProtModel.hpp"
#include "AbstractPrior.hpp"

#include <unordered_map>


///////////////////////////////////////////////////////////////////////////////
//                            DISCRETE MODEL PRIOR                           //
///////////////////////////////////////////////////////////////////////////////
class DiscreteModelPrior : public AbstractPrior
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    DiscreteModelPrior(
        std::unordered_map<ProtModel, double>model);
    // ________________________________________________________________________
    // if we have only one model this is basically a fixed prior
    virtual bool                                needsIntegration() const
    {
        assert(_modelProbs.size() > 0);
        return _modelProbs.size() > 1;
    }
    // ________________________________________________________________________
    virtual ParameterContent                    getInitialValue() const;
    // ________________________________________________________________________
    virtual log_double                          accountForMeanSubstChange(
        TreeAln&                 traln,
        const AbstractParameter* param,
        double                   myOld,
        double                   myNew) const;
    // ________________________________________________________________________
    ParameterContent                            drawFromPrior(
        Randomness&rand)  const;
    // ________________________________________________________________________
    virtual log_double                          getLogProb(
        const ParameterContent& content) const;
    // ________________________________________________________________________
    virtual void                                print(
        std::ostream&out) const;
    // ________________________________________________________________________
    double                                      getFirstDerivative(
        const AbstractParameter& param) const {assert(0); return 0; }
    // ________________________________________________________________________
    virtual log_double                          getUpdatedValue(
        double                   oldRawVal,
        double                   newRawVal,
        const AbstractParameter* param) const
    {assert(0); return log_double(); }
    // ________________________________________________________________________
    virtual AbstractPrior*                      clone() const
    {return new  DiscreteModelPrior(*this); }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::unordered_map<ProtModel, double> _modelProbs;
};


#endif
