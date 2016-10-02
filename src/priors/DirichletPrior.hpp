
#ifndef _DIRICHLET_PRIOR
#define _DIRICHLET_PRIOR

#include "AbstractPrior.hpp"

///////////////////////////////////////////////////////////////////////////////
//                              DIRICHLET PRIOR                              //
///////////////////////////////////////////////////////////////////////////////
class DirichletPrior : public AbstractPrior
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    DirichletPrior(
        std::vector<double>a)
        : alphas(a)
    {}
    // ________________________________________________________________________
    virtual ParameterContent                    drawFromPrior(
        Randomness&rand)  const {assert(0); return ParameterContent{};};
    // ________________________________________________________________________
    virtual log_double                          getLogProb(
        const ParameterContent& content) const;
    // ________________________________________________________________________
    virtual void                                print(
        std::ostream& out) const;
    // ________________________________________________________________________
    virtual ParameterContent                    getInitialValue() const;
    // ________________________________________________________________________
    virtual bool                                needsIntegration() const
    {return true;}
    // ________________________________________________________________________
    virtual log_double                          accountForMeanSubstChange(
        TreeAln&                 traln,
        const AbstractParameter* param,
        double                   myOld,
        double                   myNew) const
    {
        assert(0);
        return log_double::fromAbs(0);
    }
    // ________________________________________________________________________
    virtual AbstractPrior*                      clone() const
    {
        return new  DirichletPrior(*this);
    }
    // ________________________________________________________________________
    double                                      getFirstDerivative(
        const AbstractParameter& param) const
    {assert(0); return 0;}                       //
    // ________________________________________________________________________
    virtual log_double                          getUpdatedValue(
        double                   oldRawVal,
        double                   newRawVal,
        const AbstractParameter& param) const
    {assert(0); return log_double();}

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<double>alphas;
};


#endif
