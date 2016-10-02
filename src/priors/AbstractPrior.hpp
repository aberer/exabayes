#ifndef _PRIORS_H
#define _PRIORS_H

#include "Density.hpp"

#include "Randomness.hpp"
#include "TreeAln.hpp"

#include "ParameterContent.hpp"

#include <vector>
#include <limits>
#include <cassert>
#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//                               ABSTRACT PRIOR                              //
///////////////////////////////////////////////////////////////////////////////
class AbstractPrior
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    using UPtr = std::unique_ptr<AbstractPrior>;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    AbstractPrior()
        : _keepInitData{false}
    {}
    // ________________________________________________________________________
    virtual ~AbstractPrior(){}
    // ________________________________________________________________________
    /**
     *  @brief obtains a pre-defined initial value, depending on the
     *  prior.
     */
    virtual ParameterContent                    getInitialValue() const = 0;
    // ________________________________________________________________________
    virtual log_double                          accountForMeanSubstChange(
        TreeAln&                 traln,
        const AbstractParameter* param,
        double                   myOld,
        double                   myNew) const = 0;
    // ________________________________________________________________________
    virtual ParameterContent                    drawFromPrior(
        Randomness&rand)  const = 0;
    // ________________________________________________________________________
    virtual log_double                          getLogProb(
        const ParameterContent&content) const = 0;
    // ________________________________________________________________________
    virtual void                                print(
        std::ostream&out) const = 0;
    // ________________________________________________________________________
    virtual bool                                needsIntegration() const = 0;
    // ________________________________________________________________________
    virtual AbstractPrior*                      clone() const = 0;
    // ________________________________________________________________________
    virtual double                              getFirstDerivative(
        const AbstractParameter& param) const = 0;
    // ________________________________________________________________________
    // only for internal branch lengths; this is very ugly, however we need
    // this for maintaining numerical stability with > 300 partitions
    virtual auto                                getUpdatedValue(
        double                   oldRawVal,
        double                   newRawVal,
        const AbstractParameter& param) const
        ->log_double = 0;
    // ________________________________________________________________________
    friend auto                                 operator<<(
        std::ostream&  out,
        AbstractPrior* rhs)
        ->std::ostream
    &
    {
        rhs->print(out);
        return out;
    }
    // ________________________________________________________________________
    bool                                        isKeepInitData() const
    {
        return _keepInitData;
    }
    // ________________________________________________________________________
    void                                        setKeepInitData()
    {
        _keepInitData = true;
    }

    ///////////////////////////////////////////////////////////////////////////
    //                             PROTECTED DATA                            //
    ///////////////////////////////////////////////////////////////////////////
protected:
    bool _keepInitData;
};


#endif
