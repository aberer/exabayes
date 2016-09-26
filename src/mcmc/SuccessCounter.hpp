#ifndef _SUCCESSCTR_H
#define _SUCCESSCTR_H

#include "Serializable.hpp"

#include <iostream>
#include <list>

#define SIZE_OF_LAST 100

// TODO  an success interval would be nice


///////////////////////////////////////////////////////////////////////////////
//                              SUCCESS COUNTER                              //
///////////////////////////////////////////////////////////////////////////////
class SuccessCounter : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    SuccessCounter();
    // ________________________________________________________________________
    SuccessCounter(
        const SuccessCounter& rhs) = default;
    // ________________________________________________________________________
    SuccessCounter(
        SuccessCounter&& rhs) = default;
    // ________________________________________________________________________
    SuccessCounter&                        operator=(
        const SuccessCounter& rhs)  = default;
    // ________________________________________________________________________
    SuccessCounter&                        operator=(
        SuccessCounter&& rhs)  = default;
    // ________________________________________________________________________
    void                                   accept();
    // ________________________________________________________________________
    void                                   reject();
    // ________________________________________________________________________
    int                                    getRecentlySeen() const
    {
        return localAcc + localRej;
    }
    // ________________________________________________________________________
    double                                 getRatioInLastInterval() const;
    // ________________________________________________________________________
    double                                 getRatioOverall() const;
    // ________________________________________________________________________
    void                                   nextBatch();
    // ________________________________________________________________________
    int                                    getBatch() const
    {return batch; }
    // ________________________________________________________________________
    nat                                    getTotalSeen() const
    {return globalAcc + globalRej; }
    // ________________________________________________________________________
    virtual void                           deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                           serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    SuccessCounter                         operator+(
        const SuccessCounter&rhs) const;
    // ________________________________________________________________________
    friend void                            swap(
        SuccessCounter&elemA,
        SuccessCounter&elemB);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                                   reset();
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&        rhs,
        const SuccessCounter&b);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    nat globalAcc;
    nat globalRej;
    nat localAcc;
    nat localRej;
    nat batch;
};


#endif
