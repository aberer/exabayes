#ifndef _COMPLEX_TUNER_HPP
#define _COMPLEX_TUNER_HPP

#include "SuccessCounter.hpp"
#include "Serializable.hpp"


/**
 * This is a style of tuning that is necessary for tuning the branch
 * length distribution proposal.
 *
 */


///////////////////////////////////////////////////////////////////////////////
//                               COMPLEX TUNER                               //
///////////////////////////////////////////////////////////////////////////////
class ComplexTuner : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    //
    // ________________________________________________________________________
    // implements Serializable
    //
    virtual void                    deserialize(
        std::istream&in);
    // ________________________________________________________________________
    // implements Serializable
    //
    virtual void                    serialize(
        std::ostream&out) const;

    // ________________________________________________________________________
    ComplexTuner(
        double parameter,
        double minBound,
        double maxBound,
        double maxDelta,
        bool   logScale)
        : _parameter{parameter}
        , _sctr{}
        , _tuneUp{true}
        , _prevSuccess{0.}
        , _minBound{minBound}
        , _maxBound{maxBound}
        , _maxDelta{maxDelta}
        , _logScale{logScale}
    {}
    // ________________________________________________________________________
    void                            accept(){_sctr.accept(); }
    // ________________________________________________________________________
    void                            reject(){_sctr.reject(); }
    // ________________________________________________________________________
    double                          getRecentlySeen() const
    {return _sctr.getRecentlySeen(); }
    // ________________________________________________________________________
    double                          tuneParameter(
        int    batch,
        double parameter,
        bool   increase) const;
    // ________________________________________________________________________
    void                            tune();
    // ________________________________________________________________________
    void                            setParameter(
        double p){_parameter = p; }
    // ________________________________________________________________________
    double                          getParameter() const
    {return _parameter; }
    // ________________________________________________________________________
    double                          getRatio() const
    {return _sctr.getRatioInLastInterval(); }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double         _parameter;
    SuccessCounter _sctr;
    bool           _tuneUp;
    double         _prevSuccess;

    double         _minBound;
    double         _maxBound;

    double         _maxDelta; // upper bound by how much we tune in one step

    bool           _logScale;
};


#endif
