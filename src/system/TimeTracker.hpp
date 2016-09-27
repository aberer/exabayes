#ifndef TIMETRACKER_H
#define TIMETRACKER_H

#include "Serializable.hpp"

#include <chrono>
#include <ratio>
#include <array>
#define CLOCK std::chrono


///////////////////////////////////////////////////////////////////////////////
//                                TIME TRACKER                               //
///////////////////////////////////////////////////////////////////////////////
class TimeTracker : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    virtual void                     deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                     serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    TimeTracker(
        size_t numCpu);
    // ________________________________________________________________________
    double                           getRecentElapsed()
    const;
    // ________________________________________________________________________
    void                             updateTime();
    // ________________________________________________________________________
    double                           getAccWallTime()
    const
    {return _accWallTime;}
    // ________________________________________________________________________
    double                           getAccCpuTime()
    const {return _accCpuTime;}
    // ________________________________________________________________________
    static double                    getDuration(
        CLOCK::system_clock::time_point tp);
    // ________________________________________________________________________
    static auto                      getTimePoint()
        ->CLOCK::system_clock::time_point;
    // ________________________________________________________________________
    static auto                      formatForWatch(
        double time)
        ->std::array<double, 3>;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    CLOCK::system_clock::time_point _timePoint;

    double _accWallTime;
    double _accCpuTime;

    size_t _numCpu;
};


#endif /* TIMETRACKER_H */
