#ifndef _LOG_DOUBLE
#define _LOG_DOUBLE

#include "GlobalVariables.hpp"

#include <cmath>
#include <ostream>
#include <cassert>

#include <limits>

///////////////////////////////////////////////////////////////////////////////
//                                 LOG DOUBLE                                //
///////////////////////////////////////////////////////////////////////////////
class log_double
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    log_double()
        : _val{0}
        , _error{0}
    {}
    // ________________________________________________________________________
    static log_double                      fromAbs(
        double val)
    {
        if (val <= 0)
        {
            tout << "error: tried to create a log-double of " <<  val
                 << std::endl;
            assert(val > 0);
        }

        auto result = log_double();
        result._val = log(val);
        return result;
    }
    // ________________________________________________________________________
    static log_double                      fromLog(
        double val)
    {
        auto result = log_double();
        result._val = val;
        return result;
    }
    // ________________________________________________________________________
    // could overload numeric_limits<log_double>::min()
    static log_double                      lowest()
    {
        auto result = log_double::fromAbs(1.);
        result._val = std::numeric_limits<double>::lowest();
        return result;
    }
    // ________________________________________________________________________
    bool                                   isNegativeInfinity() const
    {return isInfinity() && _val < 0.;}
    // ________________________________________________________________________
    bool                                   isInfinity() const
    {return std::isinf(_val);}
    // ________________________________________________________________________
    bool                                   isNaN() const
    {return std::isnan(_val);}

    // maybe allow this again, if we can be sure that log_double is used
    // carefully.

    // operator double()
    // {
    //   auto tmp = exp(_val);
    //   return tmp;
    // }
    // ________________________________________________________________________
    double                                 toAbs() const
    {return exp(_val);}
    // ________________________________________________________________________
    static log_double                      negativeInfinity()
    {
        auto result = log_double();
        result._val = 0;
        return result;
    }
    // ________________________________________________________________________
    friend bool                            operator<(
        const log_double& lhs,
        const log_double& rhs)
    {return lhs._val < rhs._val;}
    // ________________________________________________________________________
    friend void                            operator /=(
        log_double&      lhs,
        const log_double&rhs)
    {lhs = lhs / rhs;}
    // ________________________________________________________________________
    friend void                            operator*=(
        log_double&       lhs,
        const log_double& rhs)
    {lhs = lhs * rhs;}
    // ________________________________________________________________________
    friend log_double                      operator/(
        const log_double&lhs,
        const log_double&rhs);
    // ________________________________________________________________________
    friend log_double                      operator*(
        const log_double& lhs,
        const log_double& rhs);

    // ________________________________________________________________________
    // this is only meant for communication with PLL
    double                                 getRawLog() const
    {return _val;}
    // ________________________________________________________________________
    friend log_double                      exponentiate(
        log_double val,
        double     exponent);
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&     out,
        const log_double& rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double _val;
    double _error;
};


#endif
