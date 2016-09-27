#ifndef _MY_EXTENSIONS_HPP
#define _MY_EXTENSIONS_HPP

#include <vector>
#include <type_traits>

#include <memory>
#include <iostream>

#include "log_double.hpp"

#include "common.h"

#define PRINT(x) (tout << "#x=x"  << std::endl)


#if !defined(HAVE_CXX14)

namespace std {
// START definitions
// ____________________________________________________________________________
template<typename T, typename ... Args>
std::unique_ptr<T>                    make_unique(
    Args&& ... args);
// END

// still missing in c++11.
// ____________________________________________________________________________
template<typename T, typename ... Args>
std::unique_ptr<T>                    make_unique(
    Args&& ... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args) ...));
}
}
#endif


// a bit more number formatting

///////////////////////////////////////////////////////////////////////////////
//                            THOUSANDS SEPARATOR                            //
///////////////////////////////////////////////////////////////////////////////
template<typename T>
class ThousandsSeparator : public std::numpunct<T>
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ThousandsSeparator(
        T Separator)
        : m_Separator(Separator)
    {}

    ///////////////////////////////////////////////////////////////////////////
    //                          PROTECTED INTERFACE                          //
    ///////////////////////////////////////////////////////////////////////////
protected:
    // ________________________________________________________________________
    std::string                       do_grouping() const
    {return "\03";}

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    T m_Separator;
};


// ____________________________________________________________________________
void                                  formatRange(
    std::ostream&          out,
    const std::vector<nat>&values);


///////////////////////////////////////////////////////////////////////////////
//                               SHARED POD PTR                              //
///////////////////////////////////////////////////////////////////////////////
template<typename T>
struct shared_pod_ptr
{
    // ________________________________________________________________________
    shared_pod_ptr(
        T*arg = nullptr)
        : _impl{arg, [](T* bla){free(bla);}}
    {}
    // ________________________________________________________________________
    T*                                get(){return _impl.get();}
private:
    std::shared_ptr<T>               _impl;
};


// ____________________________________________________________________________
template<typename T, int ALIGN>
T*                                    aligned_malloc(
    size_t size);


#include "AlignedAllocator.hpp"

// ____________________________________________________________________________
template<typename T>
struct aligned_vector
{
    typedef std::vector<T, AlignedAllocator<T, EXA_ALIGN> >type;
};


// ____________________________________________________________________________
// that's mostly for development. accessing an environment variable is easier
// than
std::string                           getEnvironmentVariable(
    const std::string& key);


#endif
