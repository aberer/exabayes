#ifndef _CHECKPOINTABLE_H
#define _CHECKPOINTABLE_H

#include "common.h"

#include <iostream>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <limits>

///////////////////////////////////////////////////////////////////////////////
//                                SERIALIZABLE                               //
///////////////////////////////////////////////////////////////////////////////
class Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    virtual void                    deserialize(
        std::istream&in)  = 0;
    // ________________________________________________________________________
    virtual void                    serialize(
        std::ostream&out) const = 0;
    // ________________________________________________________________________
    Serializable()
        : DELIM('&')
        , checkpointIsBinary(true)
    {}
    // ________________________________________________________________________
    virtual ~Serializable(){}
    // ________________________________________________________________________
    void                            getOfstream(
        std::string   name,
        std::ofstream&result);
    // ________________________________________________________________________
    void                            getIfstream(
        std::string   name,
        std::ifstream&result);
    // ________________________________________________________________________
    template<typename T>
    void                            cWrite(
        std::ostream&out,
        const T&     toWrite) const;
    // ________________________________________________________________________
    template<typename T>
    T                               cRead(
        std::istream&in);
    // ________________________________________________________________________
    std::string                     readString(
        std::istream&in);
    // ________________________________________________________________________
    void                            writeString(
        std::ostream&out,
        std::string  toWrite) const;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                            readDelimiter(
        std::istream&in);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    char DELIM;
    bool checkpointIsBinary;
};


#endif
